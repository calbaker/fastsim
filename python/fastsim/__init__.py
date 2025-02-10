from importlib.metadata import version
from pathlib import Path
from typing import Any, List, Union, Dict 
from typing_extensions import Self
import numpy as np
import re
import inspect
import pandas as pd  # type: ignore[import-untyped]
import polars as pl

import fastsim
from .fastsim import *  # noqa: F403
from .fastsim import Cycle  # type: ignore[attr-defined]

DEFAULT_LOGGING_CONFIG = dict(
    format="%(asctime)s.%(msecs)03d | %(filename)s:%(lineno)s | %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def package_root() -> Path:
    """Returns the package root directory."""
    return Path(__file__).parent


def resources_root() -> Path:
    """
    Returns the resources root directory.
    """
    path = package_root() / "resources"
    return path


__version__ = version("fastsim")


def set_param_from_path(
    model: Any,
    path: str,
    value: Any,
) -> Any:
    """
    Set parameter `value` on `model` for `path` to parameter

    Example usage:
    todo
    """
    path_list = path.split(".")

    def _get_list(path_elem, container):
        list_match = re.match(r"([\w\d]+)\[(\d+)\]", path_elem)
        if list_match is not None:
            list_name = list_match.group(1)
            index = int(list_match.group(2))
            lst = container.__getattribute__(list_name).tolist()
            return lst, list_name, index
        else:
            return None, None, None

    containers = [model]
    lists = [(None, None, None)] * len(path_list)
    for i, path_elem in enumerate(path_list):
        container = containers[-1]

        list_attr, list_name, list_index = _get_list(path_elem, container)
        if list_attr is not None:
            attr = list_attr[list_index]
            # save for when we repack the containers
            lists[i] = (list_attr, list_name, list_index)
        else:
            attr = container.__getattribute__(path_elem)

        if i < len(path_list) - 1:
            containers.append(attr)

    prev_container = value

    # iterate through remaining containers, inner to outer
    for list_tuple, container, path_elem in zip(
        lists[-1::-1], containers[-1::-1], path_list[-1::-1]
    ):
        if list_tuple is not None:
            list_attr, list_name, list_index = list_tuple
            list_attr[list_index] = prev_container

            container.__setattr__(list_name, list_attr)
        else:
            container.__setattr__(f"__{path_elem}", prev_container)

        prev_container = container

    return model


def __array__(self):
    return np.array(self.tolist())


# creates a list of all python classes from rust structs that need variable_path_list() and
# history_path_list() added as methods
ACCEPTED_RUST_STRUCTS = [
    attr for attr in fastsim.__dir__() if not attr.startswith("__") and isinstance(getattr(fastsim, attr), type) and
    attr[0].isupper() and ("fastsim" in str(
        inspect.getmodule(getattr(fastsim, attr))))
]


def cyc_keys() -> List[str]:
    import json
    cyc = Cycle.from_resource("udds.csv")
    cyc_dict = json.loads(cyc.to_json())
    cyc_keys = [
        key for key, val in cyc_dict.items() if isinstance(val, list) and len(val) == cyc.len()
    ]

    return cyc_keys


CYC_KEYS = cyc_keys()

setattr(
    Pyo3VecWrapper,  # type: ignore[name-defined]  # noqa: F405
    "__array__",
    __array__
)  

# TODO connect to crate features
data_formats = [
    'yaml',
    'msg_pack',
    # 'toml',
    'json',
]


def to_pydict(self, data_fmt: str = "msg_pack", flatten: bool = False) -> Dict:
    """
    Returns self converted to pure python dictionary with no nested Rust objects
    # Arguments
    - `flatten`: if True, returns dict without any hierarchy
    - `data_fmt`: data format for intermediate conversion step
    """
    data_fmt = data_fmt.lower()
    assert data_fmt in data_formats, f"`data_fmt` must be one of {data_formats}"
    match data_fmt:
        case "msg_pack":
            import msgpack  # type: ignore[import-untyped]
            pydict = msgpack.loads(self.to_msg_pack())
        case "yaml":
            from yaml import load  # type: ignore[import-untyped]
            try:
                from yaml import CLoader as Loader
            except ImportError:
                from yaml import Loader
            pydict = load(self.to_yaml(), Loader=Loader)
        case "json":
            from json import loads
            pydict = loads(self.to_json())

    if not flatten:
        return pydict
    else:
        return next(iter(pd.json_normalize(pydict, sep=".").to_dict(orient='records')))


@classmethod  # type: ignore[misc]
def from_pydict(cls, pydict: Dict, data_fmt: str = "msg_pack", skip_init: bool = False) -> Self:  # type: ignore[misc]
    """
    Instantiates Self from pure python dictionary 
    # Arguments
    - `pydict`: dictionary to be converted to FASTSim object
    - `data_fmt`: data format for intermediate conversion step  
    - `skip_init`: passed to `SerdeAPI` methods to control whether initialization
      is skipped
    """
    data_fmt = data_fmt.lower()
    assert data_fmt in data_formats, f"`data_fmt` must be one of {data_formats}"
    match data_fmt.lower():
        case "yaml":
            import yaml
            obj = cls.from_yaml(yaml.dump(pydict), skip_init=skip_init)
        case "msg_pack":
            import msgpack
            obj = cls.from_msg_pack(
                msgpack.packb(pydict), skip_init=skip_init)
        case "json":
            from json import dumps
            obj = cls.from_json(dumps(pydict), skip_init=skip_init)

    return obj


def is_cyc_key(k):
    is_cyc_key = any(
        cyc_key for cyc_key in CYC_KEYS if cyc_key == k.split(".")[-1]) and "cyc" in k
    return is_cyc_key


def to_dataframe(self, pandas: bool = False, allow_partial: bool = False) -> Union[pd.DataFrame, pl.DataFrame]:
    """
    Returns time series results from fastsim object as a Polars or Pandas dataframe.

    # Arguments
    - `pandas`: returns pandas dataframe if True; otherwise, returns polars dataframe by default
    - `allow_partial`: returns dataframe of length equal to solved time steps if simulation fails early
    """
    obj_dict = self.to_pydict(flatten=True)
    history_dict: Dict[str, Any] = {}

    if len(history_dict) > 0:
        val0 = next(iter(history_dict.values()))
    else:
        val0 = None
    history_keys = ['.history.']
    for k, v in obj_dict.items():
        if is_cyc_key(k) or any(hk in k for hk in history_keys) or (val0 is not None and len(v) == len(val0)):
            history_dict[k] = v

    if allow_partial:
        cutoff = min([len(val) for val in history_dict.values()])

        if not pandas:
            df = pl.DataFrame({col: val[:cutoff]
                              for col, val in history_dict.items()})
        else:
            df = pd.DataFrame({col: val[:cutoff]
                              for col, val in history_dict.items()})
    else:
        if not pandas:
            try:
                df = pl.DataFrame(history_dict)
            except Exception as err:
                raise Exception(
                    f"{err}\nTry passing `allow_partial=True` to `to_dataframe`")
        else:
            try:
                df = pd.DataFrame(history_dict)
            except Exception as err:
                raise Exception(
                    f"{err}\nTry passing `allow_partial=True` to `to_dataframe`")
    return df


# adds variable_path_list() and history_path_list() as methods to all classes in
# ACCEPTED_RUST_STRUCTS
for item in ACCEPTED_RUST_STRUCTS:
    setattr(getattr(fastsim, item), "to_pydict", to_pydict)
    setattr(getattr(fastsim, item), "from_pydict", from_pydict)
    setattr(getattr(fastsim, item), "to_dataframe", to_dataframe)
