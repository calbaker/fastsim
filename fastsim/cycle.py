"""Module containing classes and methods for 
cycle data. For example usage, see ../README.md"""

# Import necessary python modules
from dataclasses import dataclass
import copy
import cmath
import os
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import re
import sys
from pathlib import Path
from copy import deepcopy
import types
from typing import Optional

# local modules
from . import parameters as params
from . import inspect_utils
from .rustext import RUST_AVAILABLE

if RUST_AVAILABLE:
    import fastsimrust as fsr

THIS_DIR = Path(__file__).parent
CYCLES_DIR = THIS_DIR / 'resources' / 'cycles'

# dicts for switching between legacy and new cycle keys
OLD_TO_NEW = {
    'cycSecs': 'time_s',
    'cycMps': 'mps',
    'cycGrade': 'grade',
    'cycRoadType': 'road_type',
    'name': 'name'
}
NEW_TO_OLD = {val: key for key, val in OLD_TO_NEW.items()}
STANDARD_CYCLE_KEYS = OLD_TO_NEW.values()


@dataclass
class Cycle(object):
    """Object for containing time, speed, road grade, and road charging vectors 
    for drive cycle.  Instantiate with the `from_file` or `from_dict` method.  
    """

    time_s: np.ndarray(1, dtype=float)
    mps: np.ndarray(1, dtype=float)
    grade: np.ndarray(1, dtype=float)
    road_type: np.ndarray(1, dtype=float)
    name: str

    @classmethod
    def from_file(cls, filename: str):
        """
        Load cycle from filename (str).
        Can be absolute or relative path.  If relative, looks in working dir first
        and then in `fastsim/resources/cycles`.  

        File must contain columns for:
        -- `cycSecs` or `time_s`
        -- `cycMps` or `mps`
        -- `cycGrade` or `grade` (optional)
        -- `cycRoadType` or `road_type` (optional)
        """
        filename = str(filename)

        if not filename.endswith('.csv'):
            filename += ".csv"
        if not Path(filename).exists() and (CYCLES_DIR / filename).exists():
            filename = CYCLES_DIR / filename
        elif Path(filename).exists():
            filename = Path(filename)
        else:
            raise ValueError("Invalid cycle filename.")

        cyc_df = pd.read_csv(filename)
        cyc_dict = cyc_df.to_dict(orient='list')
        cyc_dict = {key: np.array(val, dtype=float)
                    for key, val in cyc_dict.items()}
        cyc_dict['name'] = filename.stem

        return cls.from_dict(cyc_dict)

    @classmethod
    def from_dict(cls, cyc_dict: dict):
        """
        Load cycle from dict, which must contain keys for:
        -- `cycSecs` or `time_s`
        -- `cycMps` or `mps`
        -- `cycGrade` or `grade` (optional)
        -- `cycRoadType` or `road_type` (optional)
        """
        new_cyc_dict = {}
        for key, val in cyc_dict.items():
            # generate keys from legacy or current mapping
            new_cyc_dict[OLD_TO_NEW.get(key, key)] = val
        new_cyc_dict['name'] = cyc_dict.get('name', '')
        # set zeros if not provided
        new_cyc_dict['grade'] = new_cyc_dict.get(
            'grade', np.zeros(len(new_cyc_dict['time_s'])))
        new_cyc_dict['road_type'] = new_cyc_dict.get(
            'road_type', np.zeros(len(new_cyc_dict['time_s'])))

        # check for invalid keys
        assert len(set(new_cyc_dict.keys()) - set(STANDARD_CYCLE_KEYS)) == 0
        return cls(**new_cyc_dict)

    def get_numba_cyc(self):
        """Deprecated."""
        raise NotImplementedError(
            "This method has been deprecated.  Use get_rust_cyc instead.")

    # Properties

    def get_mph(self) -> np.ndarray:
        return self.mps * params.MPH_PER_MPS

    def set_mph(self, new_value):
        self.mps = new_value / params.MPH_PER_MPS

    mph = property(get_mph, set_mph)

    # time step deltas
    @property
    def dt_s(self) -> np.ndarray:
        return np.array(np.diff(self.time_s, prepend=0), dtype=float)

    # distance at each time step
    @property
    def dist_m(self) -> np.ndarray:
        return self.mps * self.dt_s

    @property
    def avg_mps(self):
        return np.append(0.0, 0.5 * (self.mps[1:] + self.mps[:-1]))

    # distance at each time step using average speed of the step
    @property
    def dist_v2_m(self):
        return self.avg_mps * self.dt_s

    @property
    def delta_elev_m(self) -> np.ndarray:
        """
        Cumulative elevation change w.r.t. to initial
        """
        return (self.dist_m * self.grade).cumsum()

    @property
    def len(self) -> int:
        "return cycle length"
        return len(self.time_s)

    def get_cyc_dict(self) -> dict:
        """Returns cycle as dict rather than class instance."""
        keys = STANDARD_CYCLE_KEYS

        cyc = {}
        for key in keys:
            cyc[key] = copy.deepcopy(self.__getattribute__(key))

        return cyc

    def to_rust(self):
        return copy_cycle(self, return_type='rust', deep=False)

    def copy(self):
        """
        Return a copy of this Cycle instance.
        """
        cyc_dict = {
            'time_s': self.time_s.copy(),
            'mps': self.mps.copy(),
            'cycGrade': self.grade.copy(),
            'cycRoadType': self.road_type.copy(),
            'name': self.name,
        }
        return Cycle.from_dict(cyc_dict)

    def average_grade_over_range(self, distance_start_m, delta_distance_m):
        """
        Returns the average grade over the given range of distances
        - distance_start_m: non-negative-number, the distance at start of evaluation area (m)
        - delta_distance_m: non-negative-number, the distance traveled from distance_start_m (m)
        RETURN: number, the average grade (rise over run) over the given distance range
        """
        if ((self.grade == 0.0).all()):
            return 0.0
        distances_m = self.dist_m.cumsum()
        if delta_distance_m <= 0.0:
            return np.interp(distance_start_m, xp=distances_m, fp=self.grade)
        elevations_m = self.delta_elev_m
        e0 = np.interp(distance_start_m, xp=distances_m, fp=elevations_m)
        e1 = np.interp(distance_start_m + delta_distance_m,
                       xp=distances_m, fp=elevations_m)
        return np.tan(np.arcsin((e1 - e0) / delta_distance_m))

    def calc_distance_to_next_stop_from(self, distance_m):
        """
        Calculate the distance to next stop from `distance_m`
        - distance_m: non-negative-number, the current distance from start (m)
        RETURN: -1 or non-negative-integer
        - if there are no more stops ahead, return -1
        - else returns the distance to the next stop from distance_m
        """
        distances_of_stops_m = np.unique(
            self.dist_v2_m.cumsum()[self.mps < 1e-6])
        remaining_stops_m = distances_of_stops_m[distances_of_stops_m > (
            distance_m + 1e-6)]
        if len(remaining_stops_m) > 0:
            return remaining_stops_m[0] - distance_m
        return distances_of_stops_m[-1] - distance_m

    def modify_by_const_jerk_trajectory(self, idx, n, jerk_m__s3, accel0_m__s2):
        """
        Modifies the cycle using the given constant-jerk trajectory parameters
        - idx: non-negative integer, the point in the cycle to initiate
        modification (note: THIS point is modified since trajectory should be calculated from idx-1)
        - jerk_m__s3: number, the "Jerk" associated with the trajectory (m/s3)
        - accel0_m__s2: number, the initial acceleration (m/s2)
        NOTE:
        - modifies cyc in place to hit any critical rendezvous_points by a trajectory adjustment
        - CAUTION: NOT ROBUST AGAINST VARIABLE DURATION TIME-STEPS
        RETURN: Number, final modified speed (m/s)
        """
        num_samples = len(self.time_s)
        v0 = self.mps[idx-1]
        dt = self.dt_s[idx]
        v = v0
        for ni in range(1, int(n)+1):
            idx_to_set = (int(idx) - 1) + ni
            if idx_to_set >= num_samples:
                break
            v = speed_for_constant_jerk(ni, v0, accel0_m__s2, jerk_m__s3, dt)
            self.mps[idx_to_set] = max(v, 0.0)
        return v

    def modify_with_braking_trajectory(self, brake_accel_m__s2, idx):
        """
        Add a braking trajectory that would cover the same distance as the given constant brake deceleration
        - brake_accel_m__s2: negative number, the braking acceleration (m/s2)
        - idx: non-negative integer, the index where to initiate the stop trajectory, start of the step (i in FASTSim)
        RETURN: (non-negative-number, positive-integer)
        - the final speed of the modified trajectory (m/s) 
        - the number of time-steps required to complete the braking maneuver
        NOTE:
        - modifies the cycle in place for the braking trajectory
        """
        assert brake_accel_m__s2 < 0.0
        i = int(idx)
        if i >= len(self.time_s):
            return self.mps[-1], 0
        v0 = self.mps[i-1]
        dt = self.dt_s[i]
        # distance-to-stop (m)
        dts_m = -0.5 * v0 * v0 / brake_accel_m__s2
        # time-to-stop (s)
        tts_s = -v0 / brake_accel_m__s2
        # number of steps to take
        n = int(np.round(tts_s / dt))
        if n < 2:
            # need at least 2 steps
            n = 2
        jerk_m__s3, accel_m__s2 = calc_constant_jerk_trajectory(
            n, 0.0, v0, dts_m, 0.0, dt)
        return self.modify_by_const_jerk_trajectory(i, n, jerk_m__s3, accel_m__s2), n


class LegacyCycle(object):
    """
    Implementation of Cycle with legacy keys.
    """

    def __init__(self, cycle: Cycle):
        """
        Given cycle, returns legacy cycle.
        """
        for key, val in NEW_TO_OLD.items():
            self.__setattr__(val, copy.deepcopy(cycle.__getattribute__(key)))


ref_cyc = Cycle.from_file('udds')


def copy_cycle(cyc: Cycle, return_type: str = None, deep: bool = True) -> Cycle:
    """Returns copy of Cycle.
    Arguments:
    cyc: instantianed Cycle or CycleJit
    return_type: 
        default: infer from type of cyc
        'dict': dict
        'cycle': Cycle 
        'legacy': LegacyCycle
        'rust': RustCycle
    deep: if True, uses deepcopy on everything
    """

    cyc_dict = {}

    for key in inspect_utils.get_attrs(ref_cyc):
        val_to_copy = cyc.__getattribute__(key)
        array_types = [np.ndarray] if not RUST_AVAILABLE else [np.ndarray, fsr.Pyo3ArrayF64]
        if type(val_to_copy) in array_types:
            # has to be float or time_s will get converted to int
            cyc_dict[key] = copy.deepcopy(np.array(
                val_to_copy, dtype=float) if deep else val_to_copy)
        else:
            cyc_dict[key] = copy.deepcopy(val_to_copy) if deep else val_to_copy

    if return_type is None:
        if RUST_AVAILABLE and type(cyc) == fsr.RustCycle:
            return_type = 'rust'
        elif type(cyc) == Cycle:
            return_type = 'cycle'
        elif type(cyc) == LegacyCycle:
            return_type = "legacy"
        else:
            raise NotImplementedError(
                "Only implemented for rust, cycle, or legacy.")

    if return_type == 'dict':
        return cyc_dict
    elif return_type == 'cycle':
        return Cycle.from_dict(cyc_dict)
    elif return_type == 'legacy':
        return LegacyCycle(cyc_dict)
    elif RUST_AVAILABLE and return_type == 'rust':
        return fsr.RustCycle(**cyc_dict)
    else:
        raise ValueError(f"Invalid return_type: '{return_type}'")


def cyc_equal(a: Cycle, b: Cycle) -> bool:
    "Return True if a and b are equal"
    if a is b:
        return True
    a_dict = copy_cycle(a, 'dict')
    b_dict = copy_cycle(b, 'dict')
    return equals(a_dict, b_dict, verbose=False)


def to_microtrips(cycle, stop_speed_m__s=1e-6, keep_name=False):
    """
    Split a cycle into an array of microtrips with one microtrip being a start
    to subsequent stop plus any idle (stopped time).

    Arguments:
    ----------
    cycle: drive cycle converted to dictionary by cycle.get_cyc_dict()
    stop_speed_m__s: speed at which vehicle is considered stopped for trip
        separation
    keep_name: (optional) bool, if True and cycle contains "name", adds
        that name to all microtrips
    """
    microtrips = []
    ts = np.array(cycle['time_s'])
    vs = np.array(cycle['mps'])
    gs = np.array(cycle['grade'])
    rs = np.array(cycle['road_type'])
    mt = make_cycle([], [])
    moving = False
    for idx, (t, v, g, r) in enumerate(zip(ts, vs, gs, rs)):
        if v > stop_speed_m__s and not moving:
            if len(mt['time_s']) > 1:
                temp = make_cycle(
                    [mt['time_s'][-1]],
                    [mt['mps'][-1]],
                    [mt['grade'][-1]],
                    [mt['road_type'][-1]])
                mt['time_s'] = mt['time_s'] - mt['time_s'][0]
                for k in mt:
                    mt[k] = np.array(mt[k])
                microtrips.append(mt)
                mt = temp
        mt['time_s'] = np.append(mt['time_s'], [t])
        mt['mps'] = np.append(mt['mps'], [v])
        mt['grade'] = np.append(mt['grade'], [g])
        mt['road_type'] = np.append(mt['road_type'], [r])
        moving = v > stop_speed_m__s
    if len(mt['time_s']) > 0:
        mt['time_s'] = mt['time_s'] - mt['time_s'][0]
        microtrips.append(mt)
    if keep_name and "name" in cycle:
        for m in microtrips:
            m["name"] = cycle["name"]
    return microtrips


def make_cycle(ts, vs, gs=None, rs=None) -> dict:
    """
    (Array Num) (Array Num) (Array Num)? -> Dict
    Create a cycle from times, speeds, and grades. If grades is not
    specified, it is set to zero.
    Arguments:
    ----------
    ts: array of times [s]
    vs: array of vehicle speeds [mps]
    gs: array of grades
    rs: array of road types (charging or not)
    """
    assert len(ts) == len(vs)
    if gs is None:
        gs = np.zeros(len(ts))
    else:
        assert len(ts) == len(gs)
    if rs is None:
        rs = np.zeros(len(ts))
    else:
        assert len(ts) == len(rs)
    return {'time_s': np.array(ts),
            'mps': np.array(vs),
            'grade': np.array(gs),
            'road_type': np.array(rs)}


def equals(c1, c2, verbose=True) -> bool:
    """
    Dict Dict -> Bool
    Returns true if the two cycles are equal, false otherwise
    Arguments:
    ----------
    c1: cycle as dictionary from get_cyc_dict()
    c2: cycle as dictionary from get_cyc_dict()
    verbose: Bool, optional (default: True), if True, prints why not equal
    """
    if c1.keys() != c2.keys():
        if verbose:
            c2missing = set(c1.keys()) - set(c2.keys())
            c1missing = set(c2.keys()) - set(c1.keys())
            if len(c1missing) > 0:
                print('c2 keys not contained in c1: {}'.format(c1missing))
            if len(c2missing) > 0:
                print('c1 keys not contained in c2: {}'.format(c2missing))
        return False
    for k in c1.keys():
        if len(c1[k]) != len(c2[k]):
            if verbose:
                print(k + ' has a length discrepancy.')
            return False
        if np.any(np.array(c1[k]) != np.array(c2[k])):
            if verbose:
                print(k + ' has a value discrepancy.')
            return False
    return True


def concat(cycles, name=None):
    """
    Concatenates cycles together one after another into a single dictionary
    (Array Dict) String -> Dict
    Arguments:
    ----------
    cycles: (Array Dict)
    name: (optional) string or None, if a string, adds the "name" key to the output
    """
    final_cycle = {'time_s': np.array([]),
                   'mps': np.array([]),
                   'grade': np.array([]),
                   'road_type': np.array([])}
    keys = [k for k in final_cycle.keys()]
    first = True
    for cycle in cycles:
        if first:
            for k in keys:
                final_cycle[k] = np.array(cycle[k])
            first = False
        # if len(final_cycle['time_s']) == 0: # ???
        #     t0 = 0.0
        else:
            for k in keys:
                if k == 'time_s':
                    t0 = final_cycle[k][-1]
                    final_cycle[k] = np.concatenate([
                        final_cycle[k], np.array(cycle[k][1:]) + t0])
                else:
                    final_cycle[k] = np.concatenate([
                        final_cycle[k], np.array(cycle[k][1:])])
    if name is not None:
        final_cycle["name"] = name
    return final_cycle


def resample(cycle: Cycle, new_dt=None, start_time=None, end_time=None,
             hold_keys=None, rate_keys=None):
    """
    Cycle new_dt=?Real start_time=?Real end_time=?Real -> Cycle
    Resample a cycle with a new delta time from start time to end time.

    - cycle: Dict with keys
        'time_s': numpy.array Real giving the elapsed time
    - new_dt: Real, optional
        the new delta time of the sampling. Defaults to the
        difference between the first two times of the cycle passed in
    - start_time: Real, optional
        the start time of the sample. Defaults to 0.0 seconds
    - end_time: Real, optional
        the end time of the cycle. Defaults to the last time of the passed in
        cycle.
    - hold_keys: None or (Set String), if specified, yields values that
                 should be interpolated step-wise, holding their value until
                 an explicit change (i.e., NOT interpolated)
    - rate_keys: None or (Set String), if specified, yields values that maintain
                 the interpolated value of the given rate. So, for example,
                 if a speed, will set the speed such that the distance traveled
                 is consistent. Note: using rate keys for mps may result in
                 non-zero starting and ending speeds
    Resamples all non-time metrics by the new sample time.
    """
    if new_dt is None:
        new_dt = cycle['time_s'][1] - cycle['time_s'][0]
    if start_time is None:
        start_time = 0.0
    if end_time is None:
        end_time = cycle['time_s'][-1]
    new_cycle = {}
    eps = new_dt / 10.0
    new_cycle['time_s'] = np.arange(start_time, end_time + eps, step=new_dt)
    for k in cycle:
        if k == 'time_s':
            continue
        elif hold_keys is not None and k in hold_keys:
            f = interp1d(cycle['time_s'], cycle[k], 0)
            new_cycle[k] = f(new_cycle['time_s'])
            continue
        elif rate_keys is not None and k in rate_keys:
            rate_var = np.array(cycle[k])
            # below is same as [(rate_var[i+1] + rate_var[i])/2.0 for i in range(len(rate_var) - 1)]
            avg_rate_var = (rate_var[1:] + rate_var[:-1]) / 2.0
            dts_orig = np.diff(cycle['time_s'])
            step_averages = np.diff(
                np.interp(
                    x=new_cycle['time_s'],
                    xp=cycle['time_s'],
                    fp=np.insert((avg_rate_var * dts_orig).cumsum(), 0, 0.0)
                ),
            ) / new_dt
            step_averages = np.append(step_averages[0], step_averages)
            step_averages = np.append(step_averages, step_averages[-1])
            midstep_times = np.concatenate(
                (
                    [0.0],
                    np.arange(
                        start_time + (0.5 * new_dt), end_time - (0.5 * new_dt) + eps, step=new_dt),
                    [end_time],
                ))
            new_cycle[k] = np.interp(
                x=new_cycle['time_s'],
                xp=midstep_times,
                fp=step_averages
            )
            continue
        try:
            new_cycle[k] = np.interp(
                new_cycle['time_s'], cycle['time_s'], cycle[k])
        except:
            # if the value can't be interpolated, it must not be a numerical
            # array. Just add it back in as is.
            new_cycle[k] = deepcopy(cycle[k])
    return new_cycle


def clip_by_times(cycle, t_end, t_start=0):
    """
    Cycle Number Number -> Cycle
    INPUT:
    - cycle: Dict, a legitimate driving cycle
    - t_start: Number, time to start
    - t_end: Number, time to end
    RETURNS: Dict, the cycle with fields snipped
        to times >= t_start and <= t_end
    Clip the cycle to the given times and return
    """
    idx = np.logical_and(np.array(cycle['time_s']) >= t_start,
                         np.array(cycle['time_s']) <= t_end)
    new_cycle = {}
    for k in cycle:
        try:
            new_cycle[k] = np.array(cycle[k])[idx]
        except:
            new_cycle[k] = cycle[k]

    # reset time to start at zero
    new_cycle['time_s'] -= new_cycle['time_s'][0]
    return new_cycle


def accelerations(cycle):
    """
    Cycle -> Real
    Return the acceleration of the given cycle
    INPUTS:
    - cycle: Dict, a legitimate driving cycle
    OUTPUTS: Real, the maximum acceleration
    """
    accels = (np.diff(np.array(cycle['mps']))
              / np.diff(np.array(cycle['time_s'])))
    return accels


def peak_acceleration(cycle):
    """
    Cycle -> Real
    Return the maximum acceleration of the given cycle
    INPUTS:
    - cycle: Dict, a legitimate driving cycle
    OUTPUTS: Real, the maximum acceleration
    """
    return np.max(accelerations(cycle))


def peak_deceleration(cycle):
    """
    Cycle -> Real
    Return the minimum acceleration (maximum deceleration) of the given cycle
    INPUTS:
    - cycle: Dict, a legitimate driving cycle
    OUTPUTS: Real, the maximum acceleration
    """
    return np.min(accelerations(cycle))


def calc_constant_jerk_trajectory(n: int, D0: float, v0: float, Dr: float, vr: float, dt: float) -> tuple:
    """
    Num Num Num Num Num Int -> (Dict 'jerk_m__s3' Num 'accel_m__s2' Num)
    INPUTS:
    - n: Int, number of time-steps away from rendezvous
    - D0: Num, distance of simulated vehicle (m/s)
    - v0: Num, speed of simulated vehicle (m/s)
    - Dr: Num, distance of rendezvous point (m)
    - vr: Num, speed of rendezvous point (m/s)
    - dt: Num, step duration (s)
    RETURNS: (Tuple 'jerk_m__s3': Num, 'accel_m__s2': Num)
    Returns the constant jerk and acceleration for initial time step.
    """
    assert n > 1, f"n = {n}"
    assert Dr > D0, f"Dr = {Dr}; D0 = {D0}"
    dDr = Dr - D0
    dvr = vr - v0
    k = (dvr - (2.0 * dDr / (n * dt)) + 2.0 * v0) / (
        0.5 * n * (n - 1) * dt
        - (1.0 / 3) * (n - 1) * (n - 2) * dt
        - 0.5 * (n - 1) * dt * dt
    )
    a0 = (
        (dDr / dt)
        - n * v0
        - ((1.0 / 6) * n * (n - 1) * (n - 2) *
           dt + 0.25 * n * (n - 1) * dt * dt) * k
    ) / (0.5 * n * n * dt)
    return (k, a0)


def accel_for_constant_jerk(n, a0, k, dt):
    """
    Calculate the acceleration n timesteps away
    INPUTS:
    - n: Int, number of times steps away to calculate
    - a0: Num, initial acceleration (m/s2)
    - k: Num, constant jerk (m/s3)
    - dt: Num, time-step duration in seconds
    NOTE:
    - this is the constant acceleration over the time-step from sample n to sample n+1
    RETURN: Num, the acceleration n timesteps away (m/s2)
    """
    return a0 + (n * k * dt)


def speed_for_constant_jerk(n, v0, a0, k, dt):
    """
    Int Num Num Num Num -> Num
    Calculate speed (m/s) n timesteps away
    INPUTS:
    - n: Int, numer of timesteps away to calculate
    - v0: Num, initial speed (m/s)
    - a0: Num, initial acceleration (m/s2)
    - k: Num, constant jerk
    - dt: Num, duration of a timestep (s)
    NOTE:
    - this is the speed at sample n
    - if n == 0, speed is v0
    - if n == 1, speed is v0 + a0*dt, etc.
    RETURN: Num, the speed n timesteps away (m/s)
    """
    return v0 + (n * a0 * dt) + (0.5 * n * (n - 1) * k * dt)


def dist_for_constant_jerk(n, d0, v0, a0, k, dt):
    """
    Calculate distance (m) after n timesteps
    INPUTS:
    - n: Int, numer of timesteps away to calculate
    - d0: Num, initial distance (m)
    - v0: Num, initial speed (m/s)
    - a0: Num, initial acceleration (m/s2)
    - k: Num, constant jerk
    - dt: Num, duration of a timestep (s)
    NOTE:
    - this is the distance traveled from start (i.e., n=0) measured at sample point n
    RETURN: Num, the distance at n timesteps away (m)
    """
    term1 = dt * (
        (n * v0)
        + (0.5 * n * (n - 1) * a0 * dt)
        + ((1.0 / 6.0) * k * dt * (n - 2) * (n - 1) * n)
    )
    term2 = 0.5 * dt * dt * ((n * a0) + (0.5 * n * (n - 1) * k * dt))
    return d0 + term1 + term2


def detect_passing(cyc: Cycle, cyc0: Cycle, i: int, tol: float=0.1, vmin: Optional[float] = None, dtmax:Optional[float] = None) -> dict:
    """
    Reports back information of the first point where cyc passes cyc0, starting at step i until the next stop of cyc.
    - cyc: fastsim.Cycle, the proposed cycle of the vehicle under simulation
    - cyc0: fastsim.Cycle, the reference/lead vehicle/shadow cycle to compare with
    - i: int, the time-step index to consider
    - tol: float, the distance tolerance away from lead vehicle to be seen as "deviated" from the reference/shadow trace (m)
    - vmin: float, the minimum speed to consider for re-rendezvous (m/s)
    - dtmax: float, how far ahead in time to look (s)
    RETURNS: dict with keys:
        "has_collision": bool, True if cyc passes cyc0
        "idx": int, the index where cyc passes cyc0
        "num_steps": int, the number of time-steps until idx from i
        "distance_m": float, the distance (m) traveled of cyc0 when cyc passes
        "speed_m_per_s": float, the speed (m/s) of cyc0 when cyc passes
    """
    v0 = cyc.mps[i-1]
    d0 = cyc.dist_v2_m[:i].sum()
    v0_lv = cyc0.mps[i-1]
    d0_lv = cyc0.dist_v2_m[:i].sum()
    d = d0
    d_lv = d0_lv
    dt_total = 0.0
    rendezvous_idx = None
    rendezvous_num_steps = 0
    rendezvous_distance_m = 0
    rendezvous_speed_m_per_s = 0
    for di in range(len(cyc.mps) - i):
        idx = i + di
        v = cyc.mps[idx]
        if v == 0.0:
            break
        v_lv = cyc0.mps[idx]
        vavg = (v + v0) * 0.5
        vavg_lv = (v_lv + v0_lv) * 0.5
        dd = vavg * cyc.dt_s[idx]
        dd_lv = vavg_lv * cyc0.dt_s[idx]
        dt_total += cyc0.dt_s[idx]
        d += dd
        d_lv += dd_lv
        dtlv = d_lv - d
        v0 = v
        v0_lv = v_lv
        if di > 0 and (vmin is None or v_lv >= vmin) and dtlv < -tol:
            rendezvous_idx = idx
            rendezvous_num_steps = di + 1
            rendezvous_distance_m = d_lv
            rendezvous_speed_m_per_s = v_lv
            break
        if dtmax is not None and dt_total > dtmax:
            break
    return {
        "has_collision": rendezvous_idx is not None and rendezvous_distance_m > d0,
        "idx": rendezvous_idx,
        "num_steps": rendezvous_num_steps,
        "start_distance_m": d0,
        "distance_m": rendezvous_distance_m,
        "start_speed_m_per_s": cyc.mps[i-1],
        "time_step_duration_s": cyc.dt_s[i],
        "speed_m_per_s": rendezvous_speed_m_per_s,
    }