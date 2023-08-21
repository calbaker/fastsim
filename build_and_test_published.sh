# build and test with published version of `fastsim-proc-macros`
(cd rust/fastsim-core/ && cargo test) && \
(cd rust/fastsim-cli/ && cargo test) && \
pip install -qe ".[dev]" && \
pytest -v python/fastsim/tests/