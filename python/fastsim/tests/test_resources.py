"""Test getting resource lists via Rust API"""
import fastsim as fsim


def test_list_resources_for_cycle():
    """
    Assert list_resources works for Cycle
    """
    cyc = fsim.Cycle.from_resource("udds.csv")
    resources = cyc.list_resources()
    assert len(resources) > 0
