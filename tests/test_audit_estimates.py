"""The audit's composition rule matches the pipeline's tolerance and element-set requirement."""
from scripts_path import *  # noqa: F401,F403  (adds scripts/ to the path)
import audit_estimates as A


def test_within_rule():
    assert A.within_rule("CdAg3", "Ag0.7424Cd0.2576")      # distance 0.015
    assert A.within_rule("Zn3Cu", "Cu0.255Zn0.745")        # distance 0.010
    assert not A.within_rule("ZnCuAu2", "Cu0.30Au0.50Zn0.20")  # distance 0.10
    assert not A.within_rule("V12P7", "V5P3C")             # different element set
    assert not A.within_rule("YbSm3S4", "")                # no supported composition
