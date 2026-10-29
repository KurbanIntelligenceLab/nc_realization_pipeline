"""Fig. 1 draws an audited record as an open circle, and draws none without an audit file."""
import numpy as np

from ncr import figures as FIG


def _res():
    # one release, two later records 1.211 and 2.2123 yr after 2014.2849315068493 (offsets rounded to 1e-4)
    return {"fig1": {"ledger": {"R (2014)": {"date": 2014.2849315068493, "n": 10, "known_years": [1990.5, 2001.5],
                                            "later_offsets": [1.211, 2.2123], "t_max": 10, "age_median": 18}}}}


def _open_points(fig):
    n = 0
    for ax in fig.axes:
        for coll in ax.collections:
            fc = coll.get_facecolors()
            if len(fc) and np.allclose(fc[0][:3], 1.0):
                n += len(coll.get_offsets())
    return n


def test_audited_record_is_an_open_circle(tmp_path):
    audit = [{"release": "R (2014)", "record_date": "2016-07-01", "in_cohort_10yr": "yes",
              "finding": "earlier experimental report"}]
    f = FIG.fig1(_res(), str(tmp_path / "f.pdf"), audit)
    assert _open_points(f) == 1
    assert f._ledger_counts["R (2014)"] == (2, 2)          # the audited record is still counted


def test_no_audit_means_no_open_circles(tmp_path):
    f = FIG.fig1(_res(), str(tmp_path / "f.pdf"))
    assert _open_points(f) == 0
