from rep2struct.schema import FoldJob
from rep2struct import grouping


def _job(cid, mhc_class=1, has_tcr=True, species="human", output_needed="structure"):
    return FoldJob(clonotype_id=cid, construct_fasta=">A\nAAAA",
                   mhc_class=mhc_class, has_tcr=has_tcr, species=species,
                   output_needed=output_needed)


def test_group_key_is_stable_and_descriptive():
    assert grouping.group_key(_job("x")) == "c1_tcr_human_structure"
    assert grouping.group_key(_job("y", mhc_class=2, has_tcr=False,
                                    species="mouse", output_needed="binding_score")) \
        == "c2_notcr_mouse_binding_score"


def test_partition_splits_and_stamps_group_id():
    jobs = [_job("a"), _job("b"), _job("c", mhc_class=2)]
    groups = grouping.partition(jobs)
    assert set(groups) == {"c1_tcr_human_structure", "c2_tcr_human_structure"}
    assert len(groups["c1_tcr_human_structure"]) == 2
    assert all(j.group_id == "c1_tcr_human_structure" for j in groups["c1_tcr_human_structure"])
    assert groups["c2_tcr_human_structure"][0].group_id == "c2_tcr_human_structure"
