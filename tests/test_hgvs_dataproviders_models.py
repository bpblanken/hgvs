"""Tests for the data provider models"""

import pickle
from dataclasses import asdict

import pytest

from hgvs.dataproviders.models import TranscriptIdentityInfo, TranscriptInfo, from_mapping


@pytest.fixture
def tx_info():
    return TranscriptInfo(
        tx_ac="NM_000051.3",
        alt_ac="AC_000143.1",
        alt_aln_method="splign",
        cds_start_i=385,
        cds_end_i=9556,
        hgnc="ATM",
    )


def test_transcript_info_is_frozen(tx_info):
    with pytest.raises(AttributeError):
        tx_info.cds_start_i = 0


def test_transcript_info_round_trips_through_pickle(tx_info):
    assert pickle.loads(pickle.dumps(tx_info)) == tx_info  # noqa: S301 - trusted, built above


def test_transcript_info_has_no_dict_access(tx_info):
    with pytest.raises(TypeError):
        tx_info["cds_start_i"]


def test_from_mapping_ignores_undeclared_keys(tx_info):
    row = {**asdict(tx_info), "tx_exon_set_id": 42}
    assert from_mapping(TranscriptInfo, row) == tx_info


def test_from_mapping_reports_missing_field(tx_info):
    row = {k: v for k, v in asdict(tx_info).items() if k != "hgnc"}
    with pytest.raises(KeyError, match="hgnc"):
        from_mapping(TranscriptInfo, row)


def test_identity_info_lengths_are_a_tuple():
    info = TranscriptIdentityInfo("NM_1.1", "NM_1.1", "transcript", 0, 9, [3, 3, 3], None)
    assert info.lengths == (3, 3, 3)
    assert sum(info.lengths) == 9
