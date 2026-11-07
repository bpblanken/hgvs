"""Tests for the abstract data provider interface and the version 4 bridge"""

import pytest

from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.models import TranscriptExon, TranscriptIdentityInfo, TranscriptInfo
from hgvs.exceptions import HGVSDataNotAvailableError

TX_INFO_DICT = {
    "hgnc": "ATM",
    "cds_start_i": 385,
    "cds_end_i": 9556,
    "tx_ac": "NM_000051.3",
    "alt_ac": "AC_000143.1",
    "alt_aln_method": "splign",
}
TX_INFO = TranscriptInfo(**TX_INFO_DICT)
# a version 4 exon row carries UTA-internal columns the model does not declare
EXON_DICT = {
    "tx_ac": "NM_000051.3",
    "alt_ac": "AC_000143.1",
    "alt_aln_method": "splign",
    "alt_strand": 1,
    "ord": 0,
    "tx_start_i": 0,
    "tx_end_i": 100,
    "alt_start_i": 1000,
    "alt_end_i": 1100,
    "cigar": "100=",
    "hgnc": "ATM",
    "tx_exon_id": 7,
    "exon_aln_id": 8,
    "tx_aseq": None,
}
TX_IDENTITY_DICT = {
    "tx_ac": "NM_000051.3",
    "alt_ac": "NM_000051.3",
    "alt_aln_method": "transcript",
    "cds_start_i": 385,
    "cds_end_i": 9556,
    "lengths": [100, 200],
    "hgnc": "ATM",
}
TX_INFO_ARGS = ("NM_000051.3", "AC_000143.1", "splign")


class V4Provider(Interface):
    """A provider written against interface version 4: get_tx_info returns a dict."""

    required_version = "1.1"
    url = "fake://v4"

    def __init__(self, tx_info_result):
        self._tx_info_result = tx_info_result
        self.calls = []
        super().__init__()

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self.calls.append(("get_tx_info", tx_ac, alt_ac, alt_aln_method))
        return self._tx_info_result

    def data_version(self):
        return "fake"

    def schema_version(self):
        return "1.1"

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        self.calls.append(("get_tx_exons", tx_ac, alt_ac, alt_aln_method))
        return [dict(EXON_DICT, tx_ac=tx_ac, alt_ac=alt_ac, alt_aln_method=alt_aln_method)]

    def get_tx_identity_info(self, tx_ac):
        return {**TX_IDENTITY_DICT, "tx_ac": tx_ac, "alt_ac": tx_ac}

    def get_tx_mapping_options(self, tx_ac):
        return [{"tx_ac": tx_ac, "alt_ac": "AC_000143.1", "alt_aln_method": "splign"}]

    # the remaining abstract methods are not exercised here
    def get_acs_for_protein_seq(self, seq):
        raise NotImplementedError

    def get_assembly_map(self, assembly_name):
        raise NotImplementedError

    def get_gene_info(self, gene):
        raise NotImplementedError

    def get_pro_ac_for_tx_ac(self, tx_ac):
        raise NotImplementedError

    def get_seq(self, ac, start_i=None, end_i=None):
        raise NotImplementedError

    def get_similar_transcripts(self, tx_ac):
        raise NotImplementedError

    def get_tx_for_gene(self, gene):
        raise NotImplementedError

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError


def test_bridge_is_a_version_5_interface():
    hdp = Interface.from_v4(V4Provider(dict(TX_INFO_DICT)))
    assert isinstance(hdp, Interface)
    assert hdp.interface_version() == 5


def test_bridge_converts_dict_tx_info_to_model():
    provider = V4Provider(dict(TX_INFO_DICT))
    hdp = Interface.from_v4(provider)
    assert hdp.get_tx_info(*TX_INFO_ARGS) == TX_INFO
    assert provider.calls == [("get_tx_info", *TX_INFO_ARGS)]


def test_bridge_turns_none_tx_info_into_data_not_available():
    hdp = Interface.from_v4(V4Provider(None))
    with pytest.raises(
        HGVSDataNotAvailableError, match=r"V4Provider\.get_tx_info.*returned no data"
    ):
        hdp.get_tx_info(*TX_INFO_ARGS)


def test_bridge_reads_only_declared_keys_and_reports_missing_ones():
    hdp = Interface.from_v4(V4Provider({**TX_INFO_DICT, "tx_exon_set_id": 1}))
    assert hdp.get_tx_info(*TX_INFO_ARGS) == TX_INFO
    hdp = Interface.from_v4(V4Provider({k: v for k, v in TX_INFO_DICT.items() if k != "hgnc"}))
    with pytest.raises(KeyError, match="hgnc"):
        hdp.get_tx_info(*TX_INFO_ARGS)


def test_bridge_delegates_unchanged_methods_as_is():
    hdp = Interface.from_v4(V4Provider(dict(TX_INFO_DICT)))
    assert hdp.schema_version() == "1.1"
    assert hdp.data_version() == "fake"


def test_bridge_converts_list_results_and_ignores_undeclared_keys():
    hdp = Interface.from_v4(V4Provider(dict(TX_INFO_DICT)))
    exons = hdp.get_tx_exons(*TX_INFO_ARGS)
    assert exons == [
        TranscriptExon(
            tx_ac="NM_000051.3",
            alt_ac="AC_000143.1",
            alt_aln_method="splign",
            alt_strand=1,
            ord=0,
            tx_start_i=0,
            tx_end_i=100,
            alt_start_i=1000,
            alt_end_i=1100,
            cigar="100=",
            hgnc="ATM",
        )
    ]
    assert [o.alt_ac for o in hdp.get_tx_mapping_options("NM_000051.3")] == ["AC_000143.1"]


def test_bridge_converts_identity_info_with_tuple_lengths():
    hdp = Interface.from_v4(V4Provider(dict(TX_INFO_DICT)))
    info = hdp.get_tx_identity_info("NM_000051.3")
    assert isinstance(info, TranscriptIdentityInfo)
    assert info.lengths == (100, 200)


def test_bridge_exposes_provider_attributes():
    provider = V4Provider(dict(TX_INFO_DICT))
    hdp = Interface.from_v4(provider)
    assert hdp.url == "fake://v4"
    assert hdp.required_version == "1.1"
    with pytest.raises(AttributeError):
        _ = hdp.no_such_attribute
