"""Tests for the abstract data provider interface"""

import pickle

import pytest

from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.models import TranscriptInfo
from hgvs.decorators.lru_cache import _make_key
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


class FakeProvider(Interface):
    """Minimal provider whose get_tx_info returns whatever it is told to."""

    required_version = "1.1"

    def __init__(self, tx_info_result, mode=None, cache=None):
        self._tx_info_result = tx_info_result
        self.calls = 0
        super().__init__(mode=mode, cache=cache)

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):  # noqa: ARG002
        self.calls += 1
        return self._tx_info_result

    def data_version(self):
        return "fake"

    def schema_version(self):
        return "1.1"

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

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        raise NotImplementedError

    def get_tx_for_gene(self, gene):
        raise NotImplementedError

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        raise NotImplementedError

    def get_tx_identity_info(self, tx_ac):
        raise NotImplementedError

    def get_tx_mapping_options(self, tx_ac):
        raise NotImplementedError


def test_model_result_passes_through_without_warning(recwarn):
    hdp = FakeProvider(TX_INFO)
    assert hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign") is TX_INFO
    assert not [w for w in recwarn if issubclass(w.category, DeprecationWarning)]


def test_legacy_dict_result_is_coerced_with_warning():
    hdp = FakeProvider(dict(TX_INFO_DICT))
    with pytest.warns(DeprecationWarning, match="FakeProvider.get_tx_info returned a dict"):
        result = hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")
    assert result == TX_INFO


def test_legacy_none_result_raises_data_not_available():
    hdp = FakeProvider(None)
    with pytest.raises(HGVSDataNotAvailableError, match="returned no data"):
        hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")


def test_unexpected_result_type_raises_type_error():
    hdp = FakeProvider(("NM_000051.3", 385, 9556))
    with pytest.raises(TypeError, match=r"expected hgvs\.dataproviders\.models\.TranscriptInfo"):
        hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")


def test_coerced_result_is_cached_as_model():
    hdp = FakeProvider(dict(TX_INFO_DICT))
    with pytest.warns(DeprecationWarning, match="returned a dict"):
        first = hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")
    second = hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")
    assert second is first
    assert hdp.calls == 1
    assert hdp.get_tx_info.cache_info().hits == 1


def _write_legacy_cache(path, tx_info_value):
    """Write a persistent cache as a pre-2.0 hgvs would have: dict results."""
    entries = {
        _make_key("schema_version", (), {}, False, ()): "1.1",
        _make_key(
            "get_tx_info", ("NM_000051.3", "AC_000143.1", "splign"), {}, False, ()
        ): tx_info_value,
    }
    with path.open("wb") as f:
        pickle.dump(entries, f, 4)


def test_legacy_cache_dict_is_converted_as_it_leaves_the_cache(tmp_path, recwarn):
    cache = tmp_path / "cache.hdp"
    _write_legacy_cache(cache, dict(TX_INFO_DICT))
    hdp = FakeProvider(None, mode="run", cache=str(cache))
    result = hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")
    assert result == TX_INFO
    assert hdp.calls == 0, "run mode must be served from the cache"
    assert not [w for w in recwarn if issubclass(w.category, DeprecationWarning)]


def test_legacy_cache_dict_verifies_against_a_model_result(tmp_path):
    cache = tmp_path / "cache.hdp"
    _write_legacy_cache(cache, dict(TX_INFO_DICT))
    hdp = FakeProvider(TX_INFO, mode="verify", cache=str(cache))
    assert hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign") == TX_INFO
    assert hdp.calls == 1, "verify mode must call the provider"


def test_legacy_cache_dict_with_unknown_keys_fails_on_read(tmp_path):
    cache = tmp_path / "cache.hdp"
    _write_legacy_cache(cache, {**TX_INFO_DICT, "bogus": 1})
    hdp = FakeProvider(None, mode="run", cache=str(cache))
    with pytest.raises(TypeError, match="bogus"):
        hdp.get_tx_info("NM_000051.3", "AC_000143.1", "splign")
