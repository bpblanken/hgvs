"""Defines the abstract data provider interface"""

import abc
import dataclasses
import functools
import logging
import os
import typing
import warnings
from typing import ClassVar

import hgvs

from ..decorators.lru_cache import LEARN, RUN, VERIFY, lru_cache
from ..exceptions import HGVSDataNotAvailableError
from ..utils.PersistentDict import PersistentDict
from .models import TranscriptInfo

_logger = logging.getLogger(__name__)


def _cached_method_name(key):
    """Return the provider method name recorded in an lru_cache key, or None."""
    try:
        parts = list(key)
        return parts[parts.index("__func__") + 1]
    except (TypeError, ValueError, IndexError):
        return None


def _upgrade_legacy_cache(cache):
    """Convert cached dict results, recorded before data models existed, to models.

    Runs once when a persistent cache is opened, so caches written by earlier
    hgvs releases keep working in every cache mode.  In learn mode the
    converted entries are written back to the file on the next sync.
    """
    for key, value in list(cache.items()):
        model = _RESULT_MODELS.get(_cached_method_name(key))
        if model is None or isinstance(value, model):
            continue
        if not isinstance(value, dict):
            msg = f"cached result for {key!r} is {type(value).__name__}; expected {model.__name__}"
            raise TypeError(msg)
        try:
            cache[key] = model(**value)
        except TypeError as e:
            msg = f"cannot convert cached result for {key!r} to {model.__name__}: {e}"
            raise TypeError(msg) from e


def _coerce_to_model(method, model):
    """Wrap a provider method so that a legacy dict result becomes ``model``.

    Before data models existed, providers returned plain dicts keyed by
    the UTA column names, and some returned ``None`` for "not found".
    Providers that still do so keep working for now: dicts are converted
    with a :class:`DeprecationWarning`, and ``None`` raises
    :class:`HGVSDataNotAvailableError` as the mappers used to.  This
    adapter will be removed in a future release.
    """

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        result = method(*args, **kwargs)
        if isinstance(result, model):
            return result
        provider = type(method.__self__).__name__
        if result is None:
            msg = f"{provider}.{method.__name__}{args} returned no data"
            raise HGVSDataNotAvailableError(msg)
        if isinstance(result, dict):
            warnings.warn(
                f"{provider}.{method.__name__} returned a dict; data providers must return"
                f" hgvs.dataproviders.models.{model.__name__}. Dict results are deprecated"
                " and support will be removed in a future release.",
                DeprecationWarning,
                # attribute the warning to the caller of hdp.<method>(), skipping
                # this wrapper and the lru_cache wrapper around it
                stacklevel=3,
            )
            return model(**result)
        msg = (
            f"{provider}.{method.__name__} returned {type(result).__name__};"
            f" expected hgvs.dataproviders.models.{model.__name__}"
        )
        raise TypeError(msg)

    return wrapper


class Interface(metaclass=abc.ABCMeta):
    """Variant mapping and validation requires access to external data,
    specifically exon structures, transcript alignments, and protein
    accessions.  In order to isolate the hgvs package from the myriad
    choices and tradeoffs, these data are provided through an
    implementation of the (abstract) HGVS Data Provider Interface.

    As of June 2014, the only available data provider implementation
    uses the `Universal Transcript Archive (UTA)`_, a sister project
    that provides access to transcripts and genome-transcript
    alignments.  `Invitae`_ provides a public UTA database instance
    that is used by default; see the `UTA`_ page for instructions on
    installing your own PostgreSQL or SQLite version.  In the future,
    other implementations may be availablefor other data sources.

    Pure virtural class for the HGVS Data Provider Interface.  Every
    data provider implementation should be a subclass (possibly
    indirect) of this class.

    .. _`Universal Transcript Archive (UTA)`: https://github.com/biocommons/uta
    .. _Invitae: http://invitae.com/
    """

    def interface_version(self):
        return 5

    def __init__(self, mode=None, cache=None):
        """
        :param mode: cache mode (None[default lru cache], 'learn', 'run', 'verify')
        :type mode: str
        :param cache: local cache file name
        :type cache: str
        """
        self.mode = None
        if mode == "learn":
            self.mode = LEARN
        elif mode == "run":
            self.mode = RUN
        elif mode == "verify":
            self.mode = VERIFY

        self.cache = None
        if self.mode is not None:
            if self.mode == LEARN:
                self.cache = PersistentDict(cache, flag="c")
            else:
                self.cache = PersistentDict(cache, flag="r")
            _upgrade_legacy_cache(self.cache)

        maxsize = hgvs.global_config.lru_cache.maxsize
        if "PYTEST_CURRENT_TEST" in os.environ:
            maxsize = None
            _logger.info("%s: Using unlimited cache size", __file__)

        self.data_version = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.data_version
        )
        self.schema_version = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.schema_version
        )
        self.get_acs_for_protein_seq = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_acs_for_protein_seq
        )
        self.get_gene_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_gene_info
        )
        self.get_pro_ac_for_tx_ac = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_pro_ac_for_tx_ac
        )
        self.get_seq = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(self.get_seq)
        self.get_similar_transcripts = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_similar_transcripts
        )
        self.get_tx_exons = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_exons
        )
        self.get_tx_for_gene = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_for_gene
        )
        self.get_tx_for_region = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_for_region
        )
        self.get_tx_identity_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_identity_info
        )
        self.get_tx_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            _coerce_to_model(self.get_tx_info, _RESULT_MODELS["get_tx_info"])
        )
        self.get_tx_mapping_options = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_mapping_options
        )

        def _split_version_string(v):
            versions = list(map(int, v.split(".")))
            if len(versions) < 2:
                versions += [0]
            assert len(versions) == 2
            return versions

        assert self.required_version is not None

        rv = _split_version_string(self.required_version)
        av = _split_version_string(self.schema_version())

        if av[0] == rv[0] and av[1] >= rv[1]:
            return

        msg = f"Incompatible versions: {type(self).__name__} requires schema version {self.required_version}, but {self.url} provides version {self.schema_version()}"
        raise RuntimeError(msg)

    # required_version: what version of the remote schema is required
    # by the subclass? This value is compared to the result of
    # schema_version, which must be implemented by the class.
    required_version: ClassVar[str | None] = None

    @abc.abstractmethod
    def data_version(self):
        pass

    @abc.abstractmethod
    def schema_version(self):
        pass

    @abc.abstractmethod
    def get_acs_for_protein_seq(self, seq):
        pass

    @abc.abstractmethod
    def get_assembly_map(self, assembly_name):
        pass

    @abc.abstractmethod
    def get_gene_info(self, gene):
        pass

    @abc.abstractmethod
    def get_pro_ac_for_tx_ac(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_seq(self, ac, start_i=None, end_i=None):
        pass

    @abc.abstractmethod
    def get_similar_transcripts(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        pass

    @abc.abstractmethod
    def get_tx_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        pass

    @abc.abstractmethod
    def get_tx_identity_info(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_tx_info(self, tx_ac: str, alt_ac: str, alt_aln_method: str) -> TranscriptInfo:
        pass

    @abc.abstractmethod
    def get_tx_mapping_options(self, tx_ac):
        pass


def _result_models(cls):
    """Map each abstract method whose return annotation is a data model to that model."""
    models = {}
    for name in cls.__abstractmethods__:
        rtype = typing.get_type_hints(getattr(cls, name)).get("return")
        if isinstance(rtype, type) and dataclasses.is_dataclass(rtype):
            models[name] = rtype
    return models


# Provider methods whose results are data models, read from the annotations above.
_RESULT_MODELS = _result_models(Interface)


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
