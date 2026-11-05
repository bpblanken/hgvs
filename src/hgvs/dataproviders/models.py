"""Data models returned by HGVS data providers.

These models are the contract between the hgvs library and any data
provider implementation.  A provider is responsible for mapping its own
storage schema onto these models; nothing in hgvs outside a provider
should depend on the provider's native representation.

All coordinates are interbase (0-based, half-open).  The ``_i`` suffix
marks an interbase coordinate.

"""

import warnings
from dataclasses import dataclass
from typing import Any


class _DeprecatedItemAccess:
    """Transitional dict-style access for models that replaced dict results.

    Data providers returned plain dicts before these models existed.
    ``model["field"]`` and ``model.get("field")`` still work, but emit a
    :class:`DeprecationWarning`; use attribute access instead.  This shim
    will be removed in a future release.
    """

    __slots__ = ()

    def __getitem__(self, key: str) -> Any:
        warnings.warn(
            f"{type(self).__name__}[{key!r}] is deprecated; use .{key} instead",
            DeprecationWarning,
            stacklevel=2,
        )
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(key) from None

    def get(self, key: str, default: Any = None) -> Any:
        warnings.warn(
            f"{type(self).__name__}.get({key!r}) is deprecated; use .{key} instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(self, key, default)


@dataclass(frozen=True, slots=True)
class TranscriptInfo(_DeprecatedItemAccess):
    """A transcript as aligned to a reference sequence, with its CDS bounds.

    :param tx_ac: transcript accession with version (e.g., ``NM_000051.3``)
    :param alt_ac: accession of the reference sequence the transcript is
        aligned to (e.g., ``NC_000011.9``)
    :param alt_aln_method: alignment method (e.g., ``splign``, ``blat``)
    :param cds_start_i: CDS start on the transcript, or ``None`` for a
        non-coding transcript
    :param cds_end_i: CDS end on the transcript, or ``None`` for a
        non-coding transcript
    :param hgnc: HGNC gene symbol, or ``None`` if not known
    """

    tx_ac: str
    alt_ac: str
    alt_aln_method: str
    cds_start_i: int | None
    cds_end_i: int | None
    hgnc: str | None


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
