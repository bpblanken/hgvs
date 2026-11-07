"""Data models returned by HGVS data providers.

These models are the contract between the hgvs library and any data
provider implementation.  A provider is responsible for mapping its own
storage schema onto these models; nothing in hgvs outside a provider
should depend on the provider's native representation.

All coordinates are interbase (0-based, half-open).  The ``_i`` suffix
marks an interbase coordinate.

"""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TranscriptInfo:
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


@dataclass(frozen=True, slots=True)
class TranscriptIdentityInfo:
    """A transcript described on its own sequence: CDS bounds and exon lengths.

    :param tx_ac: transcript accession with version
    :param alt_ac: reference the transcript is defined on; for a transcript
        described on itself this is ``tx_ac``
    :param alt_aln_method: alignment method; ``transcript`` for the identity
        alignment
    :param cds_start_i: CDS start on the transcript, or ``None`` if non-coding
    :param cds_end_i: CDS end on the transcript, or ``None`` if non-coding
    :param lengths: exon lengths in transcript order; their sum is the
        transcript length
    :param hgnc: HGNC gene symbol, or ``None`` if not known
    """

    tx_ac: str
    alt_ac: str
    alt_aln_method: str
    cds_start_i: int | None
    cds_end_i: int | None
    lengths: tuple[int, ...]
    hgnc: str | None


@dataclass(frozen=True, slots=True)
class TranscriptExon:
    """One exon of a transcript aligned to a reference sequence.

    :param tx_ac: transcript accession with version
    :param alt_ac: reference sequence accession
    :param alt_aln_method: alignment method
    :param alt_strand: strand of the alignment on the reference, ``1`` or ``-1``
    :param ord: exon ordinal in transcript order, starting at 0
    :param tx_start_i: exon start on the transcript
    :param tx_end_i: exon end on the transcript
    :param alt_start_i: exon start on the reference
    :param alt_end_i: exon end on the reference
    :param cigar: CIGAR string of the exon alignment, transcript to reference
    :param hgnc: HGNC gene symbol, or ``None`` if not known
    """

    tx_ac: str
    alt_ac: str
    alt_aln_method: str
    alt_strand: int
    ord: int
    tx_start_i: int
    tx_end_i: int
    alt_start_i: int
    alt_end_i: int
    cigar: str
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
