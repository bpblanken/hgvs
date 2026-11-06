"""Data models returned by HGVS data providers.

These models are the contract between the hgvs library and any data
provider implementation.  A provider is responsible for mapping its own
storage schema onto these models; nothing in hgvs outside a provider
should depend on the provider's native representation.

All coordinates are interbase (0-based, half-open).  The ``_i`` suffix
marks an interbase coordinate.

"""

import datetime
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


@dataclass(frozen=True, slots=True)
class TranscriptMappingOption:
    """A reference sequence and method a transcript has an alignment for.

    :param tx_ac: transcript accession with version
    :param alt_ac: reference sequence accession the transcript is aligned to
    :param alt_aln_method: alignment method (e.g., ``splign``, ``blat``)
    """

    tx_ac: str
    alt_ac: str
    alt_aln_method: str


@dataclass(frozen=True, slots=True)
class TranscriptAlignment:
    """The span a transcript's alignment covers on a reference sequence.

    :param tx_ac: transcript accession with version
    :param alt_ac: reference sequence accession the transcript is aligned to
    :param alt_aln_method: alignment method (e.g., ``splign``, ``blat``)
    :param alt_strand: strand of the alignment on the reference, ``1`` or ``-1``
    :param start_i: alignment start on the reference
    :param end_i: alignment end on the reference
    """

    tx_ac: str
    alt_ac: str
    alt_aln_method: str
    alt_strand: int
    start_i: int
    end_i: int


@dataclass(frozen=True, slots=True)
class SimilarTranscript:
    """How two transcripts compare structurally.

    ``es`` is exon structure, ``fp`` fingerprint, ``eq`` equal.

    :param tx_ac1: the transcript the comparison was made for
    :param tx_ac2: the transcript it is compared to
    :param hgnc_eq: same gene symbol
    :param cds_eq: identical CDS sequences
    :param es_fp_eq: identical exon structures, UTRs included
    :param cds_es_fp_eq: identical exon structures clipped to the CDS
    :param cds_exon_lengths_fp_eq: identical CDS-clipped exon lengths
    """

    tx_ac1: str
    tx_ac2: str
    hgnc_eq: bool
    cds_eq: bool
    es_fp_eq: bool
    cds_es_fp_eq: bool
    cds_exon_lengths_fp_eq: bool


@dataclass(frozen=True, slots=True)
class GeneInfo:
    """Basic information about a gene.

    :param hgnc: HGNC gene symbol
    :param symbol: display symbol, usually the same as ``hgnc``
    :param gene_id: gene identifier in the provider's source
    :param maploc: cytogenetic map location (e.g., ``11q22-q23``)
    :param descr: description
    :param summary: longer summary
    :param aliases: comma-separated alias symbols
    :param type: gene type, if known
    :param xrefs: cross-references, if known
    :param added: when the record was added to the provider's source
    """

    hgnc: str
    symbol: str | None
    gene_id: str | None
    maploc: str | None
    descr: str | None
    summary: str | None
    aliases: str | None
    type: str | None
    xrefs: str | None
    added: datetime.datetime | None


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
