"""Local DNA read assembly and variant calling API.

POST /assemble with a JSON body {"reads": [...]} returns the shortest
sequence containing every read (lexicographically smallest on ties) and
the earliest 0-based position of each read, in the original input order.

POST /variants with a reference and aligned forward reads returns
confident single-base substitutions as {"variants": [...]}.

POST /align with a reference and 1..20 reads returns the best full-length
semi-global placement of each read or its reverse complement, allowing
substitutions and at most one single-base insertion/deletion.

POST /digest with an A/C/G/T sequence, a linear/circular topology and
1..5 restriction templates returns the de-duplicated cut positions and
the fragments between consecutive cuts.
"""

from __future__ import annotations

from typing import Literal

from fastapi import FastAPI
from pydantic import BaseModel, Field, field_validator, model_validator

from .aligner import align_reads
from .assembler import assemble
from .digester import DigestTemplate, digest
from .variants import AlignedRead, call_variants

MIN_READS = 2
MAX_READS = 12
MIN_READ_LEN = 1
MAX_READ_LEN = 80
ALPHABET = frozenset("ACGT")

MIN_REFERENCE_LEN = 1
MAX_REFERENCE_LEN = 1000
MIN_ALIGNED_READS = 1
MAX_ALIGNED_READS = 100
MIN_QUALITY = 2
MAX_QUALITY = 40

MIN_ALIGN_READS = 1
MAX_ALIGN_READS = 20
MIN_MAX_EDITS = 0
MAX_MAX_EDITS = 5

MIN_DIGEST_SEQUENCE_LEN = 12
MAX_DIGEST_SEQUENCE_LEN = 1000
MIN_DIGEST_TEMPLATES = 1
MAX_DIGEST_TEMPLATES = 5
MIN_PATTERN_LEN = 2
MAX_PATTERN_LEN = 12
PATTERN_ALPHABET = frozenset("ACGTNRY")

app = FastAPI(title="DNA Assembly API", version="1.0.0")


class AssembleRequest(BaseModel):
    reads: list[str] = Field(min_length=MIN_READS, max_length=MAX_READS)

    @field_validator("reads")
    @classmethod
    def _validate_reads(cls, reads: list[str]) -> list[str]:
        for read in reads:
            if not MIN_READ_LEN <= len(read) <= MAX_READ_LEN:
                raise ValueError(
                    f"each read must be {MIN_READ_LEN}-{MAX_READ_LEN} bases long"
                )
            if not ALPHABET.issuperset(read):
                raise ValueError("reads may only contain A, C, G, T")
        return reads


class AssembleResponse(BaseModel):
    sequence: str
    length: int
    positions: list[int]


@app.post("/assemble", response_model=AssembleResponse)
def assemble_endpoint(request: AssembleRequest) -> AssembleResponse:
    sequence, positions = assemble(request.reads)
    return AssembleResponse(
        sequence=sequence, length=len(sequence), positions=positions
    )


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


# ---------------------------------------------------------------------------
# POST /align
# ---------------------------------------------------------------------------


class AlignRequest(BaseModel):
    reference: str = Field(
        min_length=MIN_REFERENCE_LEN, max_length=MAX_REFERENCE_LEN
    )
    reads: list[str] = Field(
        min_length=MIN_ALIGN_READS, max_length=MAX_ALIGN_READS
    )
    max_edits: int = Field(
        ge=MIN_MAX_EDITS, le=MAX_MAX_EDITS, strict=True
    )

    @field_validator("reference")
    @classmethod
    def _validate_reference(cls, reference: str) -> str:
        if not ALPHABET.issuperset(reference):
            raise ValueError("reference may only contain A, C, G, T")
        return reference

    @field_validator("reads")
    @classmethod
    def _validate_reads(cls, reads: list[str]) -> list[str]:
        for index, read in enumerate(reads):
            if not MIN_READ_LEN <= len(read) <= MAX_READ_LEN:
                raise ValueError(
                    f"reads[{index}]: each read must be "
                    f"{MIN_READ_LEN}-{MAX_READ_LEN} bases long"
                )
            if not ALPHABET.issuperset(read):
                raise ValueError(
                    f"reads[{index}]: reads may only contain A, C, G, T"
                )
        return reads


class AlignmentResponse(BaseModel):
    start: int
    end: int
    strand: str
    edits: int
    cigar: str


class AlignResponse(BaseModel):
    alignments: list[AlignmentResponse | None]


@app.post("/align", response_model=AlignResponse)
def align_endpoint(request: AlignRequest) -> AlignResponse:
    aligned = align_reads(request.reference, request.reads, request.max_edits)
    return AlignResponse(
        alignments=[
            None
            if alignment is None
            else AlignmentResponse(
                start=alignment.start,
                end=alignment.end,
                strand=alignment.strand,
                edits=alignment.edits,
                cigar=alignment.cigar,
            )
            for alignment in aligned
        ]
    )


# ---------------------------------------------------------------------------
# POST /variants
# ---------------------------------------------------------------------------


class AlignedReadRequest(BaseModel):
    start: int = Field(ge=0)
    sequence: str = Field(min_length=MIN_READ_LEN, max_length=MAX_READ_LEN)
    qualities: list[int] = Field(
        min_length=MIN_READ_LEN, max_length=MAX_READ_LEN
    )

    @field_validator("sequence")
    @classmethod
    def _validate_sequence(cls, sequence: str) -> str:
        if not ALPHABET.issuperset(sequence):
            raise ValueError("sequence may only contain A, C, G, T")
        return sequence

    @field_validator("qualities")
    @classmethod
    def _validate_qualities(cls, qualities: list[int]) -> list[int]:
        for quality in qualities:
            if not MIN_QUALITY <= quality <= MAX_QUALITY:
                raise ValueError(
                    f"each quality must be between {MIN_QUALITY} "
                    f"and {MAX_QUALITY}"
                )
        return qualities


class VariantsRequest(BaseModel):
    reference: str = Field(
        min_length=MIN_REFERENCE_LEN, max_length=MAX_REFERENCE_LEN
    )
    reads: list[AlignedReadRequest] = Field(
        min_length=MIN_ALIGNED_READS, max_length=MAX_ALIGNED_READS
    )

    @field_validator("reference")
    @classmethod
    def _validate_reference(cls, reference: str) -> str:
        if not ALPHABET.issuperset(reference):
            raise ValueError("reference may only contain A, C, G, T")
        return reference

    @model_validator(mode="after")
    def _validate_reads(self) -> "VariantsRequest":
        for index, read in enumerate(self.reads):
            if len(read.sequence) != len(read.qualities):
                raise ValueError(
                    f"reads[{index}]: sequence and qualities must have "
                    "the same length"
                )
            if read.start + len(read.sequence) > len(self.reference):
                raise ValueError(
                    f"reads[{index}]: read extends past the end of the "
                    "reference"
                )
        return self


class VariantResponse(BaseModel):
    position: int
    ref: str
    alt: str
    depth: int


class VariantsResponse(BaseModel):
    variants: list[VariantResponse]


@app.post("/variants", response_model=VariantsResponse)
def variants_endpoint(request: VariantsRequest) -> VariantsResponse:
    aligned = [
        AlignedRead(
            start=read.start,
            sequence=read.sequence,
            qualities=read.qualities,
        )
        for read in request.reads
    ]
    called = call_variants(request.reference, aligned)
    return VariantsResponse(
        variants=[
            VariantResponse(
                position=variant.position,
                ref=variant.ref,
                alt=variant.alt,
                depth=variant.depth,
            )
            for variant in called
        ]
    )


# ---------------------------------------------------------------------------
# POST /digest
# ---------------------------------------------------------------------------


class TemplateRequest(BaseModel):
    pattern: str = Field(
        min_length=MIN_PATTERN_LEN, max_length=MAX_PATTERN_LEN
    )
    offset: int = Field(ge=0, strict=True)

    @field_validator("pattern")
    @classmethod
    def _validate_pattern(cls, pattern: str) -> str:
        if not PATTERN_ALPHABET.issuperset(pattern):
            raise ValueError(
                "pattern may only contain A, C, G, T, N, R, Y"
            )
        return pattern

    @model_validator(mode="after")
    def _validate_offset(self) -> "TemplateRequest":
        if self.offset > len(self.pattern):
            raise ValueError(
                f"offset must be between 0 and the pattern length "
                f"({len(self.pattern)})"
            )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(
        min_length=MIN_DIGEST_SEQUENCE_LEN,
        max_length=MAX_DIGEST_SEQUENCE_LEN,
    )
    topology: Literal["linear", "circular"]
    templates: list[TemplateRequest] = Field(
        min_length=MIN_DIGEST_TEMPLATES, max_length=MAX_DIGEST_TEMPLATES
    )

    @field_validator("sequence")
    @classmethod
    def _validate_sequence(cls, sequence: str) -> str:
        if not ALPHABET.issuperset(sequence):
            raise ValueError("sequence may only contain A, C, G, T")
        return sequence

    @model_validator(mode="after")
    def _validate_patterns_fit(self) -> "DigestRequest":
        for index, template in enumerate(self.templates):
            if len(template.pattern) > len(self.sequence):
                raise ValueError(
                    f"templates[{index}]: pattern must not be longer "
                    "than the sequence"
                )
        return self


class FragmentResponse(BaseModel):
    start: int
    end: int
    sequence: str


class DigestResponse(BaseModel):
    cuts: list[int]
    fragments: list[FragmentResponse]


@app.post("/digest", response_model=DigestResponse)
def digest_endpoint(request: DigestRequest) -> DigestResponse:
    templates = [
        DigestTemplate(pattern=template.pattern, offset=template.offset)
        for template in request.templates
    ]
    cuts, fragments = digest(request.sequence, request.topology, templates)
    return DigestResponse(
        cuts=cuts,
        fragments=[
            FragmentResponse(
                start=fragment.start,
                end=fragment.end,
                sequence=fragment.sequence,
            )
            for fragment in fragments
        ],
    )
