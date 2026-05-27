"""Detect whether RNA and ATAC modalities are paired multiome or separate datasets."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any


SUFFIX_PATTERNS = [
    re.compile(r"-\d+$"),        # e.g. AAACGAA-1
    re.compile(r"_[A-Z0-9]+$"),  # e.g. AAACGAA_LIBRARY
]


def _normalize(bc: str) -> list[str]:
    """Return candidate normalized forms of a barcode (original + stripped variants)."""
    variants = [bc]
    for pat in SUFFIX_PATTERNS:
        m = pat.search(bc)
        if m:
            variants.append(bc[: m.start()])
    return variants


@dataclass
class PairingResult:
    status: str           # paired | separate | ambiguous | rna_only | atac_only
    confidence: str       # high | medium | low
    method: str           # pairing.{exact_barcode_match|prefix_suffix_normalized|
                          #           translation_table|single_file_multiome|
                          #           rna_only_input|atac_only_input|ambiguous_overlap|no_match}
    overlap: float        # |intersection|/|union| (0.0 for single-modality inputs)
    n_rna: int
    n_atac: int
    n_shared: int
    normalization: str | None = None
    assumptions: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "confidence": self.confidence,
            "method": self.method,
            "overlap": round(self.overlap, 6),
            "n_rna": self.n_rna,
            "n_atac": self.n_atac,
            "n_shared": self.n_shared,
            "normalization": self.normalization,
            "assumptions": self.assumptions,
        }


def detect_pairing(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    *,
    single_file_multiome: bool = False,
) -> PairingResult:
    n_rna, n_atac = len(rna_barcodes), len(atac_barcodes)

    # Single-modality inputs — declared by absence of the other barcode set.
    # This is a valid user-declared workflow, not a failure.
    if n_rna == 0 and n_atac == 0:
        raise ValueError("detect_pairing: both RNA and ATAC barcode sets are empty")
    if n_atac == 0:
        return PairingResult(
            status="rna_only",
            confidence="high",
            method="pairing.rna_only_input",
            overlap=0.0,
            n_rna=n_rna,
            n_atac=0,
            n_shared=0,
            assumptions=["Only RNA input provided; no ATAC modality to pair against."],
        )
    if n_rna == 0:
        return PairingResult(
            status="atac_only",
            confidence="high",
            method="pairing.atac_only_input",
            overlap=0.0,
            n_rna=0,
            n_atac=n_atac,
            n_shared=0,
            assumptions=["Only ATAC input provided; no RNA modality to pair against."],
        )

    if single_file_multiome:
        # Both modalities share a Cell Ranger ARC .h5 -> paired by construction.
        shared = rna_barcodes & atac_barcodes
        return PairingResult(
            status="paired",
            confidence="high",
            method="pairing.single_file_multiome",
            overlap=1.0,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=len(shared),
            assumptions=["RNA and ATAC came from the same Cell Ranger ARC .h5"],
        )

    # Strategy 2: exact match
    shared = rna_barcodes & atac_barcodes
    union = rna_barcodes | atac_barcodes
    overlap = len(shared) / max(len(union), 1)
    if overlap >= 0.99:
        return PairingResult(
            status="paired",
            confidence="high",
            method="pairing.exact_barcode_match",
            overlap=overlap,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=len(shared),
        )

    # Strategy 3: prefix/suffix normalization
    rna_norm = {v for bc in rna_barcodes for v in _normalize(bc)}
    atac_norm = {v for bc in atac_barcodes for v in _normalize(bc)}
    shared_n = rna_norm & atac_norm
    union_n = rna_norm | atac_norm
    overlap_n = len(shared_n) / max(len(union_n), 1)
    if overlap_n >= 0.99:
        return PairingResult(
            status="paired",
            confidence="medium",
            method="pairing.prefix_suffix_normalized",
            overlap=overlap_n,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=len(shared_n),
            normalization="strip -N / _LIBRARY suffixes",
        )

    # Intermediate -> ambiguous
    if 0.30 <= overlap < 0.99 or 0.30 <= overlap_n < 0.99:
        return PairingResult(
            status="ambiguous",
            confidence="low",
            method="pairing.ambiguous_overlap",
            overlap=max(overlap, overlap_n),
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=len(shared) if overlap >= overlap_n else len(shared_n),
        )

    # Low overlap -> separate datasets (valid branch, not a failure)
    return PairingResult(
        status="separate",
        confidence="high",
        method="pairing.no_match",
        overlap=overlap,
        n_rna=n_rna,
        n_atac=n_atac,
        n_shared=len(shared),
    )


def pairing_via_translation(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    translation: dict[str, str],
) -> PairingResult:
    """Re-check pairing after rewriting ATAC barcodes into RNA-space via a translation table.

    `translation` maps ATAC barcode -> RNA barcode (one-to-one cell pairing).
    Unmapped ATAC barcodes are excluded from the comparison; the unmapped
    count is recorded in `assumptions` so S0 can record coverage in
    `validation_report.json`.

    Status is `paired` when the *translated* overlap reaches ≥0.99 of the
    union (matching the threshold used by `pairing.exact_barcode_match`);
    otherwise `separate` is returned with the actual overlap. Callers
    decide whether to fall through to the next ladder rung.
    """
    n_rna = len(rna_barcodes)
    n_atac_raw = len(atac_barcodes)
    translated: set[str] = set()
    n_unmapped = 0
    for bc in atac_barcodes:
        rna = translation.get(bc)
        if rna is None:
            n_unmapped += 1
        else:
            translated.add(rna)
    shared = rna_barcodes & translated
    union = rna_barcodes | translated
    overlap = len(shared) / max(len(union), 1)
    status = "paired" if overlap >= 0.99 else "separate"
    return PairingResult(
        status=status,
        confidence="high" if status == "paired" else "medium",
        method="pairing.translation_table",
        overlap=overlap,
        n_rna=n_rna,
        n_atac=n_atac_raw,
        n_shared=len(shared),
        normalization=f"atac->rna via user-supplied translation table ({len(translation)} pairs)",
        assumptions=[
            f"{n_unmapped} of {n_atac_raw} ATAC barcodes had no translation entry "
            "(dropped from the comparison)."
        ],
    )
