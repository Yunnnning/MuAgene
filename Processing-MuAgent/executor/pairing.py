"""Detect whether RNA and ATAC modalities are paired multiome or separate datasets."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

# Jaccard |∩|/|∪| threshold for declaring paired via exact or suffix-normalized overlap.
PAIRING_OVERLAP_THRESHOLD = 0.80
# Fraction of the smaller modality's barcodes that must appear in the other set.
SUBSET_COVERAGE_THRESHOLD = 0.80
# Jaccard overlap band that requires user input when no translation table resolves it.
AMBIGUOUS_OVERLAP_LOW = 0.30

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


def jaccard_overlap(rna_barcodes: set[str], atac_barcodes: set[str]) -> tuple[float, int, int, int]:
    """Return (|∩|/|∪|, n_shared, n_rna, n_atac)."""
    shared = rna_barcodes & atac_barcodes
    union = rna_barcodes | atac_barcodes
    overlap = len(shared) / max(len(union), 1)
    return overlap, len(shared), len(rna_barcodes), len(atac_barcodes)


@dataclass
class SubsetRelation:
    """How two barcode sets relate when one modality calls fewer cells than the other."""

    relation: str  # none | atac_subset_of_rna | rna_subset_of_atac
    coverage: float  # |∩| / |smaller modality| when relation != none
    n_shared: int
    n_rna: int
    n_atac: int
    atac_in_rna_fraction: float  # |∩| / |ATAC|
    rna_in_atac_fraction: float  # |∩| / |RNA|

    def as_dict(self) -> dict[str, Any]:
        return {
            "relation": self.relation,
            "coverage": round(self.coverage, 6),
            "n_shared": self.n_shared,
            "n_rna": self.n_rna,
            "n_atac": self.n_atac,
            "atac_in_rna_fraction": round(self.atac_in_rna_fraction, 6),
            "rna_in_atac_fraction": round(self.rna_in_atac_fraction, 6),
        }


def detect_subset_relation(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    *,
    threshold: float = SUBSET_COVERAGE_THRESHOLD,
) -> SubsetRelation:
    """Classify whether ATAC or RNA barcodes are largely contained in the other set.

    Uses coverage of the *smaller* modality when both pass the threshold (typical multiome
    case: fewer ATAC cells than RNA cells at ingest, same barcodes).
    """
    n_rna, n_atac = len(rna_barcodes), len(atac_barcodes)
    shared = rna_barcodes & atac_barcodes
    n_shared = len(shared)
    atac_frac = n_shared / n_atac if n_atac else 0.0
    rna_frac = n_shared / n_rna if n_rna else 0.0

    atac_is_subset = n_atac > 0 and atac_frac >= threshold
    rna_is_subset = n_rna > 0 and rna_frac >= threshold

    if atac_is_subset and (not rna_is_subset or n_atac <= n_rna):
        return SubsetRelation(
            relation="atac_subset_of_rna",
            coverage=atac_frac,
            n_shared=n_shared,
            n_rna=n_rna,
            n_atac=n_atac,
            atac_in_rna_fraction=atac_frac,
            rna_in_atac_fraction=rna_frac,
        )
    if rna_is_subset:
        return SubsetRelation(
            relation="rna_subset_of_atac",
            coverage=rna_frac,
            n_shared=n_shared,
            n_rna=n_rna,
            n_atac=n_atac,
            atac_in_rna_fraction=atac_frac,
            rna_in_atac_fraction=rna_frac,
        )
    return SubsetRelation(
        relation="none",
        coverage=max(atac_frac, rna_frac),
        n_shared=n_shared,
        n_rna=n_rna,
        n_atac=n_atac,
        atac_in_rna_fraction=atac_frac,
        rna_in_atac_fraction=rna_frac,
    )


def is_atac_subset_of_rna(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    *,
    threshold: float = SUBSET_COVERAGE_THRESHOLD,
) -> bool:
    """True when ≥threshold of ATAC barcodes appear in the RNA barcode set."""
    return detect_subset_relation(rna_barcodes, atac_barcodes, threshold=threshold).relation == "atac_subset_of_rna"


def is_rna_subset_of_atac(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    *,
    threshold: float = SUBSET_COVERAGE_THRESHOLD,
) -> bool:
    """True when ≥threshold of RNA barcodes appear in the ATAC barcode set."""
    return detect_subset_relation(rna_barcodes, atac_barcodes, threshold=threshold).relation == "rna_subset_of_atac"


@dataclass
class PairingResult:
    status: str           # paired | separate | ambiguous | rna_only | atac_only
    confidence: str       # high | medium | low
    method: str           # pairing.{exact_barcode_match|prefix_suffix_normalized|
                          #           atac_subset_of_rna|rna_subset_of_atac|
                          #           translation_table|single_file_multiome|...}
    overlap: float        # |intersection|/|union| (0.0 for single-modality inputs)
    n_rna: int
    n_atac: int
    n_shared: int
    normalization: str | None = None
    subset_relation: str | None = None
    subset_coverage: float | None = None
    assumptions: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
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
        if self.subset_relation is not None:
            out["subset_relation"] = self.subset_relation
        if self.subset_coverage is not None:
            out["subset_coverage"] = round(self.subset_coverage, 6)
        return out


def _paired_from_subset(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    rel: SubsetRelation,
    *,
    normalization: str | None = None,
) -> PairingResult:
    overlap, n_shared, n_rna, n_atac = jaccard_overlap(rna_barcodes, atac_barcodes)
    label = rel.relation.replace("_", " ")
    return PairingResult(
        status="paired",
        confidence="high" if rel.coverage >= 0.95 else "medium",
        method=f"pairing.{rel.relation}",
        overlap=overlap,
        n_rna=n_rna,
        n_atac=n_atac,
        n_shared=n_shared,
        normalization=normalization,
        subset_relation=rel.relation,
        subset_coverage=rel.coverage,
        assumptions=[
            f"{label}: {rel.coverage:.1%} of the relevant modality's barcodes match "
            f"({n_shared} shared; Jaccard overlap={overlap:.3f}).",
        ],
    )


def detect_pairing(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    *,
    single_file_multiome: bool = False,
) -> PairingResult:
    n_rna, n_atac = len(rna_barcodes), len(atac_barcodes)

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
        shared = rna_barcodes & atac_barcodes
        overlap, _, _, _ = jaccard_overlap(rna_barcodes, atac_barcodes)
        rel = detect_subset_relation(rna_barcodes, atac_barcodes)
        assumptions = [
            "RNA and ATAC came from the same Cell Ranger ARC run; "
            "fragments.tsv.gz uses the GEX barcode namespace (all droplets), "
            "while the filtered matrix carries the cell-called subset.",
        ]
        return PairingResult(
            status="paired",
            confidence="high",
            method="pairing.single_file_multiome",
            overlap=overlap,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=len(shared),
            subset_relation=rel.relation if rel.relation != "none" else None,
            subset_coverage=rel.coverage if rel.relation != "none" else None,
            assumptions=assumptions,
        )

    overlap, n_shared, _, _ = jaccard_overlap(rna_barcodes, atac_barcodes)
    if overlap >= PAIRING_OVERLAP_THRESHOLD:
        return PairingResult(
            status="paired",
            confidence="high",
            method="pairing.exact_barcode_match",
            overlap=overlap,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=n_shared,
        )

    rel = detect_subset_relation(rna_barcodes, atac_barcodes)
    if rel.relation != "none":
        return _paired_from_subset(rna_barcodes, atac_barcodes, rel)

    rna_norm = {v for bc in rna_barcodes for v in _normalize(bc)}
    atac_norm = {v for bc in atac_barcodes for v in _normalize(bc)}
    overlap_n, n_shared_n, _, _ = jaccard_overlap(rna_norm, atac_norm)
    if overlap_n >= PAIRING_OVERLAP_THRESHOLD:
        return PairingResult(
            status="paired",
            confidence="medium",
            method="pairing.prefix_suffix_normalized",
            overlap=overlap_n,
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=n_shared_n,
            normalization="strip -N / _LIBRARY suffixes",
        )

    rel_n = detect_subset_relation(rna_norm, atac_norm)
    if rel_n.relation != "none":
        return _paired_from_subset(
            rna_barcodes,
            atac_barcodes,
            rel_n,
            normalization="strip -N / _LIBRARY suffixes (subset check on normalized barcodes)",
        )

    if (
        AMBIGUOUS_OVERLAP_LOW <= overlap < PAIRING_OVERLAP_THRESHOLD
        or AMBIGUOUS_OVERLAP_LOW <= overlap_n < PAIRING_OVERLAP_THRESHOLD
    ):
        return PairingResult(
            status="ambiguous",
            confidence="low",
            method="pairing.ambiguous_overlap",
            overlap=max(overlap, overlap_n),
            n_rna=n_rna,
            n_atac=n_atac,
            n_shared=n_shared if overlap >= overlap_n else n_shared_n,
            assumptions=[
                f"Subset check: {rel.as_dict()} (exact barcodes); "
                f"normalized subset: {rel_n.as_dict()}.",
            ],
        )

    return PairingResult(
        status="separate",
        confidence="high",
        method="pairing.no_match",
        overlap=overlap,
        n_rna=n_rna,
        n_atac=n_atac,
        n_shared=n_shared,
        assumptions=[
            f"Jaccard overlap={overlap:.3f}; subset coverage below {SUBSET_COVERAGE_THRESHOLD:.0%}.",
        ],
    )


def pairing_via_translation(
    rna_barcodes: set[str],
    atac_barcodes: set[str],
    translation: dict[str, str],
) -> PairingResult:
    """Re-check pairing after rewriting ATAC barcodes into RNA-space via a translation table."""
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
    overlap, n_shared, _, _ = jaccard_overlap(rna_barcodes, translated)
    if overlap >= PAIRING_OVERLAP_THRESHOLD:
        return PairingResult(
            status="paired",
            confidence="high",
            method="pairing.translation_table",
            overlap=overlap,
            n_rna=n_rna,
            n_atac=n_atac_raw,
            n_shared=n_shared,
            normalization=f"atac->rna via user-supplied translation table ({len(translation)} pairs)",
            assumptions=[
                f"{n_unmapped} of {n_atac_raw} ATAC barcodes had no translation entry "
                "(dropped from the comparison).",
            ],
        )

    rel = detect_subset_relation(rna_barcodes, translated)
    if rel.relation != "none":
        result = _paired_from_subset(
            rna_barcodes,
            translated,
            rel,
            normalization=f"atac->rna via translation table ({len(translation)} pairs)",
        )
        result.n_atac = n_atac_raw
        result.method = "pairing.translation_table"
        if n_unmapped:
            result.assumptions.append(
                f"{n_unmapped} of {n_atac_raw} ATAC barcodes had no translation entry."
            )
        return result

    return PairingResult(
        status="separate",
        confidence="medium",
        method="pairing.translation_table",
        overlap=overlap,
        n_rna=n_rna,
        n_atac=n_atac_raw,
        n_shared=n_shared,
        normalization=f"atac->rna via user-supplied translation table ({len(translation)} pairs)",
        assumptions=[
            f"{n_unmapped} of {n_atac_raw} ATAC barcodes had no translation entry "
            "(dropped from the comparison).",
            f"Post-translation Jaccard={overlap:.3f}; subset={rel.as_dict()}.",
        ],
    )
