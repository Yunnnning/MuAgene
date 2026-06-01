"""Shared constants/helpers for SnapATAC2 spectral embedding in the ATAC path."""
from __future__ import annotations

from typing import Any

import numpy as np

# Primary latent from snap.tl.spectral (SnapATAC2 default for knn / umap / leiden).
ATAC_LATENT_KEY = "X_spectral"
# Backward-compat alias kept in S5+ exports for Signac-oriented consumers.
ATAC_LATENT_ALIAS = "X_lsi"


def get_atac_latent(obsm: dict[str, Any]) -> np.ndarray | None:
    """Return the ATAC clustering latent, preferring X_spectral over legacy X_lsi."""
    if ATAC_LATENT_KEY in obsm:
        return np.asarray(obsm[ATAC_LATENT_KEY])
    if ATAC_LATENT_ALIAS in obsm:
        return np.asarray(obsm[ATAC_LATENT_ALIAS])
    return None
