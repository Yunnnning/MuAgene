from types import SimpleNamespace
from unittest import mock

import pandas as pd

from executor.qc_summary import _final_section
from executor.run_paths import RunPaths


def test_final_section_reads_canonical_paired_deliverable(tmp_path):
    paths = RunPaths(tmp_path)
    paths.ensure()
    paths.deliv_results.mkdir(parents=True, exist_ok=True)
    paths.processed_h5mu.write_text("placeholder")
    adata = SimpleNamespace(
        n_obs=1,
        obs=pd.DataFrame(index=["cell-1"]),
        obs_names=pd.Index(["cell-1"]),
    )
    mdata = SimpleNamespace(mod={"rna": adata, "atac": adata})

    with mock.patch("mudata.read_h5mu", return_value=mdata) as read_h5mu:
        section = _final_section(
            tmp_path,
            {"rna_final": 1, "atac_final": 1, "rna_post_doublet": 1, "atac_post_doublet": 1},
        )

    read_h5mu.assert_called_once_with(str(paths.processed_h5mu))
    assert "paired" in section
