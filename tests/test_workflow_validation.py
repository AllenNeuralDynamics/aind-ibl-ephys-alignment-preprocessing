"""Packed results must retain both discovery decisions and consumer files."""

import json

import numpy as np
import pandas as pd
import pytest
from aind_ephys_ibl_gui_conversion.conversion_plan import inventory_digest

from aind_ibl_ephys_alignment_preprocessing.workflow_validation import audit_package, record_output_manifests


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def packed(tmp_path):
    root = tmp_path / "858398"
    folder = root / "rec/ProbeA"
    folder.mkdir(parents=True)
    inventory = dict(recordings=[], sortings=[])
    plan = dict(
        inventory_digest="unit-inventory",
        main_segment=dict(duration_s=10),
        analyzers=[dict(unit_count=1)],
        shank_coverage=dict(recorded=[2], sorted=[2], missing=[]),
    )
    _write(root / "ephys_inventory.json", inventory)
    _write(
        root / "ephys_conversion_plan.json",
        dict(inventory_digest=inventory_digest(inventory), units={"rec__ProbeA": plan}),
    )
    _write(
        folder / "ephys_extraction.json",
        dict(
            plan=plan, execution_enabled=True, status="complete", blocks=[dict(selected_main=True)], skipped_blocks=[]
        ),
    )
    _write(root / "track_data/pick.json", [])
    _write(
        root / "datapackage.json",
        dict(
            schema_version="4.1.0",
            external_assets={},
            probes=dict(
                rec=dict(
                    ProbeA=dict(
                        ephys=dict(asset=None, path="rec/ProbeA"),
                        xyz_picks=[dict(image=dict(asset=None, path="track_data/pick.json"))],
                    )
                )
            ),
        ),
    )
    arrays = {
        "spikes.times": [0.1, 0.3],
        "spikes.clusters": np.array([0, 0], dtype=int),
        "spikes.amps": [1, 2],
        "spikes.depths": [20, 20],
        "spike_shank_indices": [0, 0],
        "unit_shank_indices": [0],
        "clusters.channels": np.array([2], dtype=int),
        "channels.localCoordinates": [[0, 30], [0, 10], [0, 20], [0, 0]],
        "channels.shankInd": [2, 2, 2, 2],
        "channels.rawInd": np.arange(4),
        "channels.contactId": ["s2e0", "s2e1", "s2e2", "s2e3"],
        "clusters.peakToTrough": [0.001],
        "clusters.waveforms": np.zeros((1, 3, 1)),
        "clusters.waveformTimes": [-0.001, 0, 0.001],
        "clusters.waveformChannels": np.array([[2]], dtype=int),
        "_iblqc_ephysSpectralDensityLF.power": np.ones((3, 4)),
        "_iblqc_ephysSpectralDensityLF.freqs": [0, 1, 2],
    }
    for band in ("AP", "LF"):
        arrays[f"_iblqc_ephysTimeRms{band}.rms"] = np.ones((2, 4))
        arrays[f"_iblqc_ephysTimeRms{band}.timestamps"] = [0.1, 0.2]
    for name, array in arrays.items():
        np.save(folder / f"{name}.npy", array)
    pd.DataFrame(dict(cluster_id=[0], analyzer_id=["a"], source_unit_id=[7])).to_csv(
        folder / "clusters.metrics.csv", index=False
    )
    _write(
        folder / "band_corr/row_channels.json", dict(matrix_rows=[0, 1, 2, 3], shanks={"2": dict(rows=[3, 1, 2, 0])})
    )
    for band in ("delta", "theta", "alpha", "beta", "gamma"):
        for suffix in ("mean_corr", "coherency"):
            np.save(folder / "band_corr" / f"{band}_{suffix}.npy", np.eye(4))
    (folder / "optional_metric.txt").write_text("optional adapter output")
    assert record_output_manifests(tmp_path) == 1
    return root


def test_complete_package_preserves_decisions_and_sparse_shanks(packed):
    report = audit_package(packed)
    assert report["status"] == "passed"
    assert report["global_inventory_digest_verified"]
    assert report["collections"]["rec__ProbeA"]["unit_shanks"] == [0]
    assert report["collections"]["rec__ProbeA"]["physical_unit_shanks"] == [2]


@pytest.mark.parametrize("damage", ["inventory", "receipt", "qc", "reference", "missing_unit", "optional_file"])
def test_packaging_regressions_are_rejected(packed, damage):
    if damage == "inventory":
        _write(packed / "ephys_inventory.json", dict(changed=True))
    elif damage == "receipt":
        path = packed / "rec/ProbeA/ephys_extraction.json"
        receipt = json.loads(path.read_text())
        receipt["plan"]["main_segment"]["duration_s"] = 9
        _write(path, receipt)
    elif damage == "qc":
        (packed / "rec/ProbeA/band_corr/theta_coherency.npy").unlink()
    elif damage == "optional_file":
        (packed / "rec/ProbeA/optional_metric.txt").unlink()
    else:
        path = packed / "datapackage.json"
        package = json.loads(path.read_text())
        if damage == "reference":
            package["probes"]["rec"]["ProbeA"]["xyz_picks"][0]["image"]["path"] = "../outside"
        else:
            package["probes"] = {}
        _write(path, package)
    with pytest.raises((ValueError, FileNotFoundError)):
        audit_package(packed)
