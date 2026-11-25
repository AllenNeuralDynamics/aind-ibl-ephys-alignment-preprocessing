"""Inventory transport and discovery/Flatten/worker/pack contracts."""

import io
import json
import shutil
from pathlib import Path
from unittest.mock import Mock

from aind_ibl_ephys_alignment_preprocessing import stages
from aind_ibl_ephys_alignment_preprocessing.ephys_inventory import (
    S3MetadataStore,
    carry_inventory,
    discover_inventory,
    write_unit_inventory,
)
from aind_ibl_ephys_alignment_preprocessing.types import PipelineConfig

STREAM = "Record Node 103#Neuropix-PXI-111.ProbeA"


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def _config(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    manifest = data / "manifest.csv"
    manifest.write_text(
        "mouseid,sorted_recording,probe_file,histology_track_id,ephys_collection\n123,rec_sorted,tracks,trackA,ProbeA\n"
    )
    return PipelineConfig(
        data_root=data, results_root=tmp_path / "results", manifest_csv=manifest, neuroglancer_file=data / "tracks.json"
    )


def _sources(config):
    raw = config.data_root / "raw"
    name = f"experiment7_{STREAM}.zarr"
    _write(
        raw / "ecephys_compressed" / name / ".zattrs",
        {
            "num_segments": 1,
            "sampling_frequency": 30000,
            "annotations": {"experiment_name": "experiment7"},
        },
    )
    _write(raw / "ecephys_compressed" / name / "traces_seg0/.zarray", {"shape": [30000, 2]})
    _write(
        raw / "ecephys_clipped/Record Node 103/experiment7/recording3/structure.oebin",
        {
            "continuous": [{"folder_name": "Neuropix-PXI-111.ProbeA/", "num_channels": 2, "sample_rate": 30000}],
        },
    )
    sort = config.data_root / "sorted"
    _write(
        sort / f"postprocessed/experiment7_{STREAM}_recording3.zarr/.zmetadata",
        {
            "metadata": {
                "recording_info/.zattrs": {
                    "recording_attributes": {
                        "sampling_frequency": 30000,
                        "num_samples": [30000],
                        "num_channels": 2,
                    }
                },
                "sorting/.zattrs": {"num_segments": 1},
                "sorting/unit_ids/.zarray": {"shape": [3]},
            }
        },
    )
    return [
        {"key": "raw-id", "id": "raw-id", "name": "raw", "role": "raw", "local_mount": "raw"},
        {
            "key": "sorted-id",
            "id": "sorted-id",
            "name": "renamed-sort",
            "aliases": ["rec_sorted"],
            "role": "sorted",
            "local_mount": "sorted",
            "source_keys": ["raw-id"],
        },
    ]


def test_discover_flatten_worker_and_pack_retain_inventory(tmp_path, monkeypatch):
    config = _config(tmp_path)
    descriptors = _sources(config)
    monkeypatch.setattr(stages, "_probe_viability", lambda *args: (True, None, []))
    paths = stages.stage_discover(config, inventory_sources=descriptors)
    assert len(paths) == 1
    canonical = json.loads((config.results_root / "ephys_inventory.json").read_text())
    assert canonical["metadata_complete"] is True
    assert len(canonical["recordings"]) == 1

    # Flatten receives only this directory; no discovery root is available.
    flattened = tmp_path / "flattened"
    shutil.copytree(paths[0].parent, flattened)
    cfg = stages.find_ephys_stream_config(flattened)
    assert cfg["inventory"] == "ephys_inventory.json"
    sliced = json.loads(Path(cfg["_inventory_path"]).read_text())
    assert len(sliced["recordings"]) == 1
    assert sliced["scope"]["sorted_recording"] == "rec_sorted"
    worker = config.model_copy(update={"results_root": tmp_path / "worker"})
    extract = Mock()
    monkeypatch.setattr(stages, "run_ephys_for_stream", extract)
    stages.stage_ephys(worker, stream_config=cfg)
    extract.assert_called_once()
    assert (worker.results_root / "rec__ProbeA/ephys_inventory.json").is_file()

    # Pack reads canonical discovery observations, before intermediate cleanup.
    pack = config.model_copy(
        update={"manifest_csv": config.results_root / "manifest.csv", "results_root": tmp_path / "pack"}
    )
    carry_inventory(pack)
    assert json.loads((pack.results_root / "123/ephys_inventory.json").read_text()) == canonical


def test_rejected_rows_do_not_erase_inventory(tmp_path, monkeypatch):
    import pytest

    from aind_ibl_ephys_alignment_preprocessing.coverage import CoverageError

    config = _config(tmp_path)
    descriptors = _sources(config)
    monkeypatch.setattr(stages, "_probe_viability", lambda *args: (False, "missing track", []))
    with pytest.raises(CoverageError):
        stages.stage_discover(config, inventory_sources=descriptors)
    observed = json.loads((config.results_root / "ephys_inventory.json").read_text())
    assert len(observed["sortings"]) == 1
    assert len(observed["recordings"]) == 1


def test_old_discovery_outputs_remain_usable(tmp_path):
    config = _config(tmp_path)
    carry_inventory(config)
    assert not (config.results_root / "123/ephys_inventory.json").exists()


def test_unmounted_source_is_incomplete_not_empty_success(tmp_path):
    config = _config(tmp_path)
    inventory = discover_inventory(config, [], [{"key": "raw-id", "role": "raw"}])
    assert not inventory["metadata_complete"]
    assert inventory["findings"][0]["code"] == "source_unavailable"


def test_s3_transport_uses_bounded_reads_and_one_level_listing(monkeypatch):
    client = Mock()
    monkeypatch.setattr("aind_ibl_ephys_alignment_preprocessing.ephys_inventory.boto3.client", lambda *a, **k: client)
    client.get_paginator.return_value.paginate.return_value = [{"CommonPrefixes": [{"Prefix": "nested/raw/ecephys/"}]}]
    client.get_object.return_value = {"Body": io.BytesIO(b'{"num_segments": 1}')}
    store = S3MetadataStore("s3://bucket/nested/raw", anonymous=True)
    assert store.directories(".") == ["ecephys"]
    client.get_paginator.return_value.paginate.assert_called_once_with(
        Bucket="bucket",
        Prefix="nested/raw/",
        Delimiter="/",
    )
    assert store.read_json(".zattrs") == {"num_segments": 1}
    assert store.read_json(".zattrs") == {"num_segments": 1}
    client.get_object.assert_called_once_with(Bucket="bucket", Key="nested/raw/.zattrs")


def test_s3_timestamp_read_uses_exact_byte_range(monkeypatch):
    client = Mock()
    monkeypatch.setattr("aind_ibl_ephys_alignment_preprocessing.ephys_inventory.boto3.client", lambda *a, **k: client)
    client.get_object.return_value = {"Body": io.BytesIO(b"12345678")}
    store = S3MetadataStore("s3://bucket/nested/raw")
    assert store.read_bytes("recording.zarr/times_seg0/0", offset=16, size=8) == b"12345678"
    client.get_object.assert_called_once_with(
        Bucket="bucket", Key="nested/raw/recording.zarr/times_seg0/0", Range="bytes=16-23"
    )
    client.get_object.return_value = {"Body": io.BytesIO(b"range ignored")}
    import pytest

    with pytest.raises(ValueError, match="requested byte range"):
        store.read_bytes("recording.zarr/times_seg0/0", offset=16, size=8)


def test_s3_compressed_timestamp_read_is_bounded(monkeypatch):
    from aind_ibl_ephys_alignment_preprocessing.ephys_inventory import MAX_METADATA_BYTES

    client = Mock()
    monkeypatch.setattr("aind_ibl_ephys_alignment_preprocessing.ephys_inventory.boto3.client", lambda *a, **k: client)
    body = Mock()
    body.__enter__ = Mock(return_value=body)
    body.__exit__ = Mock(return_value=False)
    body.read.return_value = b"compressed timestamp bytes"
    client.get_object.return_value = {"Body": body}
    store = S3MetadataStore("s3://bucket/raw")
    assert store.read_chunk("recording.zarr/times_seg0/0") == b"compressed timestamp bytes"
    body.read.assert_called_once_with(MAX_METADATA_BYTES + 1)


def test_slicing_never_takes_other_sessions(tmp_path):
    config = _config(tmp_path)
    inventory = discover_inventory(config, [], _sources(config))
    other = dict(inventory["recordings"][0], id="other-rec", asset="other-raw")
    inventory["recordings"].append(other)
    inventory["assets"].append({"key": "other-raw", "name": "other", "role": "raw", "metadata_status": "complete"})
    path = tmp_path / "unit/config.json"
    _write(path, {"ephys_collection": "ProbeA", "sorted_recording": "rec_sorted"})
    write_unit_inventory(inventory, path)
    sliced = json.loads((path.parent / "ephys_inventory.json").read_text())
    assert len(sliced["recordings"]) == 1
    assert all(a["key"] != "other-raw" for a in sliced["assets"])


def test_per_sort_launcher_keeps_the_inventory_slice(tmp_path):
    config = _config(tmp_path)
    sources = _sources(config)
    inventory = discover_inventory(config, [], sources)
    _write(config.data_root / "ephys_inventory.json", inventory)
    (config.data_root / "sorted/spikesorted").mkdir()
    _write(config.data_root / "sorted/data_description.json", {"source_data": ["rec"]})
    paths = stages.stage_ephys_launch(config)
    assert len(paths) == 1
    cfg = json.loads(paths[0].read_text())
    assert cfg["inventory"] == "ephys_inventory.json"
    sliced = json.loads((paths[0].parent / cfg["inventory"]).read_text())
    assert len(sliced["sortings"]) == 1
