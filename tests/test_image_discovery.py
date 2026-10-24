"""Image assets are identified by source location and layout, not modality names."""

import json
from pathlib import Path

import pytest
from aind_zarr_utils.neuroglancer import get_image_sources

from aind_ibl_ephys_alignment_preprocessing.discovery import image_asset_path


@pytest.mark.parametrize("fixture", ["nested_image_sources.json", "transformed_image_sources.json"])
@pytest.mark.parametrize("prefix", ["SmartSPIM_", ""])
def test_real_image_sources_resolve_to_the_asset_mount(tmp_path, fixture, prefix):
    text = (Path(__file__).parent / "fixtures" / fixture).read_text()
    state = json.loads(text.replace("SmartSPIM_", prefix))
    sources = get_image_sources(state, remove_zarr_protocol=True)
    for uri in sources.values():
        name = uri.split("/image_tile_fusing/")[0].rsplit("/", 1)[1]
        assert image_asset_path(tmp_path, uri) == tmp_path / name


def test_local_image_source_preserves_its_existing_asset_root(tmp_path):
    asset = tmp_path / "elsewhere" / "my-brain"
    asset.mkdir(parents=True)
    source = asset / "image_tile_fusing" / "OMEZarr" / "channel.ome.zarr" / "0"
    assert image_asset_path(tmp_path / "data", str(source)) == asset


def test_bad_image_source_reports_the_uri(tmp_path):
    with pytest.raises(ValueError, match="s3://bucket/unrecognised/channel.zarr"):
        image_asset_path(tmp_path, "s3://bucket/unrecognised/channel.zarr")
