"""Tests for :mod:`pymzlib.isobaric`, the projection of mzLib's ``IsobaricMassTag``.

Offline except for the live canary at the end. ``isobaric_kits.json`` (every kit) and
``isobaric_kits_TMT18.json`` (``--kit TMT18``) are recorded from the real bridge. The m/z values are
mzLib's; what is pinned here is that the table and the per-kit grouping agree, in order, with the
labels on the right channels.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pymzlib import _bridge, isobaric

FIXTURES = Path(__file__).parent / "fixtures"


def _load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))


@pytest.fixture()
def replay(monkeypatch):
    seen: dict = {}

    def use(name: str):
        def fake_invoke(*args, **kwargs):
            seen["args"] = args
            return _load(name)

        monkeypatch.setattr(_bridge, "invoke", fake_invoke)
        return seen

    return use


def test_every_kit_groups_the_table(replay):
    seen = replay("isobaric_kits.json")
    result = isobaric.kits()

    assert seen["args"] == ("isobaric", "kits")
    assert result.kit is None
    names = [k.name for k in result.kits]
    assert names == ["TMT6", "TMT10", "TMT11", "TMT16", "TMT18", "iTRAQ4", "iTRAQ8", "diLeu4", "diLeu12"]
    assert result.kit_count == len(names)
    assert result.record_count == sum(k.channel_count for k in result.kits) == len(result.channels)
    for kit in result.kits:
        assert len(kit.channels) == kit.channel_count
        assert [c.index for c in kit.channels] == list(range(kit.channel_count))
        assert kit.reporter_ion_mzs == sorted(kit.reporter_ion_mzs)
    assert result.absolute_tolerance == 0.003
    assert len(result.caveats) == 5


def test_tmtpro_16_is_the_first_sixteen_of_18(replay):
    replay("isobaric_kits.json")
    result = isobaric.kits()
    tmt16 = result.kit_named("TMT16").reporter_ion_mzs
    tmt18 = result.kit_named("TMT18").reporter_ion_mzs
    assert tmt16 == tmt18[:16]


def test_itraq8_has_no_120_channel(replay):
    replay("isobaric_kits.json")
    assert isobaric.kits().kit_named("iTRAQ8").labels == [
        "113", "114", "115", "116", "117", "118", "119", "121"
    ]


def test_one_kit_and_its_window(replay):
    seen = replay("isobaric_kits_TMT18.json")
    result = isobaric.kits("TMT18")

    assert seen["args"] == ("isobaric", "kits", "--kit", "TMT18")
    assert result.kit == "TMT18"
    (kit,) = result.kits
    assert kit.labels[:3] == ["126", "127N", "127C"] and kit.labels[-1] == "135N"
    for c in kit.channels:
        assert c.kit == "TMT18"
        assert c.mz_min == pytest.approx(c.reporter_ion_mz - result.absolute_tolerance, abs=1e-9)
        assert c.mz_max == pytest.approx(c.reporter_ion_mz + result.absolute_tolerance, abs=1e-9)


def test_kit_named_refuses_a_kit_not_listed(replay):
    replay("isobaric_kits_TMT18.json")
    with pytest.raises(_bridge.UsageError, match="TMT6"):
        isobaric.kits("TMT18").kit_named("TMT6")


@pytest.mark.parametrize("kit", ["", "   ", 10])
def test_a_blank_or_non_string_kit_is_refused_before_the_bridge(monkeypatch, kit):
    monkeypatch.setattr(_bridge, "invoke", lambda *a, **k: pytest.fail("the bridge was called"))
    with pytest.raises(_bridge.UsageError):
        isobaric.kits(kit)


# ---- live: the real bridge ----------------------------------------------------------------------


@pytest.fixture()
def built_bridge():
    # As test_proteins_live: skip in a source checkout with no built bridge, fail under CI.
    try:
        _bridge.bridge_path()
    except _bridge.BridgeNotFoundError as exc:
        if os.environ.get("CI"):
            pytest.fail(f"the mzLib bridge is not built under CI, so isobaric kits went untested: {exc}")
        pytest.skip(f"bridge binary not built: {exc}")


def test_live_kits_match_the_recording_and_refuse_a_partial_name(built_bridge):
    assert isobaric.kits().columns == _load("isobaric_kits.json")["columns"]
    assert isobaric.kits("iTRAQ-4plex on K").kits[0].name == "iTRAQ4"
    with pytest.raises(_bridge.UsageError, match="TMT10plex"):
        isobaric.kits("TMT10plex")
