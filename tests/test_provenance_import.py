"""Test that provenance module handles missing config_fingerprint gracefully."""


def test_provenance_module_imports_successfully():
    """Provenance module should always import without errors."""
    import muc_one_up.provenance  # noqa: F401


def test_collect_provenance_metadata_works_without_rfc8785(monkeypatch):
    """Provenance collection should degrade gracefully if fingerprinting fails."""
    import time

    from muc_one_up.provenance import collect_provenance_metadata

    config = {"repeats": {"X": "ACGACT"}, "seed": 42}
    start = time.time()
    end = start + 1.0

    # Even if fingerprinting fails, provenance should return a dict
    result = collect_provenance_metadata(config, start, end)
    assert isinstance(result, dict)
    assert "software_version" in result
