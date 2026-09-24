"""Golden md5 of seeded legacy amplicon FASTQs with real pbsim3/ccs/samtools.

The 0.45 compatibility promise: simulations without a read profile produce the
same reads. The md5s below were measured with origin/main eabb721 (MucOneUp
0.45.0) and are unchanged by #132. Exact procedure (repository ``config.json``
with the two ``model_file`` entries made absolute and ``threads`` set to 4 in
``pacbio_params`` and ``ont_amplicon_params``)::

    muconeup --config config.json --log-level WARNING simulate --out-base dupc \\
        --seed 42 --fixed-lengths 40 --mutation-name dupC --mutation-targets 1,20
    muconeup --config legacy_config.json --log-level WARNING reads amplicon \\
        dupc.001.simulated.fa --platform ont --coverage 100 --seed 42 --no-align \\
        --out-dir out_ont --out-base legacy        # legacy_amplicon_ont.fastq
    # same with --platform pacbio                  # legacy_amplicon_hifi.fastq

Measured with pbsim3 and ccs from the ``env_pacbio`` conda environment and
samtools from ``env_wessim``. Tool locations: ``MUCONEUP_PBSIM3``,
``MUCONEUP_CCS``, ``MUCONEUP_SAMTOOLS`` (default: ``pbsim``/``pbsim3``, ``ccs``,
``samtools`` on PATH) and ``MUCONEUP_PBSIM3_MODELS`` (directory with
``QSHMM-ONT-HQ.model`` and ``ERRHMM-SEQUEL.model``; default
``reference/pbsim3``). The test skips with a reason when any is missing.

Run: ``pytest -m integration tests/e2e/test_legacy_amplicon_golden.py --no-cov``
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli

REPO = Path(__file__).resolve().parents[2]
GOLDEN_MD5 = {
    "ont": ("legacy_amplicon_ont.fastq", "7ec5c90da9afe8ff0c1cda30d682b899"),
    "pacbio": ("legacy_amplicon_hifi.fastq", "221a5a49c9d7d0e197ff4214199d7c2b"),
}
MODELS = {"ont_amplicon_params": "QSHMM-ONT-HQ.model", "pacbio_params": "ERRHMM-SEQUEL.model"}
THREADS = 4


def _tool(env: str, *names: str) -> str:
    value = os.environ.get(env)
    if value:
        return value
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    pytest.skip(f"{names[0]} not found: put it on PATH or set {env}")


def _config(tmp_path: Path) -> Path:
    models = Path(os.environ.get("MUCONEUP_PBSIM3_MODELS", REPO / "reference" / "pbsim3"))
    config = json.loads((REPO / "config.json").read_text())
    for section, model in MODELS.items():
        if not (models / model).exists():
            pytest.skip(f"pbsim3 model {model} not found in {models}: set MUCONEUP_PBSIM3_MODELS")
        config[section]["model_file"] = str(models / model)
        config[section]["threads"] = THREADS
    config["tools"].update(
        pbsim3=_tool("MUCONEUP_PBSIM3", "pbsim", "pbsim3"),
        ccs=_tool("MUCONEUP_CCS", "ccs"),
        samtools=_tool("MUCONEUP_SAMTOOLS", "samtools"),
    )
    path = tmp_path / "legacy_config.json"
    path.write_text(json.dumps(config))
    return path


def _invoke(args: list[str]) -> None:
    result = CliRunner().invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output


@pytest.mark.integration
@pytest.mark.e2e
@pytest.mark.parametrize("platform", ["ont", "pacbio"])
def test_seeded_legacy_amplicon_fastq_matches_golden_md5(tmp_path: Path, platform: str) -> None:
    config = _config(tmp_path)
    _invoke(
        [
            "--config",
            str(REPO / "config.json"),
            "--log-level",
            "WARNING",
            "simulate",
            "--out-dir",
            str(tmp_path),
            "--out-base",
            "dupc",
            "--seed",
            "42",
            "--fixed-lengths",
            "40",
            "--mutation-name",
            "dupC",
            "--mutation-targets",
            "1,20",
        ]
    )
    out_dir = tmp_path / f"out_{platform}"
    _invoke(
        [
            "--config",
            str(config),
            "--log-level",
            "WARNING",
            "reads",
            "amplicon",
            str(tmp_path / "dupc.001.simulated.fa"),
            "--platform",
            platform,
            "--coverage",
            "100",
            "--seed",
            "42",
            "--no-align",
            "--out-dir",
            str(out_dir),
            "--out-base",
            "legacy",
        ]
    )
    name, expected = GOLDEN_MD5[platform]
    assert hashlib.md5((out_dir / name).read_bytes()).hexdigest() == expected
