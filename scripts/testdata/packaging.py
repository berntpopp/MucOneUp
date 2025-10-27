"""Tarball packaging for test datasets."""

import logging
import tarfile
from pathlib import Path

from .metadata import MetadataGenerator


class Packager:
    """Creates compressed tarballs for distribution."""

    def __init__(self, version: str):
        """
        Initialize packager.

        Args:
            version: Version string for naming
        """
        self.version = version
        self.logger = logging.getLogger(__name__)

    def create_tarball(self, dataset_dir: Path, output_name: str | None = None) -> Path:
        """
        Create compressed tarball of dataset.

        Args:
            dataset_dir: Directory to package
            output_name: Optional custom output name (default: auto-generated)

        Returns:
            Path to created tarball
        """
        if output_name is None:
            output_name = f"testdata_40-70_{self.version}.tar.gz"

        tarball_path = dataset_dir.parent / output_name

        self.logger.info(f"Creating tarball: {output_name}")

        # Create tarball with maximum compression
        with tarfile.open(tarball_path, "w:gz", compresslevel=9) as tar:
            tar.add(dataset_dir, arcname=dataset_dir.name)

        # Calculate tarball checksum
        checksum = MetadataGenerator._calculate_checksum(tarball_path)

        # Write checksum file
        checksum_file = tarball_path.with_suffix(".tar.gz.sha256")
        with open(checksum_file, "w") as f:
            f.write(f"{checksum}  {output_name}\n")

        size_mb = tarball_path.stat().st_size / 1024 / 1024
        self.logger.info(f"✓ Tarball created: {size_mb:.2f} MB")
        self.logger.info(f"✓ Tarball SHA256: {checksum}")

        return tarball_path
