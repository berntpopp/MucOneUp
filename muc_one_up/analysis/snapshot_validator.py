"""
SNaPshot Assay Validation Module

Provides in-silico validation of SNaPshot assay for MUC1 VNTR mutations.
Simulates complete workflow: PCR → Digest → SNaPshot Extension → Detection.

Classes:
    PCRSimulator: PCR amplification with multi-product handling
    DigestSimulator: Restriction enzyme digest selection
    SnapshotExtensionSimulator: Single-base extension simulation
    SnapshotValidator: Workflow orchestrator

Design Principles:
    - Config-based (no hardcoded values)
    - SOLID principles (SRP, DIP, OCP, LSP, ISP)
    - DRY (single source of truth in config)
    - Modular (independently testable components)
"""

import logging
from typing import Any

from Bio.Restriction import Restriction
from Bio.Seq import Seq


class PCRSimulator:
    """
    Simulates PCR amplification with multiple product handling.

    Handles non-specific primers that bind to multiple sites, producing
    multiple amplicons of different sizes. Applies size filtering and
    mutation-based prioritization.

    Attributes:
        forward_primer (str): Forward primer sequence
        reverse_primer (str): Reverse primer sequence (reverse complemented if needed)
        max_products (int): Maximum products to return
        size_min (int): Minimum amplicon size
        size_max (int): Maximum amplicon size
    """

    def __init__(self, config: dict[str, Any]):
        """
        Initialize PCR simulator with configuration.

        Args:
            config: PCR configuration containing primers and parameters
        """
        pcr_config = config["pcr"]

        self.forward_primer = pcr_config["forward_primer"]
        self.reverse_primer = pcr_config["reverse_primer"]
        self.max_products = pcr_config["max_products"]
        self.size_min = pcr_config["size_range"]["min"]
        self.size_max = pcr_config["size_range"]["max"]

        # Apply reverse complement if specified
        if pcr_config.get("reverse_needs_rc", False):
            self.reverse_primer = str(Seq(self.reverse_primer).reverse_complement())

        self.logger = logging.getLogger(__name__)

    def amplify(self, template: str, mutation_patterns: dict[str, str]) -> dict[str, Any]:
        """
        Simulate PCR amplification with size filtering and mutation prioritization.

        Finds all primer binding sites, generates all possible amplicons,
        filters by size range, and prioritizes by mutation type.

        Args:
            template: Template DNA sequence
            mutation_patterns: Dict with "mutant" and "normal" patterns

        Returns:
            Dict containing:
                - success (bool): Whether PCR succeeded
                - product_count (int): Number of products
                - products (List[Dict]): Product details
                - non_specific (bool): Whether primers are non-specific
        """
        # Find all primer binding sites
        forward_sites = self._find_binding_sites(template, self.forward_primer, reverse=False)
        reverse_sites = self._find_binding_sites(template, self.reverse_primer, reverse=True)

        if not forward_sites or not reverse_sites:
            self.logger.warning("No primer binding sites found")
            return {
                "success": False,
                "product_count": 0,
                "products": [],
                "non_specific": False,
            }

        # Generate all possible amplicons
        products_mutant = []
        products_normal = []
        products_other = []

        mutant_pattern = mutation_patterns.get("mutant", "")
        normal_pattern = mutation_patterns.get("normal", "")

        for f_pos in forward_sites:
            for r_pos in reverse_sites:
                if f_pos < r_pos:
                    # Extract amplicon
                    amplicon_seq = template[f_pos : r_pos + len(self.reverse_primer)]

                    # Filter by size
                    if self.size_min <= len(amplicon_seq) <= self.size_max:
                        # Extract flanked region (between primers)
                        flanked_start = f_pos + len(self.forward_primer)
                        flanked_end = r_pos
                        flanked_region = template[flanked_start:flanked_end]

                        product = {
                            "sequence": amplicon_seq,
                            "length": len(amplicon_seq),
                            "forward_pos": f_pos,
                            "reverse_pos": r_pos,
                            "flanked_region": flanked_region,
                            "mutation_type": "unknown",
                        }

                        # Categorize by mutation type
                        if mutant_pattern and mutant_pattern in amplicon_seq:
                            product["mutation_type"] = "mutant"
                            products_mutant.append(product)
                        elif normal_pattern and normal_pattern in amplicon_seq:
                            product["mutation_type"] = "normal"
                            products_normal.append(product)
                        else:
                            products_other.append(product)

        # Prioritize: mutant first, then normal, then unknown
        all_products = products_mutant + products_normal + products_other

        # Limit to max products
        final_products = all_products[: self.max_products]

        non_specific = len(forward_sites) > 1 or len(reverse_sites) > 1

        self.logger.info(
            f"PCR: {len(forward_sites)} x {len(reverse_sites)} = "
            f"{len(forward_sites) * len(reverse_sites)} potential products, "
            f"{len(final_products)} after size filtering ({self.size_min}-{self.size_max}bp)"
        )

        return {
            "success": True,
            "product_count": len(final_products),
            "products": final_products,
            "non_specific": non_specific,
        }

    def _find_binding_sites(self, template: str, primer: str, reverse: bool = False) -> list[int]:
        """
        Find all binding sites for a primer in template.

        Args:
            template: Template sequence
            primer: Primer sequence (already RC'd for reverse primer in __init__)
            reverse: If True, indicates this is the reverse primer (already RC'd)

        Returns:
            List of binding site positions
        """
        sites = []
        start = 0

        # Primer is already in the correct orientation (forward or RC'd)
        # Just search for it directly in template
        while True:
            pos = template.find(primer, start)
            if pos == -1:
                break
            sites.append(pos)
            start = pos + 1

        return sites


class DigestSimulator:
    """
    Simulates restriction enzyme digest selection.

    Identifies restriction enzyme recognition sites in amplicons and
    determines which amplicons survive digest (no sites) vs. get
    destroyed (have sites).

    Attributes:
        enzyme_name (str): Restriction enzyme name (e.g., "MwoI")
        enzyme: Bio.Restriction enzyme object
        recognition_site (str): Recognition site pattern (for reference)
    """

    def __init__(self, config: dict[str, Any]):
        """
        Initialize digest simulator with configuration.

        Args:
            config: Digest configuration containing enzyme info
        """
        digest_config = config["digest"]

        self.enzyme_name = digest_config["enzyme"]
        self.recognition_site = digest_config.get("recognition_site", "")

        # Import enzyme from Bio.Restriction
        try:
            self.enzyme = getattr(Restriction, self.enzyme_name)
        except AttributeError as err:
            raise ValueError(
                f"Unknown restriction enzyme: {self.enzyme_name}. "
                f"Available in Bio.Restriction: {dir(Restriction)}"
            ) from err

        self.logger = logging.getLogger(__name__)

    def digest(self, amplicons: list[dict[str, Any]]) -> dict[str, Any]:
        """
        Simulate restriction digest - products with sites are destroyed.

        Args:
            amplicons: List of PCR products from PCRSimulator

        Returns:
            Dict containing:
                - total_products (int): Total input amplicons
                - digested_count (int): Number digested
                - survivor_count (int): Number that survived
                - survivors (List[Dict]): Surviving amplicons
                - digested (List[Dict]): Digested amplicons
                - mechanism (str): Description of mechanism
        """
        survivors = []
        digested = []

        for amp in amplicons:
            seq_obj = Seq(amp["sequence"])
            sites = self.enzyme.search(seq_obj)

            if len(sites) == 0:
                # No restriction site → survives digest
                survivors.append(
                    {
                        **amp,
                        "has_restriction_site": False,
                        "site_count": 0,
                        "survives_digest": True,
                    }
                )
            else:
                # Has restriction site → gets digested
                digested.append(
                    {
                        **amp,
                        "has_restriction_site": True,
                        "site_count": len(sites),
                        "site_positions": sites,
                        "survives_digest": False,
                    }
                )

        self.logger.info(
            f"Digest ({self.enzyme_name}): {len(survivors)} survivors, " f"{len(digested)} digested"
        )

        return {
            "total_products": len(amplicons),
            "digested_count": len(digested),
            "survivor_count": len(survivors),
            "survivors": survivors,
            "digested": digested,
            "mechanism": f"{self.enzyme_name} digest selection",
        }


class SnapshotExtensionSimulator:
    """
    Simulates SNaPshot single-base extension.

    Models primer annealing to amplicon and single-base extension with
    fluorescent ddNTPs, determining which fluorescent signal is produced.

    Attributes:
        primers (Dict[str, str]): SNaPshot primer sequences
        fluorophore_map (Dict[str, Dict]): Base to fluorophore mapping
    """

    def __init__(self, config: dict[str, Any]):
        """
        Initialize SNaPshot extension simulator with configuration.

        Args:
            config: SNaPshot configuration containing primers and fluorophore map
        """
        snapshot_config = config["snapshot"]

        self.primers = snapshot_config["primers"]
        self.fluorophore_map = snapshot_config["fluorophore_map"]

        self.logger = logging.getLogger(__name__)

    def extend(self, template: str, primer_name: str) -> dict[str, Any]:
        """
        Simulate SNaPshot single-base extension.

        Args:
            template: Amplicon sequence (from digest survivors)
            primer_name: Which primer to use (key from config primers dict)

        Returns:
            Dict containing:
                - binds (bool): Whether primer binds
                - position (int): Binding position (-1 if no binding)
                - next_base (str): Base to be incorporated (A, C, G, T)
                - ddNTP (str): ddNTP that gets incorporated
                - fluorophore (str): Fluorophore color
                - fluorophore_dye (str): Fluorophore dye name
                - peak_size (int): Expected peak size (primer length + 1)
                - interpretation (str): Result interpretation
        """
        if primer_name not in self.primers:
            raise ValueError(
                f"Unknown primer '{primer_name}'. Available: {list(self.primers.keys())}"
            )

        primer_seq = self.primers[primer_name]

        # Find primer binding site
        pos = template.find(primer_seq)
        if pos == -1:
            self.logger.debug(f"Primer {primer_name} does not bind to template")
            return {"binds": False, "position": -1, "interpretation": "No binding"}

        # Get next base after primer
        extension_pos = pos + len(primer_seq)
        if extension_pos >= len(template):
            return {
                "binds": True,
                "position": pos,
                "next_base": None,
                "interpretation": "Primer binds at template end (no extension)",
            }

        next_base = template[extension_pos]
        fluor_info = self.fluorophore_map.get(next_base, {})

        # Interpretation based on base
        interpretations = {
            "C": "8C mutation detected (Black peak)",
            "A": "Incomplete extension or alternate mutation",
            "G": "Unexpected - check amplicon sequence",
            "T": "Unexpected - check amplicon sequence",
        }

        result = {
            "binds": True,
            "position": pos,
            "next_base": next_base,
            "ddNTP": f"dd{next_base}TP",
            "fluorophore": fluor_info.get("color", "Unknown"),
            "fluorophore_dye": fluor_info.get("dye", "Unknown"),
            "peak_size": len(primer_seq) + 1,
            "interpretation": interpretations.get(next_base, f"Unexpected base: {next_base}"),
        }

        self.logger.debug(
            f"Extension: {primer_name} at pos {pos}, next base: {next_base} "
            f"({result['fluorophore']})"
        )

        return result


class SnapshotValidator:
    """
    Orchestrates complete SNaPshot validation workflow.

    Coordinates PCR, digest, and SNaPshot extension to simulate complete
    assay and determine mutation presence.

    Attributes:
        mutation_name (str): Mutation being validated (e.g., "dupC")
        mutation_config (Dict): Complete config for this mutation
        pcr (PCRSimulator): PCR simulator instance
        digest (DigestSimulator): Digest simulator instance
        extension (SnapshotExtensionSimulator): Extension simulator instance
        mutation_patterns (Dict): Mutant and normal patterns
    """

    def __init__(self, config: dict[str, Any], mutation_name: str):
        """
        Initialize validator with configuration.

        Args:
            config: Complete config dictionary
            mutation_name: Mutation to validate (e.g., "dupC")

        Raises:
            ValueError: If config is missing required sections
        """
        # Validate config structure
        if "snapshot_validation" not in config:
            raise ValueError("Config missing 'snapshot_validation' section")

        if mutation_name not in config["snapshot_validation"]:
            raise ValueError(
                f"Mutation '{mutation_name}' not found in config. "
                f"Available: {list(config['snapshot_validation'].keys())}"
            )

        self.mutation_name = mutation_name
        self.mutation_config = config["snapshot_validation"][mutation_name]

        # Initialize components with config (Dependency Injection)
        self.pcr = PCRSimulator(self.mutation_config)
        self.digest = DigestSimulator(self.mutation_config)
        self.extension = SnapshotExtensionSimulator(self.mutation_config)

        # Mutation patterns for validation
        self.mutation_patterns = self.mutation_config["validation"]

        self.logger = logging.getLogger(__name__)

    def validate_complete_workflow(self, template_seq: str) -> dict[str, Any]:
        """
        Simulate complete PCR → Digest → SNaPshot workflow.

        Args:
            template_seq: Genomic template sequence

        Returns:
            Dict containing:
                - pcr_results: PCR simulation results
                - digest_results: Digest simulation results
                - snapshot_results: SNaPshot extension results
                - mutation_detected: Whether mutation was detected
                - expected_fluorescence: Expected fluorescence color
                - summary: Human-readable summary
        """
        self.logger.info(f"Starting SNaPshot validation for {self.mutation_name}")

        # Step 1: PCR Amplification
        pcr_results = self.pcr.amplify(
            template=template_seq,
            mutation_patterns={
                "mutant": self.mutation_patterns["mutant_pattern"],
                "normal": self.mutation_patterns["normal_pattern"],
            },
        )

        if not pcr_results["success"]:
            return {
                "pcr_results": pcr_results,
                "mutation_detected": False,
                "summary": "PCR failed - no products generated",
            }

        self.logger.info(f"PCR: {pcr_results['product_count']} products")

        # Step 2: Restriction Digest Selection
        digest_results = self.digest.digest(pcr_results["products"])

        self.logger.info(
            f"Digest: {digest_results['survivor_count']} survivors, "
            f"{digest_results['digested_count']} digested"
        )

        # Step 3: SNaPshot Extension (on survivors only)
        snapshot_results = []
        for survivor in digest_results["survivors"]:
            # Try primary primer (primer_7c)
            extension_result = self.extension.extend(
                template=survivor["sequence"], primer_name="primer_7c"
            )
            snapshot_results.append({"amplicon": survivor, "extension": extension_result})

        # Step 4: Determine mutation presence
        # Mutation detected if any survivor shows "C" (Black fluorescence)
        mutation_detected = any(
            r["extension"].get("next_base") == "C"
            for r in snapshot_results
            if r["extension"].get("binds", False)
        )

        expected_fluor = self.mutation_patterns.get("expected_mutant_fluorescence", "Black")

        return {
            "pcr_results": pcr_results,
            "digest_results": digest_results,
            "snapshot_results": snapshot_results,
            "mutation_detected": mutation_detected,
            "expected_fluorescence": expected_fluor,
            "summary": self._generate_summary(
                pcr_results, digest_results, snapshot_results, mutation_detected
            ),
        }

    def validate_amplicon(self, amplicon_seq: str) -> dict[str, Any]:
        """
        Validate if specific amplicon survives digest.

        Useful for testing individual amplicons without full workflow.

        Args:
            amplicon_seq: Single amplicon sequence

        Returns:
            Dict containing:
                - mutation_present: Whether mutant pattern present
                - has_normal: Whether normal pattern present
                - has_restriction_site: Whether enzyme site present
                - site_count: Number of restriction sites
                - survives_digest: Whether amplicon survives
                - amplicon_length: Length in bp
                - detectable: Whether mutation is detectable
        """
        # Check mutation patterns
        mutant_pattern = self.mutation_patterns["mutant_pattern"]
        normal_pattern = self.mutation_patterns["normal_pattern"]

        has_mutant = mutant_pattern in amplicon_seq
        has_normal = normal_pattern in amplicon_seq

        # Check restriction sites
        seq_obj = Seq(amplicon_seq)
        sites = self.digest.enzyme.search(seq_obj)

        return {
            "mutation_present": has_mutant,
            "has_normal": has_normal,
            "has_restriction_site": len(sites) > 0,
            "site_count": len(sites),
            "site_positions": sites if sites else [],
            "survives_digest": len(sites) == 0,
            "amplicon_length": len(amplicon_seq),
            "detectable": has_mutant and len(sites) == 0,
        }

    def _generate_summary(
        self,
        pcr_results: dict[str, Any],
        digest_results: dict[str, Any],
        snapshot_results: list[dict[str, Any]],
        mutation_detected: bool,
    ) -> str:
        """
        Generate human-readable summary of validation results.

        Args:
            pcr_results: PCR simulation results
            digest_results: Digest simulation results
            snapshot_results: SNaPshot extension results
            mutation_detected: Whether mutation was detected

        Returns:
            Human-readable summary string
        """
        if mutation_detected:
            fluorescence_details = []
            for result in snapshot_results:
                if result["extension"].get("binds", False):
                    ext = result["extension"]
                    fluorescence_details.append(
                        f"{ext.get('fluorophore', 'Unknown')} "
                        f"({ext.get('fluorophore_dye', 'Unknown')})"
                    )

            fluor_str = (
                ", ".join(fluorescence_details) if fluorescence_details else "No fluorescence"
            )

            return (
                f"{self.mutation_name} mutation DETECTED: "
                f"{digest_results['survivor_count']} survivor(s), "
                f"fluorescence: {fluor_str}"
            )
        else:
            if digest_results["survivor_count"] == 0:
                return (
                    f"No {self.mutation_name} mutation detected: "
                    f"all products digested (normal phenotype)"
                )
            else:
                return (
                    f"No {self.mutation_name} mutation detected: "
                    f"{digest_results['survivor_count']} survivor(s) but no mutant signal"
                )
