# muc_one_up/distribution.py

import logging
import random

from .type_defs import LengthModelDict


def sample_repeat_count(length_model: LengthModelDict) -> int:
    """Sample repeat count from configured distribution.

    Samples from a normal distribution truncated between min_repeats and max_repeats.
    Falls back to mean_repeats if distribution type is not supported.

    Args:
        length_model: Dictionary with keys 'distribution', 'min_repeats', 'max_repeats',
                     'mean_repeats', 'median_repeats'

    Returns:
        Integer representing the sampled repeat count

    Example:
        >>> model = {
        ...     "distribution": "normal",
        ...     "min_repeats": 20,
        ...     "max_repeats": 100,
        ...     "mean_repeats": 60,
        ...     "median_repeats": 60
        ... }
        >>> count = sample_repeat_count(model)
        >>> 20 <= count <= 100
        True
    """
    dist_type = length_model.get("distribution", "normal")
    min_rep = length_model["min_repeats"]
    max_rep = length_model["max_repeats"]
    mean_rep = length_model["mean_repeats"]

    if dist_type == "normal":
        while True:
            val = int(random.gauss(mean_rep, (max_rep - min_rep) / 4.0))
            if min_rep <= val <= max_rep:
                logging.debug("Sampled repeat count: %d", val)
                return val
    else:
        logging.warning("Distribution type '%s' not implemented; using mean_repeats.", dist_type)
        return int(mean_rep)  # Ensure integer return
