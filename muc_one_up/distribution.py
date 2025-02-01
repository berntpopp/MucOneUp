# muc_one_up/distribution.py

import random
import math
import logging


def sample_repeat_count(length_model):
    """
    Sample from a normal distribution truncated between min_repeats and max_repeats.

    :param length_model: Dict with keys 'distribution', 'min_repeats', 'max_repeats',
                         'mean_repeats', etc.
    :return: An integer representing the sampled repeat count.
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
        logging.warning("Distribution type '%s' not implemented; using mean_repeats.",
                        dist_type)
        return mean_rep
