# muc_one_up/distribution.py
import random
import math

def sample_repeat_count(length_model):
    """
    Example: sample from a normal distribution truncated between min_repeats and max_repeats.
    """
    dist_type = length_model.get("distribution", "normal")
    min_rep = length_model["min_repeats"]
    max_rep = length_model["max_repeats"]
    mean_rep = length_model["mean_repeats"]
    
    if dist_type == "normal":
        # naive approach: keep sampling until we fall into [min_rep, max_rep]
        while True:
            val = int(random.gauss(mean_rep, (max_rep - min_rep)/4.0))
            if val >= min_rep and val <= max_rep:
                return val
    else:
        # implement alternative distributions if needed
        return mean_rep  # fallback
