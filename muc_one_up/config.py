import json
import os

def load_config(config_path):
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :return: Python dict with configuration data.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as fh:
        config = json.load(fh)

    # TODO: Add any validation you want here, e.g. checking for keys, etc.
    # e.g., if "repeats" not in config: raise ValueError("No 'repeats' section in config")

    return config
