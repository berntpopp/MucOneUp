# muc_one_up/config.py

import json
import os
import logging


def load_config(config_path):
    """
    Load and validate the JSON config for MucOneUp.

    :param config_path: Path to the JSON config file.
    :return: Python dict with configuration data.
    :raises FileNotFoundError: if the config file does not exist.
    """
    if not os.path.exists(config_path):
        logging.error("Config file not found: %s", config_path)
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as fh:
        config = json.load(fh)

    logging.debug("Configuration loaded from %s", config_path)
    # TODO: Add any validation you want here (e.g. check for required keys)
    return config
