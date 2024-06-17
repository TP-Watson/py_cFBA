"""Logging for the Python cFBA Toolbox."""

__all__ = []

import logging


def get_logger() -> logging.Logger:
    """Return a logger with the specified name. The logger writes messages to the console.

    Returns:
        Logger
    """
    logging.captureWarnings(True)
    logger = logging.getLogger("py_cFBA")

    if not logger.hasHandlers():
        logger.setLevel(logging.DEBUG)

        formatter = logging.Formatter("%(message)s")
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        stream_handler.setLevel(logging.INFO)

        logger.addHandler(stream_handler)

    return logger


logger = get_logger()
