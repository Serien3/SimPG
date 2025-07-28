import logging
import sys

__all__ = ["set_default_logging"]


def set_default_logging(verbose: bool = False) -> None:
    logLevel = logging.INFO if verbose else logging.WARNING
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
