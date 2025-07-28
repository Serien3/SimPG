import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


from .SimPG import *

__all__ = [
    "run_SimPG",
    "simulate_Whole_Genome_Sequencing_for_population",
    "set_default_logging",
]
