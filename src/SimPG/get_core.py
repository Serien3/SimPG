"""Get the core sequence nodes of the crowd"""

import pickle
from collections import Counter
from .classes import Minigfa
from . import logger
from typing import Optional
import os
import time

__all__ = ["get_coreSeg_in_Pangenome"]


def _get_coreSeg_in_Pangenome(
    gfa_message: Minigfa, every_sample_Whole_Genome_Sequencing_filepath: str
) -> set[tuple[str, str]]:

    core_segs_only: set[tuple[str, str]] = set()
    """
    if (
        bed_message is not None
        and every_sample_Whole_Genome_Sequencing_filepath is None
        and gfa_message is None
    ):
        for _, _, _, _, list_of_segments in bed_message:
            core_segs_only.add((list_of_segments[0], "+"))
            core_segs_only.add((list_of_segments[-1], "+"))
        return core_segs_only
    """
    starttime = time.time()
    dic_path = {}
    with open(every_sample_Whole_Genome_Sequencing_filepath, "rb") as f:
        while True:
            try:
                key, path_list = pickle.load(f)
                if path_list is None:
                    logger.warning(
                        f"Find a path_list is None.This should be because there is No target_SR for {key}"
                    )
                    continue
            except EOFError:
                break
            dic_path[key] = path_list
    counter = Counter()
    idx = 1
    for path in dic_path.values():
        if path is None:
            continue

        segs = [p for p in path]
        for s in segs:
            counter[s] += 1
        # print(f"finish{idx}")
        idx += 1
    core_segs_only = {s for s, c in counter.items() if c == len(dic_path)}
    out: set[tuple[str, str]] = set()
    linear_sample = gfa_message.get_linear_reference()
    for segID, orient in core_segs_only:
        if orient == "+" and gfa_message.get_source_sample(segID) == linear_sample:
            out.add((segID, orient))
    logger.info(
        f"Finish getting the core sequence nodes in {(time.time() - starttime):.2f} seconds.There are {len(out)} nodes in the pan-genome that are core sequence nodes"
    )
    return out


def _save_something(something, s: str) -> None:
    """Save something in pickle format"""
    with open(s, "wb") as f:
        pickle.dump(something, f)


def _save_to_tmp(content, filename: str) -> str:
    """
    Create (or open) a file in the tmp folder under the working directory and write the content.
    """

    # Get the current working directory
    cwd = os.getcwd()
    # Construct the path to the tmp directory
    tmp_dir = os.path.join(cwd, "tmp")
    # If the tmp directory does not exist, create it
    os.makedirs(tmp_dir, exist_ok=True)
    # Construct the full path of the file
    file_path = os.path.join(tmp_dir, filename)
    # Open the file and write the content (if the file does not exist, create a new one, if it exists, overwrite it)
    with open(file_path, "wb") as f:
        pickle.dump(content, f)
    return file_path


def get_coreSeg_in_Pangenome(
    gfa_message: Minigfa,
    every_sample_Whole_Genome_Sequencing_filepath: Optional[str] = None,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> set[tuple[str, str]]:
    """
    Get the core sequence nodes belonging to a certain group of people
    Notice:

    Args:
        gfa_message (Minigfa): Composite data storing GFA file information.
        every_sample_Whole_Genome_Sequencing_filepath (str | None, optional): The file location of the walking route of each sample. The default is the my_walks.pl file in the tmp folder of the working directory
        is_saved_as_pickle (bool, optional): Whether to save as a pickle file for reuse. Defaults to False.
        file_path (str | None, optional): If you choose to save as a pickle file,the graph will be saved in `file_path`. By default, the file name will be `myCoreseg.pl` in folder /tmp under your working folder.


    Returns:
        set[tuple[str, str]]:Save the collection of core sequence nodes
    """
    if every_sample_Whole_Genome_Sequencing_filepath is None:
        every_sample_Whole_Genome_Sequencing_filepath = os.path.join(
            os.getcwd(), "tmp", "my_walks.pl"
        )
    out_set = _get_coreSeg_in_Pangenome(
        gfa_message=gfa_message,
        every_sample_Whole_Genome_Sequencing_filepath=every_sample_Whole_Genome_Sequencing_filepath,
    )
    if is_saved_as_pickle:
        if file_path is None:
            _save_to_tmp(out_set, "myCoreseg.pl")
        else:
            _save_something(out_set, file_path)
    return out_set


if __name__ == "__main__":
    pass
