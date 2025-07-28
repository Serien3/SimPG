"""Assemble a pan-genome for a specific population"""

import pickle
import networkx as nx
from .classes import Minibed
from . import logger
from typing import Optional
import os
import time


def _save_graph(graph, s: str) -> None:
    """Save a graph in pickle format"""
    with open(s, "wb") as f:
        pickle.dump(graph, f)


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


def _generate_sequence(start_str, end_str) -> list[tuple[str, str]]:
    """Creates a contiguous sequence of nodes from start_str to end_str

    Returns:
        list[tuple[str, str]]
    """
    prefix = start_str[0]

    # Extract the numeric part and convert to integer
    start_num = int(start_str[1:])
    end_num = int(end_str[1:])

    if start_num > end_num:
        raise ValueError("The starting number cannot be greater than the ending number")

    # Use range to generate a sequence of [start_num, end_num] and then piece together the prefix
    seq = [f"{prefix}{i}" for i in range(start_num, end_num + 1)]
    seq = [(x, "+") for x in seq]
    return seq


def simulate_Population_Pangenome(
    bed_message: Minibed,
    every_sample_Whole_Genome_Sequencing_filepath: Optional[str] = None,
    is_added_linear_reference_genome: bool = False,
    is_output_inspection_results_in_graph: bool = False,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> nx.DiGraph:
    """Simulate the pan-genome of a specific population

    Args:
        bed_message (Minibed): Composite data storing Bed file information.
        every_sample_Whole_Genome_Sequencing_filepath (str |None, optional): The file location of the walking route of each sample.The default is the `my_walks.pl` file in the tmp folder of the working directory
        is_added_linear_reference_genome(bool,optional): Whether to add a linear reference genome in new pan-genome graph.Defaults to False.
        is_output_inspection_results (bool, optional): Whether to output the key parameters of the graph to stdout. Defaults to False.
        is_saved_as_pickle (bool, optional): Whether to save as a pickle file for reuse. Defaults to False.
        file_path (str | None, optional): If you choose to save as a pickle file,the graph will be saved in `file_path`. By default, the file name will be `myPangenome.pl` in folder /tmp under your working folder.

    Returns:
        nx.DiGraph: Pan-genome graph
    """
    if every_sample_Whole_Genome_Sequencing_filepath is None:
        every_sample_Whole_Genome_Sequencing_filepath = os.path.join(
            os.getcwd(), "tmp", "my_walks.pl"
        )
    Pangenome_DiGraph = nx.DiGraph()
    sources, sinks = bed_message.get_linear_sources_and_sinks()
    # print(sources)
    # print(sources.values())
    starttime = time.time()
    with open(every_sample_Whole_Genome_Sequencing_filepath, "rb") as f:
        while True:
            try:
                key, path_list = pickle.load(f)
                if path_list is None:
                    logger.warning(
                        f"Find a path_list is None.This should be because there is No target_SR for {key}.Jump out"
                    )
                    continue
                for u, v in zip(path_list[:-1], path_list[1:]):
                    if v[0] in sources.values():
                        continue
                    Pangenome_DiGraph.add_edge(u, v)
            except EOFError:
                break
    if is_added_linear_reference_genome:
        for chr in sources.keys():
            start_segID = sources[chr]
            end_segID = sinks[chr]
            linear_list = _generate_sequence(start_segID, end_segID)
            for node_from, node_to in zip(linear_list[:-1], linear_list[1:]):
                Pangenome_DiGraph.add_edge(node_from, node_to)
    if is_saved_as_pickle:
        if file_path is None:
            _save_to_tmp(Pangenome_DiGraph, "myPangenome.pl")
        else:
            _save_graph(Pangenome_DiGraph, file_path)
    logger.info(
        "Finish simulate population pangenome in %0.2f seconds."
        % (time.time() - starttime)
    )
    if is_output_inspection_results_in_graph:
        count1 = nx.number_weakly_connected_components(Pangenome_DiGraph)
        print("=================Key information of new graph=================")
        print("The number of weakly connected components of the graph: ", count1)
        print("Number of nodes in the graph: ", nx.number_of_nodes(Pangenome_DiGraph))
        print("Number of edges in the graph: ", nx.number_of_edges(Pangenome_DiGraph))
        sources = [
            n for n in Pangenome_DiGraph.nodes() if Pangenome_DiGraph.in_degree(n) == 0
        ]
        sinks = [
            n for n in Pangenome_DiGraph.nodes() if Pangenome_DiGraph.out_degree(n) == 0
        ]
        print(
            "The number of nodes with zero in-degree(chromosome starting nodes) in the graph: ",
            len(sources),
        )
        print(
            "The number of nodes with zero in-degree(chromosome starting nodes) in the graph: ",
            len(sinks),
        )
        print("nodes with zero in-degree are:", end=" ")
        for x in sources:
            print(x, end=" ")
        print("\n")
        print("nodes with zero out-degree are:", end=" ")
        for y in sinks:
            print(y, end=" ")
        print("\n")
        print("===========================================================")
    return Pangenome_DiGraph


if __name__ == "__main__":
    pass
