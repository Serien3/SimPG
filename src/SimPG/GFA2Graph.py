"""Parse the GFA file,get the Pangenome graph"""

import pickle
import networkx as nx
from . import logger
from collections import Counter, deque
from typing import Optional
from .classes import Minigfa, Minibed
import os
import time

__all__ = ["turn_GFA_to_DiGraph"]


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


def _turn_GFA_to_DiGraph_complex(
    gfa_message: Minigfa, bed_message: Minibed
) -> nx.DiGraph:
    """turn GFA to DiGraph with bed message

    Args:
        gfa_message (Minigfa): GFA data message
        bed_message (Minibed): Bed data message

    Returns:
        nx.DiGraph:Pan-genome graph structure(from 5' to 3' we think)
    """

    def find_twice(strings):
        # Count the number of times all strings appear
        cnt = Counter(strings)
        # Only keep strings that occur 2 times
        return [s for s, c in cnt.items() if c >= 2]

    def double_pass_one(seg_in_double, seg_not_in):
        if first_pos[seg_not_in] < first_pos[seg_in_double]:
            return 1
        elif first_pos[seg_not_in] > last_pos[seg_in_double]:
            return -1
        else:
            return 0

    def reverse_orient(one_orient):
        return "+" if one_orient == "-" else "-"

    G = nx.DiGraph()
    sources: dict[str, str] = {}
    sinks: dict[str, str] = {}
    sources, sinks = bed_message.get_linear_sources_and_sinks()
    last_key = next(reversed(sinks))
    last_value = sinks[last_key]
    last_linear_int = int(last_value[1:])
    sum = 0
    sum_double: list[str] = []
    all_bedline_list: list[str] = []
    for chr in sinks.keys():
        temp_list = _generate_sequence(sources[chr], sinks[chr])
        for u, v in zip(temp_list[:-1], temp_list[1:]):
            G.add_node(u, SR=0)
            G.add_node(v, SR=0)
            G.add_edge(u, v, SR=0, weight=0)
    last_chr = object()
    for chr, is_inverved, _, _, seg_in_bubble in bed_message:
        cur_chr = chr
        if last_chr != cur_chr:
            all_bedline_list.append(sources[chr])
            last_chr = cur_chr
        all_bedline_list.extend(seg_in_bubble[1:])
        if is_inverved:
            temp = find_twice(seg_in_bubble)
            for segname in temp:
                if int(segname[1:]) > last_linear_int:
                    sum += 1
                sum_double.extend(temp)
                G.add_node((segname, "+"), SR=gfa_message.get_SRank(segname))
                G.add_node((segname, "-"), SR=gfa_message.get_SRank(segname))
    sum_double_set = set(sum_double)
    all_bedline_set = set(all_bedline_list)
    first_pos = {}
    last_pos = {}
    for idx, x in enumerate(all_bedline_list):
        # If it appears for the first time, record first_pos
        first_pos.setdefault(x, idx)
        # Update last_pos each time, and the final one is the "last time"
        last_pos[x] = idx

    num = 0
    cache_link = deque()
    for from_id, from_orient, to_id, to_orient, SRi in gfa_message.get_all_Link():
        num += 1
        if num % 200000 == 0:
            logger.debug(f"Already finish {num} Links")
        if from_id not in all_bedline_set or to_id not in all_bedline_set:
            continue
        if from_id in sum_double_set and to_id in sum_double_set:
            G.add_node((from_id, from_orient), SR=gfa_message.get_SRank(from_id))
            G.add_node((to_id, to_orient), SR=gfa_message.get_SRank(to_id))
            G.add_node(
                (from_id, reverse_orient(from_orient)),
                SR=gfa_message.get_SRank(from_id),
            )
            G.add_node(
                (to_id, reverse_orient(to_orient)), SR=gfa_message.get_SRank(to_id)
            )
            G.add_edge(
                (from_id, from_orient),
                (to_id, to_orient),
                SR=SRi,
                weight=0,
            )
            G.add_edge(
                (to_id, reverse_orient(to_orient)),
                (from_id, reverse_orient(from_orient)),
                SR=SRi,
                weight=0,
            )
        elif from_id not in sum_double_set and to_id not in sum_double_set:
            if int(from_id[1:]) <= last_linear_int:
                if from_orient == "+":
                    G.add_node((to_id, to_orient), SR=gfa_message.get_SRank(to_id))
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                if from_orient == "-":
                    G.add_node(
                        (to_id, reverse_orient(to_orient)),
                        SR=gfa_message.get_SRank(to_id),
                    )
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
            elif int(to_id[1:]) <= last_linear_int:
                if to_orient == "+":
                    G.add_node(
                        (from_id, from_orient), SR=gfa_message.get_SRank(from_id)
                    )
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                if to_orient == "-":
                    G.add_node(
                        (from_id, reverse_orient(from_orient)),
                        SR=gfa_message.get_SRank(from_id),
                    )
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
            else:
                if (from_id, from_orient) in G or (to_id, to_orient) in G:
                    G.add_node(
                        (from_id, from_orient), SR=gfa_message.get_SRank(from_id)
                    )
                    G.add_node((to_id, to_orient), SR=gfa_message.get_SRank(to_id))
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                elif (from_id, reverse_orient(from_orient)) in G or (
                    to_id,
                    reverse_orient(to_orient),
                ) in G:
                    G.add_node(
                        (from_id, reverse_orient(from_orient)),
                        SR=gfa_message.get_SRank(from_id),
                    )
                    G.add_node(
                        (to_id, reverse_orient(to_orient)),
                        SR=gfa_message.get_SRank(to_id),
                    )
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
                else:
                    # cache_link.append((from_id, from_orient, to_id, to_orient, SRi, 0))
                    # logging.info(f"error!{from_id} to {to_id}")
                    if first_pos[from_id] < first_pos[to_id]:
                        G.add_node(
                            (from_id, from_orient), SR=gfa_message.get_SRank(from_id)
                        )
                        G.add_node((to_id, to_orient), SR=gfa_message.get_SRank(to_id))
                        G.add_edge(
                            (from_id, from_orient),
                            (to_id, to_orient),
                            SR=SRi,
                            weight=0,
                        )
                    else:
                        G.add_node(
                            (from_id, reverse_orient(from_orient)),
                            SR=gfa_message.get_SRank(from_id),
                        )
                        G.add_node(
                            (to_id, reverse_orient(to_orient)),
                            SR=gfa_message.get_SRank(to_id),
                        )
                        G.add_edge(
                            (to_id, reverse_orient(to_orient)),
                            (from_id, reverse_orient(from_orient)),
                            SR=SRi,
                            weight=0,
                        )

        else:
            if from_id in sum_double:
                if int(to_id[1:]) <= last_linear_int:
                    if to_orient == "+":
                        G.add_edge(
                            (from_id, from_orient),
                            (to_id, to_orient),
                            SR=SRi,
                            weight=0,
                        )
                        continue
                    elif to_orient == "-":
                        G.add_edge(
                            (to_id, reverse_orient(to_orient)),
                            (from_id, reverse_orient(from_orient)),
                            SR=SRi,
                            weight=0,
                        )
                        continue
                else:
                    if (to_id, to_orient) in G:
                        G.add_edge(
                            (from_id, from_orient),
                            (to_id, to_orient),
                            SR=SRi,
                            weight=0,
                        )
                    elif (to_id, reverse_orient(to_orient)) in G:
                        G.add_edge(
                            (to_id, reverse_orient(to_orient)),
                            (from_id, reverse_orient(from_orient)),
                            SR=SRi,
                            weight=0,
                        )
                    else:
                        flag = double_pass_one(from_id, to_id)
                        if flag == 1:
                            G.add_node(
                                (to_id, reverse_orient(to_orient)),
                                SR=gfa_message.get_SRank(to_id),
                            )
                            G.add_edge(
                                (to_id, reverse_orient(to_orient)),
                                (from_id, reverse_orient(from_orient)),
                                SR=SRi,
                                weight=0,
                            )
                        elif flag == -1:
                            G.add_node(
                                (to_id, to_orient),
                                SR=gfa_message.get_SRank(to_id),
                            )
                            G.add_edge(
                                (from_id, from_orient),
                                (to_id, to_orient),
                                SR=SRi,
                                weight=0,
                            )
                        elif flag == 0:
                            cache_link.append(
                                (from_id, from_orient, to_id, to_orient, SRi, 1)
                            )
            elif to_id in sum_double:
                if int(from_id[1:]) <= last_linear_int:
                    if from_orient == "+":
                        G.add_edge(
                            (from_id, from_orient),
                            (to_id, to_orient),
                            SR=SRi,
                            weight=0,
                        )
                        continue
                    elif from_orient == "-":
                        G.add_edge(
                            (to_id, reverse_orient(to_orient)),
                            (from_id, reverse_orient(from_orient)),
                            SR=SRi,
                            weight=0,
                        )
                        continue
                else:
                    if (from_id, from_orient) in G:
                        G.add_edge(
                            (from_id, from_orient),
                            (to_id, to_orient),
                            SR=SRi,
                            weight=0,
                        )
                    elif (from_id, reverse_orient(from_orient)) in G:
                        G.add_edge(
                            (to_id, reverse_orient(to_orient)),
                            (from_id, reverse_orient(from_orient)),
                            SR=SRi,
                            weight=0,
                        )
                    else:
                        flag = double_pass_one(to_id, from_id)
                        if flag == 1:
                            G.add_node(
                                (from_id, from_orient),
                                SR=gfa_message.get_SRank(from_id),
                            )
                            G.add_edge(
                                (from_id, from_orient),
                                (to_id, to_orient),
                                SR=SRi,
                                weight=0,
                            )
                        elif flag == -1:
                            G.add_node(
                                (from_id, reverse_orient(from_orient)),
                                SR=gfa_message.get_SRank(from_id),
                            )
                            G.add_edge(
                                (to_id, reverse_orient(to_orient)),
                                (from_id, reverse_orient(from_orient)),
                                SR=SRi,
                                weight=0,
                            )
                        elif flag == 0:
                            cache_link.append(
                                (from_id, from_orient, to_id, to_orient, SRi, 2)
                            )
    last_n = 0
    while cache_link:
        n = len(cache_link)  # The number of original elements in this round
        if last_n == n:
            break
        for _ in range(n):
            link = cache_link.popleft()
            from_id, from_orient, to_id, to_orient, SRi, flag = link
            if flag == 0:
                if (from_id, from_orient) in G or (to_id, to_orient) in G:
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                elif (from_id, reverse_orient(from_orient)) in G or (
                    to_id,
                    reverse_orient(to_orient),
                ) in G:
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
                else:
                    cache_link.append(link)
            if flag == 1:
                if (to_id, to_orient) in G:
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                elif (to_id, reverse_orient(to_orient)) in G:
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
                else:
                    cache_link.append(link)
            if flag == 2:
                if (from_id, from_orient) in G:
                    G.add_edge(
                        (from_id, from_orient),
                        (to_id, to_orient),
                        SR=SRi,
                        weight=0,
                    )
                elif (from_id, reverse_orient(from_orient)) in G:
                    G.add_edge(
                        (to_id, reverse_orient(to_orient)),
                        (from_id, reverse_orient(from_orient)),
                        SR=SRi,
                        weight=0,
                    )
                else:
                    cache_link.append(link)
        last_n = n
        logger.info(
            f"This round of processing is completed, a total of {n} edges, and the next round is about to begin"
        )
    if len(cache_link) > 0:
        logger.warning(f"cache_link remaining:{[(x[0],x[2]) for x in cache_link]}")
    sources_graph = [n for n in G.nodes() if G.in_degree(n) == 0]
    sinks_graph = [n for n in G.nodes() if G.out_degree(n) == 0]
    chrNum = len(sources.keys())
    while len(sources_graph) > chrNum:
        for x in sources_graph[chrNum:]:
            G.remove_node(x)
        sources_graph = [n for n in G.nodes() if G.in_degree(n) == 0]
        logger.debug(sources_graph)
    while len(sinks_graph) > chrNum:
        for x in sinks_graph[chrNum:]:
            G.remove_node(x)
        sinks_graph = [n for n in G.nodes() if G.out_degree(n) == 0]
        logger.debug(sinks_graph)
    return G


def _turn_GFA_to_DiGraph_simple(gfa_message: Minigfa) -> nx.DiGraph:

    def reverse_orient(one_orient):
        return "+" if one_orient == "-" else "-"

    G = nx.DiGraph()
    for from_id, from_orient, to_id, to_orient, SRi in gfa_message.get_all_Link():
        u = (from_id, from_orient)
        v = (to_id, to_orient)
        G.add_node(u, SR=gfa_message.get_SRank(from_id))
        G.add_node(v, SR=gfa_message.get_SRank(to_id))
        G.add_edge(u, v, SR=SRi, weight=0)
        ref_u = (from_id, reverse_orient(from_orient))
        ref_v = (to_id, reverse_orient(to_orient))
        G.add_node(ref_u, SR=gfa_message.get_SRank(from_id))
        G.add_node(ref_v, SR=gfa_message.get_SRank(to_id))
        G.add_edge(ref_v, ref_u, Sr=SRi, weight=0)
    return G


# todo
def turn_GFA_to_DiGraph(
    gfa_message: Minigfa,
    bed_message: Optional[Minibed] = None,
    is_output_inspection_results_in_graph: bool = False,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> nx.DiGraph:
    """
    Convert the GFA file information and Bed file information (optional, if default, you may find some loops or paths that should not exist in your graph) into a directed graph.
    You can choose whether to output the key parameter information of the graph and whether to save it.
    Notice: It will generate a tmp folder storage under your working folder

    Args:
        gfa_maessage (Minigfa): Composite data storing GFA file information
        bed_message (Minibed | None, optional): Composite data storing Bed file information. Defaults to None.
        is_output_inspection_results (bool, optional): Whether to output the key parameters of the graph to stdout. Defaults to False.
        is_saved_as_pickle (bool, optional): Whether to save as a pickle file for reuse. Defaults to False.
        file_path (str | None, optional): If you choose to save as a pickle file,the graph will be saved in `file_path`. By default, the file name will be `myMinigraph.pl` in folder /tmp under your working folder.

    Returns:
        DiGraph: Directed graph representing the pan-genome.We think direction is 5' end to 3' end as you follow the diagram,and the ID of each node is a tuple, the first element is the ID of the segment in GFA, and the second element is the symbol "+" or "-". "-" represents the reverse complementary sequence of the connected segment sequence
    """
    starttime = time.time()
    if bed_message is None:
        Minigraph = _turn_GFA_to_DiGraph_simple(gfa_message)
    elif bed_message is not None:
        Minigraph = _turn_GFA_to_DiGraph_complex(gfa_message, bed_message)
    logger.info(
        "Finish turn GFA to Digraph in %0.2f seconds." % (time.time() - starttime)
    )
    if is_saved_as_pickle:
        if file_path is None:
            _save_to_tmp(Minigraph, "myMinigraph.pl")
        else:
            _save_graph(Minigraph, file_path)
    if is_output_inspection_results_in_graph:
        count1 = nx.number_weakly_connected_components(Minigraph)
        print("=================Key information of Minigraph=================")
        print("The number of weakly connected components of the graph: ", count1)
        print("Number of nodes in the graph: ", nx.number_of_nodes(Minigraph))
        print("Number of edges in the graph: ", nx.number_of_edges(Minigraph))
        sources = [n for n in Minigraph.nodes() if Minigraph.in_degree(n) == 0]
        sinks = [n for n in Minigraph.nodes() if Minigraph.out_degree(n) == 0]
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
    return Minigraph


if __name__ == "__main__":
    pass
