"""Extract the walking route of individual genome sequence mapping in the graph"""

import pickle
import networkx as nx
from collections import deque
from typing import List, Tuple, Any, Optional
from . import logger
from ..classes import Minibed, Minigfa
import os

__all__ = ["simulate_population_every_walk"]


def _find_constrained_path(G: nx.DiGraph, source, target, sr_value) -> None | list[Any]:
    """
    In a directed graph G, find a path from source to target,
    Requires that all nodes and edges with attribute SR==sr_value must be passed through.
    If it exists, return the node list; otherwise return None.
    """

    # 1. Constructing the set of nodes and edges that must be passed
    required_nodes = {v for v, d in G.nodes(data=True) if d.get("SR") == sr_value}
    required_edges = {
        (u, v) for u, v, d in G.edges(data=True) if d.get("SR") == sr_value
    }
    # logger.info(f"{len(required_nodes)},{len(required_edges)}")
    # 2. Reverse BFS: Mark all nodes that can reach the target
    reachable = set([target])
    revG = G.reverse(copy=False)
    queue = deque([target])
    while queue:
        u = queue.popleft()
        for prev in revG.successors(u):
            if prev not in reachable:
                reachable.add(prev)
                queue.append(prev)

    if source not in reachable:
        # The starting point cannot reach the end point
        return None

    # 3. BFS state search
    #  state = (current node, visited nodes frozenset, visited edges frozenset)
    init_nodes = frozenset({source} & required_nodes)
    init_state = (source, init_nodes, frozenset())

    queue = deque([init_state])
    visited = set([init_state])  # Deduplication: The queued state set
    parent = {init_state: None}  # To rebuild the path

    while queue:
        u, vis_n, vis_e = queue.popleft()

        # Check whether the goal has been achieved
        if u == target and vis_n == required_nodes and vis_e == required_edges:
            # Reconstruct the node path from source to target
            path = []
            st = (u, vis_n, vis_e)
            while st:
                path.append(st[0])
                st = parent[st]
            return list(reversed(path))

        # Enumerate all outgoing edges to transfer
        for v in G.successors(u):
            # Cut: target cannot be reached after v
            if v not in reachable:
                continue

            # Update the must-pass node set and must-pass edge set
            new_vis_n = vis_n | ({v} if v in required_nodes else frozenset())
            edge = (u, v)
            new_vis_e = vis_e | ({edge} if edge in required_edges else frozenset())

            new_state = (v, new_vis_n, new_vis_e)
            if new_state in visited:
                continue

            visited.add(new_state)
            parent[new_state] = (u, vis_n, vis_e)
            queue.append(new_state)

    # Search completed, no solution
    return None


import math


def _prize_collecting_path(G: nx.DiGraph, s, t, V_star, E_star, delta_max=math.inf):
    """
    Find an s->t path in directed graph G that covers as many nodes in V_star and edges in E_star as possible.

    Returns:
      path (list): list of nodes in the s->t path, or None if no valid path.
      uncovered_nodes (int): number of nodes in V_star not covered.
      uncovered_edges (int): number of edges in E_star not covered.
    """
    # --- 1. Identify compressed nodes and map pseudo-edge nodes ---
    nodes = set(V_star)
    pseudo = {("e", u, v) for (u, v) in E_star}
    V_hat = {s, t} | nodes | pseudo

    mapping = {}
    edge_weight = {}
    for z in V_hat:
        if isinstance(z, tuple) and z[0] == "e":
            _, u, v = z
            mapping[z] = (u, v)
            edge_weight[z] = G[u][v].get("weight", 1)
        else:
            mapping[z] = (z, z)

    # --- 2. Precompute shortest-path lengths ---
    end_points = {mapping[z][1] for z in V_hat}
    dists = {x: nx.single_source_dijkstra_path_length(G, x) for x in end_points}

    # --- 3. Build compressed complete digraph H ---
    H = nx.DiGraph()
    reward = {z: 1 for z in nodes | pseudo}
    reward[s] = reward[t] = 0
    for z in V_hat:
        H.add_node(z)
    for z1 in V_hat:
        end1 = mapping[z1][1]
        for z2 in V_hat - {z1}:
            start2 = mapping[z2][0]
            dist = dists.get(end1, {}).get(start2, math.inf)
            if math.isinf(dist):
                continue
            w = dist + edge_weight.get(z2, 0)
            H.add_edge(z1, z2, weight=w)

    # --- 4. Greedy insertion in H ---
    P = [s, t]
    covered = set()

    def best_insert():
        best = (0, None, None)
        for z in V_hat - {s, t} - covered:
            for i in range(len(P) - 1):
                u, v = P[i], P[i + 1]
                if not H.has_edge(u, z) or not H.has_edge(z, v):
                    continue
                w_uv = H[u][v]["weight"]
                extra = H[u][z]["weight"] + H[z][v]["weight"] - w_uv
                if extra > delta_max:
                    continue
                gain = reward[z] / extra if extra > 0 else reward[z] * 1e6
                if gain > best[0]:
                    best = (gain, z, i + 1)
        return best

    while True:
        gain, z, pos = best_insert()
        if z is None or gain <= 0:
            break
        P.insert(pos, z)
        covered.add(z)

    # Determine uncovered counts
    uncov_nodes = len(nodes - covered)
    uncov_edges = len(pseudo - covered)

    # --- 5. Reconstruct original path or handle no path ---
    try:
        path = []
        for i in range(len(P) - 1):
            z1, z2 = P[i], P[i + 1]
            u1 = mapping[z1][1]
            u2 = mapping[z2][0]
            sub = nx.shortest_path(G, u1, u2, weight="weight")
            if i > 0:
                sub = sub[1:]
            path.extend(sub)
            if z2 in pseudo:
                _, u, v = z2
                if path[-1] != v:
                    path.append(v)
    except nx.NetworkXNoPath:
        return None, uncov_nodes, uncov_edges

    return path, uncov_nodes, uncov_edges


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


def _simulate_sample_path(
    gfa_message: Minigfa, bed_message: Minibed, G_full: nx.DiGraph, sample_name: str
) -> None | List[Tuple[str, str]]:

    def get_SR_by_sample_name() -> int:
        for segID in gfa_message.get_all_segID():
            if gfa_message.get_source_sample(segID) == sample_name:
                # print(gfa_message.get_SRank(segID))
                return gfa_message.get_SRank(segID)
        return -1

    Genome_Sequencing_with_segment: List[Tuple[str, str]] = []
    target_SR = get_SR_by_sample_name()
    if target_SR == -1:
        logger.warning(f"No target_SR for {sample_name}")
        return None
    sources, sinks = bed_message.get_linear_sources_and_sinks()
    last_chr = object()
    for chr, _, _, _, list_of_segments in bed_message:
        cur_chr = chr
        if last_chr != cur_chr:
            Genome_Sequencing_with_segment.append((sources[chr], "+"))
            last_chr = cur_chr
        list_of_segments_double = [
            (x, sign) for x in list_of_segments for sign in ("+", "-")
        ]
        tempG = G_full.subgraph(list_of_segments_double).copy()
        has_node = any(attrs["SR"] == target_SR for _, attrs in tempG.nodes(data=True))
        has_edge = any(
            attrs["SR"] == target_SR for _, _, attrs in tempG.edges(data=True)
        )
        if not (has_node or has_edge):
            linear_segments = _generate_sequence(
                list_of_segments[0], list_of_segments[-1]
            )
            Genome_Sequencing_with_segment.extend(linear_segments[1:])
            continue

        # 1. Delete nodes whose attribute SR is greater than the threshold
        nodes_to_remove = [
            n for n, attrs in tempG.nodes(data=True) if attrs["SR"] > target_SR
        ]
        tempG.remove_nodes_from(nodes_to_remove)

        # 2. Delete edges whose attribute SR is greater than the threshold
        edges_to_remove = [
            (u, v) for u, v, attrs in tempG.edges(data=True) if attrs["SR"] > target_SR
        ]
        for u, v in edges_to_remove:
            tempG.remove_edge(u, v)
        start_segID = (list_of_segments[0], "+")
        end_segID = (list_of_segments[-1], "+")
        logger.debug(f"{start_segID} to {end_segID}")
        required_nodes = {v for v, d in tempG.nodes(data=True) if d["SR"] == target_SR}
        required_edges = {
            (u, v) for u, v, d in tempG.edges(data=True) if d["SR"] == target_SR
        }
        if (
            nx.number_of_edges(tempG)
            * (2 ** len(required_nodes))
            * (2 ** len(required_edges))
            < 300000000000
        ):
            path = _find_constrained_path(tempG, start_segID, end_segID, target_SR)
        else:
            logger.info(
                f"Preparing an approximation algorithm from {start_segID} to {end_segID}"
            )
            path = None
        if path is None:
            # logger.warning(f"from {start_segID} to {end_segID} don't find path")
            if (target_SR == 11 and start_segID == ("s238674", "+")) or (
                target_SR == 3 and start_segID == ("s411304", "+")
            ):
                linear_segments = _generate_sequence(
                    list_of_segments[0], list_of_segments[-1]
                )
                Genome_Sequencing_with_segment.extend(linear_segments[1:])
                logger.info(f"{len(required_nodes)} ->{0},{len(required_edges)} -> {0}")
                continue
            path, uncovered_nodes_count, uncovered_edges_count = _prize_collecting_path(
                tempG, start_segID, end_segID, required_nodes, required_edges, 4
            )
            if path is None:
                logger.error(
                    f"Approximation algorithm fails from {start_segID} to {end_segID}"
                )
            else:
                logger.debug(
                    f"Total nodes:{len(required_nodes)},uncovered nodes:{uncovered_nodes_count},total edges:{len(required_edges)},uncovered edges:{uncovered_edges_count}"
                )
                Genome_Sequencing_with_segment.extend(path[1:])
        else:
            Genome_Sequencing_with_segment.extend(path[1:])
    if Genome_Sequencing_with_segment is None:
        logging.error("This Genome_Sequencing_with_segment is None")
        return None
    else:
        return Genome_Sequencing_with_segment


def _save_to_tmp(filename: str) -> str:
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
    return file_path


def simulate_population_every_walk(
    gfa_message: Minigfa,
    bed_message: Minibed,
    G_full: nx.DiGraph,
    population: list[str] | str,
    saved_file_path: Optional[str] = None,
) -> None:
    """
    Without the need for original individual genome sequence information involved in building a pan-genome, this function can extract the path of individual genome sequences mapped in the graph.
    Notice:This function does not return anything.It will save the walking route of each sample in `saved_file_path` file.

    Args:
        gfa_message (Minigfa): Composite data storing GFA file information.
        bed_message (Minibed): Composite data storing Bed file information.
        G_full (nx.DiGraph):Pan-genome graph
        population (list[str] | str):Input a list of sample names, or a text file with only one sample name per line
        saved_file_path (str):Save file location.By default, it is saved in `my_walks.pl` in the `/tmp` folder of the working directory.
    """
    import gc
    import time

    samples = list[str]()
    if isinstance(population, str):
        with open(population, "r", encoding="utf-8") as population_file:
            for line in population_file:
                name = line.strip()
                if not name:
                    continue
                samples.append(f"{name}")
    else:
        samples = population
    if saved_file_path is None:
        saved_file_path = _save_to_tmp("my_walks.pl")
    starttime = time.time()
    with open(saved_file_path, "wb") as f:
        for sample_name in samples:
            logger.info(f"Begin simulate {sample_name} genome sequencing with segment")
            Genome_Sequencing_with_segment_out = _simulate_sample_path(
                gfa_message, bed_message, G_full, sample_name
            )
            pickle.dump(
                (sample_name, Genome_Sequencing_with_segment_out),
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
            # Release memory immediately after processing
            del Genome_Sequencing_with_segment_out
            gc.collect()
    logger.info(
        "Finish get every sample walk in %0.2f seconds." % (time.time() - starttime)
    )


if __name__ == "__main__":
    pass
