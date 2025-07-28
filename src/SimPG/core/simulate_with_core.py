import pickle
import networkx as nx
import random
from . import logger
import json
from ..classes import Minigfa
from typing import Optional
import re
import os
import time

__all__ = ["simulate_Whole_Genome_Sequencing_for_population"]


def _generate_random_edge_weights(G: nx.DiGraph) -> dict:
    """
    Generate a random weight ~ Uniform(0,1) for each edge (u,v) in the graph G,
    Returns a dictionary of the form {(u,v): float}.
    """
    edge_weight = {}
    for u, v in G.edges():
        edge_weight[(u, v)] = random.random()
    return edge_weight


def _save_edge_weights_pickle(edge_weight: dict, filepath: str) -> None:
    with open(filepath, "wb") as f:
        pickle.dump(edge_weight, f)


def _load_edge_weights_pickle(filepath: str) -> dict:
    with open(filepath, "rb") as f:
        return pickle.load(f)


def _edge_dict_to_json_serializable(edge_weight: dict) -> dict:
    json_dict = {}
    for (u, v), w in edge_weight.items():
        key_str = f"{u}|{v}"
        json_dict[key_str] = w
    return json_dict


def _json_serializable_to_edge_dict(json_dict: dict) -> dict:
    edge_weight = {}
    for key_str, w in json_dict.items():
        u_str, v_str = key_str.split("|", 1)
        try:
            u = int(u_str)
        except ValueError:
            u = u_str
        try:
            v = int(v_str)
        except ValueError:
            v = v_str
        edge_weight[(u, v)] = w
    return edge_weight


def _save_edge_weights_json(edge_weight: dict, filepath: str) -> None:
    serializable = _edge_dict_to_json_serializable(edge_weight)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=2, ensure_ascii=False)


def _load_edge_weights_json(filepath: str) -> dict:
    with open(filepath, "r", encoding="utf-8") as f:
        json_dict = json.load(f)
    return _json_serializable_to_edge_dict(json_dict)


def _build_successors_weighted(G: nx.DiGraph, edge_weight: dict) -> dict:
    """
    Construct the successor list [(v, w(u,v)), ...] for each node given the (u,v)->w dictionary.
    Return a dict：{ u: [(v1, w1), (v2, w2), ...], ... }.
    """
    successors_weighted = {}
    for u in G.nodes():
        lst = []
        for v in G.successors(u):
            lst.append((v, edge_weight[(u, v)]))
        successors_weighted[u] = lst
    return successors_weighted


def _random_walk_with_loaded_weights(
    G: nx.DiGraph, s: tuple, t: tuple, edge_weight, max_steps: Optional[int] = None
) -> list:
    """
    Use the loaded edge_weight (i.e. (u,v)->float) to walk, each step:
      - Filter out the visited nodes and only do weighted random selection in the candidate list;
      - If t is reached, return path immediately;
      - If dead-end or more than max_steps, restart a new round of walk from s.

    Args:
        G: DiGraph
        s: start
        t: end
        edge_weight: Generated dictionary {(u,v): float}
        max_steps: The maximum number of steps allowed in a single walk (optional)
    返回:
        path: An acyclic path from s to t (List), if s==t directly return [s].
    """
    # Parameter Check
    if s not in G:
        raise ValueError(f"The starting point {s!r} is not in the graph G")
    if t not in G:
        raise ValueError(f"The ending point {t!r} is not in the graph G")
    if s == t:
        return [s]
    if not nx.has_path(G, s, t):
        raise nx.NetworkXNoPath(
            f"There does not exist any path from {s!r} to {t!r} in the graph."
        )

    # Pre-construct successors_weighted to avoid repeated construction in each round of walks
    successors_weighted = _build_successors_weighted(G, edge_weight)

    # Precompute “which nodes can reach t”
    #    Using the reverse graph, do a DFS/BFS from t
    G_rev = G.reverse(copy=False)
    reachable_to_t = set(nx.descendants(G_rev, t))
    reachable_to_t.add(t)

    while True:
        cur = s
        path = [s]
        visited = {s}
        steps = 0

        while True:
            if cur == t:
                return path
            if max_steps is not None and steps >= max_steps:
                break

            # Filter out the successors that have already been visited
            # cand = [(v, w) for (v, w) in successors_weighted[cur] if v not in visited]
            cand = [
                (v, w)
                for (v, w) in successors_weighted[cur]
                if v not in visited and v in reachable_to_t
            ]
            if not cand:
                break

            # Normalized random selection by weight
            total_w = sum(w for (_, w) in cand)
            r = random.random() * total_w
            cum = 0.0
            chosen = None
            for v, w in cand:
                cum += w
                if r <= cum:
                    chosen = v
                    break
            if chosen is not None:
                path.append(chosen)
                visited.add(chosen)
            cur = chosen
            steps += 1

        # This time it was unsuccessful, so we restart from the outer while.


def _random_walk_find_path_acyclic(
    G: nx.DiGraph, s, t, max_steps: Optional[int] = None
) -> list:
    """
    In a directed graph G, start a random walk from the starting point s and record the visited nodes to avoid loop revisiting.
    Once the end point t is reached, the path is returned; if the current node does not have any unvisited outgoing edges, the current walk is considered to have failed, and a new walk starts from s.

    Args:
        G: DiGraph
        s: start
        t: end
        max_steps: he maximum number of steps allowed in a single walk (optional)

    返回:
        A list `path`, which represents the complete sequence of nodes when the random walk first reaches t, for example, [s, ..., t].
        If s == t,return [s]。
    """
    if s not in G:
        raise ValueError(f"The starting point {s!r} is not in the graph G")
    if t not in G:
        raise ValueError(f"The end point {t!r} is not in the graph G")
    if s == t:
        return [s]

    # You can quickly determine: if s cannot reach t, there is no need to start a random walk
    if not nx.has_path(G, s, t):
        raise nx.NetworkXNoPath(
            f"There does not exist any path from {s!r} to {t!r} in the graph."
        )

    while True:
        # Each round of walking starts from the starting point
        cur = s
        path = [s]
        # Record the nodes that have been visited in this walk
        visited = {s}
        step_count = 0

        while True:
            # If you reach the end point, return to the path
            if cur == t:
                return path

            # 如If the maximum number of steps is set and it exceeds the limit, the process will be considered a failure.
            if max_steps is not None and step_count >= max_steps:
                break  # This roaming failed, jump to the outer layer and start again

            # Find all "unvisited" successor nodes
            all_neighbors = list(G.successors(cur))
            # Filter out already visited nodes from successor nodes
            unvisited_neighbors = [v for v in all_neighbors if v not in visited]

            # If there are no "unvisited" outgoing edges, it means that a dead end has been reached and the walk has failed.
            if not unvisited_neighbors:
                break

            # Uniformly randomly select one of all unvisited successor nodes
            next_node = random.choice(unvisited_neighbors)

            # Update the current node, path, visited and number of steps
            path.append(next_node)
            visited.add(next_node)
            cur = next_node
            step_count += 1

        # If you reach here, it means that this walk has not reached t, return to the outer loop and try again


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


def _reverse_complement(seq: str) -> str:
    """
    Convert a DNA sequence seq to its reverse complement.

    steps:
      1. Complementarity:A↔T, C↔G (Also compatible with lowercase).
      2. Reverse the entire sequence.

    Args:
        seq: Original DNA sequence,only include A, T, C, G(Also compatible with lowercase).

    Return:
        Reverse complementary sequence (also retains the original uppercase and lowercase letter pattern).
    """
    # 1. Constructing a mapping table
    trans_table = str.maketrans("ATCGatcg", "TAGCtagc")
    # 2. Translation complementation + inversion
    return seq.translate(trans_table)[::-1]


def _get_chr_graph(
    Pangenome_graph: nx.DiGraph, weak_comps: list[list[tuple[str, str]]]
):
    for idx, c in enumerate(weak_comps, start=1):
        yield idx, Pangenome_graph.subgraph(c).copy()


def _write_vcf(fileVCF, parts: list, gfa_messsage: Minigfa):

    def split_by_predicate(lst, predicate):
        """
        Divide `lst` into segments according to `predicate` and return the `(start_idx, end_idx)` of each segment in the original list.
        - Every time predicate(item) == True is encountered, the current segment [start, i] is recorded,
          Then set start = i .
        - Finally, add the section at the end that does not contain special elements.
        """
        segments = []
        start = 0
        for i, item in enumerate(lst):
            if predicate(item):
                # Record the sublist from start to i (including special elements)
                segments.append(lst[start : i + 1])
                # Start the next paragraph from the current element.
                start = i
        # If there is still some left in the last segment (start < len-1), also record it
        if start < len(lst) - 1:
            segments.append(lst[start : len(lst)])
        return segments

    Mutation_include = split_by_predicate(
        parts,
        lambda seq: gfa_messsage.get_source_sample(seq[0])
        == gfa_messsage.get_linear_reference()
        and seq[1] == "+",
    )
    for mutation in Mutation_include:
        # Regular expression processing
        start = re.search(r"\d+", mutation[0][0])
        end = re.search(r"\d+", mutation[-1][0])
        if start and end:
            num_start = int(start.group())
            num_end = int(end.group())
        else:
            raise KeyError("re failed")
        if (len(mutation) == 1) or (len(mutation) == 2 and num_start + 1 == num_end):
            continue
        if num_start > num_end:
            logger.warning(
                f"Find start -> end , have circle in s{num_start} to s{num_end},{mutation},jump over"
            )
            continue
        fileVCF.write(f"({mutation[0][0]},{mutation[-1][0]})\t")
        lst_linear = [f"s{i}" for i in range(num_start, num_end + 1)]
        fileVCF.write("(" + lst_linear[0] + "+")
        for x in lst_linear[1:]:
            fileVCF.write("," + x + "+")
        fileVCF.write(")\t")
        fileVCF.write(f"(" + mutation[0][0] + mutation[0][1])
        for x, orient in mutation[1:]:
            fileVCF.write("," + x + orient)
        fileVCF.write(")\n")


def _ensure_dir_for_file(file_path):
    """
    Ensure that the directory above file_path exists, and create it if it does not exist (including multiple levels of directories).
    """
    dir_path = os.path.dirname(file_path)
    if dir_path and not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)


# todo
def simulate_Whole_Genome_Sequencing_for_population(
    Pangenome_graph: nx.DiGraph,
    gfa_message: Minigfa,
    coreSeg: set[tuple[str, str]],
    file_out_folder: Optional[str] = None,
    is_human: bool = False,
    population_name: Optional[str] = None,
    sim_num: int = 1,
) -> None:
    """
    Without using the method of pre-setting and saving the weight matrix, directly randomly select the next node to walk.
    Notice: This function will generate two folders under `file_out`, namely `{population_name}_simulate_fasta` and `{population_name}_simulate_vcf`. Thefasta files and vcf files will be saved in the following folders respectively.

    Args:
        Pangenome_graph (nx.DiGraph): Pan-genome graph
        gfa_message (Minigfa): Composite data storing GFA file information.
        coreSeg (set[tuple[str, str]]): A collection of core sequence nodes
        file_out_folder (str |None, optional): Result output location. By default, the output is in the working directory.
        is_human (bool, optional): Is the pan-genome a human pan-genome? Defaults to False.
        population_name (str |None, optional): Give your simulated crowd a name, which will also be used as the prefix for the output files. Defaults to "My".
        sim_num (int, optional): Number of simulations. Defaults to 1.

    Raises:
        ValueError: The start or end node is not in the graph.
        networkx.NetworkXNoPath: There is no path from the start node to the end node

    """
    if file_out_folder is None:
        file_out_folder = os.getcwd()
    if population_name is None:
        population_name = "my"
    for Number in range(1, sim_num + 1):
        logger.info(f"start simulate {Number:03d}")
        weak_comps = list(nx.weakly_connected_components(Pangenome_graph))
        _ensure_dir_for_file(
            file_path=f"{file_out_folder}/{population_name}_simulate_fasta/{population_name}_simulate{Number:03d}.fa"
        )
        _ensure_dir_for_file(
            file_path=f"{file_out_folder}/{population_name}_simulate_rvcf/{population_name}_simulate{Number:03d}.rvcf"
        )
        starttime = time.time()
        with open(
            f"{file_out_folder}/{population_name}_simulate_fasta/{population_name}_simulate{Number:03d}.fa",
            "w",
        ) as fileFa, open(
            f"{file_out_folder}/{population_name}_simulate_rvcf/{population_name}_simulate{Number:03d}.rvcf",
            "w",
        ) as fileVCF:
            for idx, chr_graph in _get_chr_graph(Pangenome_graph, weak_comps):
                if is_human:
                    if idx == 23:
                        fileFa.write(f">chrX\n")
                    elif idx == 24:
                        fileFa.write(">chrY\n")
                    else:
                        fileFa.write(f">chr{idx}\n")
                else:
                    fileFa.write(f">chr{idx}\n")
                coreSeg_inchr = coreSeg & {n for n in chr_graph.nodes()}
                sorted_coreSeg_inchr = sorted(
                    coreSeg_inchr, key=lambda item: int(item[0][1:])
                )
                fileFa.write(gfa_message.get_seq(sorted_coreSeg_inchr[0][0]))
                for u, v in zip(sorted_coreSeg_inchr[:-1], sorted_coreSeg_inchr[1:]):
                    parts = _random_walk_find_path_acyclic(chr_graph, u, v, 1000)
                    # if len(parts) == 0:
                    #     raise KeyError(f"Error, from {u} to {v} not find path")
                    for segID, orient in parts[1:]:
                        seq_out = gfa_message.get_seq(segID)
                        if orient == "-":
                            seq_out = _reverse_complement(seq_out)
                        fileFa.write(seq_out)
                    seq_linear = _generate_sequence(u[0], v[0])
                    if seq_linear != parts:
                        _write_vcf(fileVCF, parts, gfa_message)
                fileFa.write("\n")
                logger.debug(f"Finish simulate chr{idx}")
        logger.info("Finish a simulation in %0.2f seconds." % (time.time() - starttime))


if __name__ == "__main__":
    pass
