from classes import *
from core import *
from utils import *


def run_SimPG(
    GFA_file_path: str,
    BED_file_path: str,
    population: list[str] | str,
    enable_to_save_temporary_folder: bool = False,
    every_sample_Whole_Genome_Sequencing_filepath: Optional[str] = None,
    whether_to_output_graph_information_in_terminal: bool = False,
    sim_file_out_folder: Optional[str] = None,
    is_human: bool = False,
    population_name: Optional[str] = None,
    sim_num: int = 1,
    logging_verbose: bool = False,
) -> None:
    set_default_logging(logging_verbose)
    gfa_message = Minigfa(GFA_file_path)
    bed_message = Minibed(BED_file_path)
    Minigraph = turn_GFA_to_DiGraph(
        gfa_message,
        bed_message,
        whether_to_output_graph_information_in_terminal,
        enable_to_save_temporary_folder,
    )
    simulate_population_every_walk(
        gfa_message,
        bed_message,
        Minigraph,
        population,
        every_sample_Whole_Genome_Sequencing_filepath,
    )
    core_seg_set = get_coreSeg_in_Pangenome(
        gfa_message,
        every_sample_Whole_Genome_Sequencing_filepath,
        enable_to_save_temporary_folder,
    )
    Pangenome_Digraph = simulate_Population_Pangenome(
        bed_message,
        every_sample_Whole_Genome_Sequencing_filepath,
        is_output_inspection_results_in_graph=whether_to_output_graph_information_in_terminal,
        is_saved_as_pickle=enable_to_save_temporary_folder,
    )
    simulate_Whole_Genome_Sequencing_for_population(
        Pangenome_Digraph,
        gfa_message,
        core_seg_set,
        sim_file_out_folder,
        is_human,
        population_name,
        sim_num,
    )


if __name__ == "__main__":
    pass
