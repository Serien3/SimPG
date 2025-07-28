from . import *
import argparse


def cli():
    parser = argparse.ArgumentParser(
        prog="SimPG",
        description="A population genome simulation tool based on large-scale pangenomic data",
    )
    parser.add_argument("GFA_input", help="")
    parser.add_argument("BED_input", help="")
    parser.add_argument("SampleID_file_for_simulated_population", help="")
    parser.add_argument("--save_temp_file", action="store_true", help="")
    parser.add_argument("--save_walk_filepath", default=None, help="")
    parser.add_argument("--report_graph_information", action="store_true", help="")
    parser.add_argument("-o", "--simulation_result_output", default=None, help="")
    parser.add_argument("--population_name", default=None, help="")
    parser.add_argument("-n", "--sim_num", type=int, default=1, help="")
    parser.add_argument("--logging_verbose", action="store_true", help="")
    parser.add_argument("--is_human", action="store_true", help="")

    args = parser.parse_args()
    run_SimPG(
        args.GFA_input,
        args.BED_input,
        args.SampleID_file_for_simulated_population,
        args.save_temp_file,
        args.save_walk_filepath,
        args.report_graph_information,
        args.simulation_result_output,
        args.is_human,
        args.population_name,
        args.sim_num,
        args.logging_verbose,
    )


if __name__ == "__main__":
    cli()
