from SimPG import run_SimPG
import argparse


def cli():
    parser = argparse.ArgumentParser(
        prog="SimPG",
        description="A population genome simulation tool based on large-scale pangenomic data",
    )
    parser.add_argument("GFA_input", help="GFA file path")
    parser.add_argument("BED_input", help="BED file path")
    parser.add_argument(
        "SampleID_file_for_simulated_population",
        help="SampleID file for simulated population",
    )
    parser.add_argument(
        "--save_temp_file",
        action="store_true",
        help="Whether to save as a pickle file for reuse. Defaults to False",
    )
    parser.add_argument(
        "--save_walk_filepath",
        default=None,
        help="Save file location.By default, it is saved in `my_walks.pl` in the `/tmp` folder of the working directory.",
    )
    parser.add_argument(
        "--report_graph_information",
        action="store_true",
        help="Whether to output the key parameters of the graph to stdout. Defaults to False.",
    )
    parser.add_argument(
        "-o",
        "--simulation_result_output",
        default=None,
        help="Result output location. By default, the output is in the working directory.",
    )
    parser.add_argument(
        "--population_name",
        default=None,
        help="Give your simulated crowd a name, which will also be used as the prefix for the output files. Defaults to `My`.",
    )
    parser.add_argument(
        "-n",
        "--sim_num",
        type=int,
        default=1,
        help=" Number of simulations. Defaults to 1.",
    )
    parser.add_argument(
        "--logging_verbose",
        action="store_true",
        help="Whether to set the log output information level to at least `INFO` level .Default to `False`, set to `Warning` level.",
    )
    parser.add_argument(
        "--is_human",
        action="store_true",
        help="Is the pan-genome a human pan-genome? Defaults to `False`.",
    )

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
