import random
import os
from ..classes import Minibed, Minigfa
from . import logger
import time

__all__ = ["sim_part", "sim_part_for_num"]


def _keep_part_vcf(input_file: str, output_file: str, fraction: float):
    """
    Randomly and evenly retain lines from a file and output to a new file

    Args:
    input_file (str): Input file path
    output_file (str):Output file path
    fraction (float): Preserve row proportions (0.0 ~ 1.0)
    """
    # Verify the validity of the scale parameters
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("The scale must be between 0.0 and 1.0")

    # Using the with statement to safely open files
    with open(input_file, "r") as infile:
        # Read all lines
        lines = infile.readlines()

        # Special case handling if the ratio is 0 or 1
        if fraction == 0.0:
            selected_lines = []
        elif fraction == 1.0:
            selected_lines = lines
        else:
            # Calculate the number of rows to keep
            n_lines = len(lines)
            n_keep = round(n_lines * fraction)

            # Randomly sample row indices (no resampling)
            selected_indices = random.sample(range(n_lines), n_keep)
            # Select rows by index and keep original order
            selected_lines = [lines[i] for i in sorted(selected_indices)]

    # Writing output files
    with open(output_file, "w") as outfile:
        outfile.writelines(selected_lines)


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


def _get_seq_by_vcf(vcf_filrPath, bed_message: Minibed):
    # sources: list[tuple[str, str]] = []
    # sinks: list[tuple[str, str]] = []
    # for bedLine_list in bed_message.bed_line.values():
    #     sources.append((bedLine_list[0].list_of_segments[0], "+"))
    #     sinks.append((bedLine_list[-1].list_of_segments[-1], "+"))
    sources: dict[str, str] = {}
    sinks: dict[str, str] = {}
    sources, sinks = bed_message.get_linear_sources_and_sinks()
    chr_seq = []
    chr_num = 0
    for chr in sources.keys():
        chr_seq.append(_generate_sequence(sources[chr], sinks[chr]))
    with open(vcf_filrPath, "r") as file:
        lines = file.readlines()
        for line in lines:
            easy_line = line.strip().split("\t")
            s = easy_line[2]
            # Remove the left and right brackets
            inner = s.strip("()")
            # Split by comma to get a list
            items = inner.split(",")
            # Convert each string in the list to -> ("s1234","+")
            items = [(x[:-1], x[-1]) for x in items]
            if int(items[0][0][1:]) > int(chr_seq[chr_num][-1][0][1:]):
                chr_num += 1
            # 1. Find the index of the first element of items in the reference genome
            if items[0] not in chr_seq[chr_num]:
                print(f"Jump {items[0]}")
                continue
            start = chr_seq[chr_num].index(items[0])
            # 2. Find the index of the last element of items in the reference genome (starting from start)
            end = chr_seq[chr_num].index(items[-1], start)

            # 3. Replace this interval with a slice
            chr_seq[chr_num][start : end + 1] = items
    return chr_seq


def _writefa_with_vcf(
    gfa_message: Minigfa,
    chr_seq,
    filepath_fasta,
    is_human: bool = False,
):
    with open(filepath_fasta, "w") as fileFa:
        for idx, seq_list in enumerate(chr_seq, start=1):
            if is_human:
                if idx == 23:
                    fileFa.write(f">chrX\n")
                elif idx == 24:
                    fileFa.write(">chrY\n")
                else:
                    fileFa.write(f">chr{idx}\n")
            else:
                fileFa.write(f">chr{idx}\n")
            for seq in seq_list:
                seq_out = gfa_message.get_seq(seq[0])
                if seq[1] == "-":
                    seq_out = _reverse_complement(seq_out)
                # fileFa.write(seq[0] + seq[1])
                # fileFa.write("\n")
                fileFa.write(seq_out)
            fileFa.write("\n")


def _ensure_dir_for_file(file_path):
    """
    Ensure that the directory above file_path exists, and create it if it does not exist (including multiple directories).
    """
    dir_path = os.path.dirname(file_path)
    if dir_path and not os.path.exists(dir_path):
        logger.info(f"Not found {file_path},already makedir for it")
        os.makedirs(dir_path, exist_ok=True)


def sim_part(
    in_rvcf: str,
    out_fasta: str,
    out_rvcf: str,
    bed_message: Minibed,
    gfa_message: Minigfa,
    fraction: float,
    is_human: bool = False,
) -> None:
    """Keep the variation of fraction ratio in in_vcf file, and output the corresponding new fasta file & rvcf (unconverted vcf file) file

    Args:
        in_rvcf (str): Input txt(unconverted vcf file) file location
        out_fasta (str): Output fasta file location
        out_rvcf (str): Output rvcf(unconverted vcf file) file location
        bed_message (Minibed): Composite data storing Bed file information.
        gfa_message (Minigfa): Composite data storing GFA file information.
        fraction (float): The proportion of the number of retained variants
        is_human (bool, optional): Is the pan-genome a human pan-genome? Defaults to False.

    Raises:
        ValueError: The scale must be between 0.0 and 1.0. But the input parameters do not meet the requirements.
    """
    _keep_part_vcf(in_rvcf, out_rvcf, fraction)
    chr_seq_list = _get_seq_by_vcf(out_rvcf, bed_message)
    _writefa_with_vcf(gfa_message, chr_seq_list, out_fasta, is_human)


def sim_part_for_num(
    in_rvcf_folder: str,
    out_folder: str,
    bed_message: Minibed,
    gfa_message: Minigfa,
    population_name: str,
    fraction: float,
    is_human: bool = False,
    num=10,
) -> None:
    """Perform sim_part processing on the num vcf files in the in_rvcf_folder folder (the folder contains the rvcf of a population)

    Args:
        in_rvcf_folder (str): Input rvcf(unconverted vcf file) file folder location.
        out_folder (str): Output file location(will generate two folder under this folder:fasta folder and rvcf folder)
        bed_message (Minibed): Composite data storing Bed file information.
        gfa_message (Minigfa): Composite data storing GFA file information.
        population_name (str): in_vcf_folder's Population name
        fraction (float): The proportion of the number of retained variants
        is_human (bool, optional): Is the pan-genome a human pan-genome?. Defaults to False.
        num (int, optional): Number of rvcf files.The default is ten times.

    Raises:
        ValueError: The scale must be between 0.0 and 1.0. But the input parameters do not meet the requirements
    """
    for i in range(1, num + 1):
        starttime = time.time()
        _ensure_dir_for_file(
            file_path=out_folder
            + f"/{population_name}_simulate_fasta"
            + f"/{population_name}_simulate{i:03d}.fa"
        )
        _ensure_dir_for_file(
            file_path=out_folder
            + f"/{population_name}_simulate_rvcf"
            + f"/{population_name}_simulate{i:03d}.rvcf"
        )
        sim_part(
            in_rvcf_folder + f"/{population_name}_simulate{i:03d}.rvcf",
            out_folder
            + f"/{population_name}_simulate_fasta"
            + f"/{population_name}_simulate{i:03d}.fa",
            out_folder
            + f"/{population_name}_simulate_rvcf"
            + f"/{population_name}_simulate{i:03d}.rvcf",
            bed_message,
            gfa_message,
            fraction,
            is_human,
        )
        logger.info(
            "Finish a simulate_part in %0.2f seconds." % (time.time() - starttime)
        )


if __name__ == "__main__":
    pass
