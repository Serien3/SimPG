"""A script to help convert sampleID to standard"""

"""
    Use it with command:
    `python3 <filePath_in> <filePath_out> <sample_chromosome_ploidy>`
"""
import argparse


def load_samples_with_suffix(file_path: str, ploidy: int) -> list[str]:
    """
    Read sample names (one per line) from the txt file specified by file_path,
    Returns a list containing two strings with each sample name suffixed with .k
    """
    samples = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            name = line.strip()
            if not name:
                continue
            for i in range(1, ploidy + 1):
                samples.append(f"{name}.{i}")
    return samples


def write_to_file(file_path: str, sample_list: list[str]) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        for name in sample_list:
            f.write(f"{name}\n")


def turn_sampleID_standard():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("ploidy", type=int, default=2)

    args = parser.parse_args()

    sample_list = load_samples_with_suffix(args.input, args.ploidy)
    write_to_file(args.output, sample_list)


if __name__ == "__main__":
    turn_sampleID_standard()
