import gzip
import datetime
import ast

from collections import defaultdict


def extract_chr22_from_vcf(input_vcf, output_vcf):
    """
    Extracts lines with chromosome 'chr22' from a VCF file and writes to an output file.

    Parameters:
        input_vcf (str): Path to the input VCF file (can be .vcf or .vcf.gz).
        output_vcf (str): Path to the output VCF file.
    """
    try:
        # Open the input file (supports gzip)
        if input_vcf.endswith(".gz"):
            input_file = gzip.open(input_vcf, "rt")  # Open as text
        else:
            input_file = open(input_vcf, "r")

        # Open the output file
        with open(output_vcf, "w") as output_file:
            for line in input_file:
                # Write header lines directly to the output
                if line.startswith("#"):
                    output_file.write(line)
                else:
                    # Split the line into columns to check the chromosome field
                    fields = line.strip().split("\t")
                    if fields[0] == "chr22":  # Chromosome is in the first column
                        output_file.write(line)

        print(f"Filtered lines with 'chr22' written to {output_vcf}")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        input_file.close()


class pan_vcf_ope:
    def __init__(self) -> None:
        pass

    def extract_AT_info(self, input_vcf):
        var_dic = {}
        """
        Extracts 'AT' annotation from the INFO field in a VCF file.

        Parameters:
            input_vcf (str): Path to the input VCF file (can be .vcf or .vcf.gz).
            output_file (str): Path to the output file where extracted AT information will be saved.
        """
        try:
            # Open the input VCF file (supports gzip)
            if input_vcf.endswith(".gz"):
                input_file = gzip.open(input_vcf, "rt")  # Open in text mode
            else:
                input_file = open(input_vcf, "r")

            # Open the output file for writing
            for line in input_file:
                if line.startswith("#"):  # Skip header lines
                    continue

                # Split the VCF line into fields
                fields = line.strip().split("\t")

                chrom = fields[0]  # Chromosome
                pos = fields[1]  # Position
                var_id = fields[2]  # Variant ID
                info = fields[7]  # INFO field

                # Extract the AT field from the INFO column
                info_dict = {
                    item.split("=")[0]: item.split("=")[1]
                    for item in info.split(";")
                    if "=" in item
                }
                at_info = info_dict.get(
                    "AT", "NA"
                )  # Get AT value, or 'NA' if not present

                if "_" in var_id:
                    # 如果有下划线，分割并取下划线前的部分
                    var_id = var_id.split("_")[0]

                id1, va1 = self.parse_variants(var_id, at_info)
                # Write the extracted information to the output file
                # print(f"{chrom}\t{pos}\t{var_id}\t{at_info}\n")
                # print(f"{tuple(id1)}\t{va1}\n")
                if tuple(id1) in var_dic:
                    var_dic[tuple(id1)] += va1
                else:
                    var_dic[tuple(id1)] = va1

        except Exception as e:
            print(f"An error occurred: {e}")

        finally:
            input_file.close()

        return var_dic

    def parse_variants(self, main_id, at_string):
        """
        Parses main ID and AT string into structured lists.

        Parameters:
            main_id (str): Main variant ID string (e.g., '>44162124>44162127').
            at_string (str): AT annotation string (e.g., '>44162124>44162125>44162127,>44162124<44162126>44162127').

        Returns:
            tuple: A tuple containing:
                - main_list: List of integers extracted from main_id.
                - at_list: List of lists, where each sublist represents a parsed variant from at_string.
        """
        # Parse main_id into a list
        # main_list = [int(item) for item in main_id.split('>') if item]
        split_char = None
        if "<" in main_id and ">" in main_id:
            raise ValueError(
                "输入的字符串中不能同时包含 < 和 > 作为分割字符，请确认输入内容"
            )
        elif "<" in main_id:
            split_char = "<"
        elif ">" in main_id:
            split_char = ">"
        else:
            raise ValueError("输入的字符串中未找到有效的分割字符，请确认输入内容")

        parts = main_id.split(split_char)
        if split_char == "<":
            sub_list = [int(num) for num in parts if num]
            main_list = sub_list[::-1]
        else:
            main_list = [int(item) for item in parts if item]

        # Parse AT string into a nested list
        at_variants = at_string.split(",")
        at_list = []
        for variant in at_variants:
            sublist = []
            split_numbers = variant.replace("<", ",-").replace(">", ",").split(",")
            sublist = [int(num) for num in split_numbers if num]

            # for part in variant.split('>'):

            #     if '<' in part:  # Handle negative values
            #         sub_parts = part.split('<')
            #         sublist.append(int(sub_parts[0]))
            #         sublist.append(-int(sub_parts[1]))
            #     elif part:  # Handle regular positive values
            #         sublist.append(int(part))
            at_list.append(sublist)

        return main_list, at_list


def read_sim_answer_txt(file_path):
    sim_ans_dic = {}
    with open(file_path, "r") as file:
        lin_count = 0
        for line in file:
            li = line.strip().split("\t")
            sim_ans_dic[eval(li[0])] = eval(li[2])
            lin_count += 1
    return sim_ans_dic, lin_count


def S_position_build_v1_0(file_path):
    seq_dic = {}
    pos_dic = {}
    chr_dic = {}
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()  # 读取文件的所有行
            # lines_starting_with_S = [line.strip() for line in lines if line.strip().startswith('S')]  # 找到以'L'开头的行并去除空格和换行符
            for line in lines:
                temp_s = line.strip()
                if temp_s.startswith("S"):
                    easy_line = line.strip().split("\t")
                    # if len(easy_line)>= 6:
                    so_i_pos = easy_line[5]
                    so_i_pos1 = so_i_pos.split(":")
                    position = int(so_i_pos1[2])
                    # sequence = easy_line[2]
                    pos_dic[easy_line[1]] = position

                    chr_pso = easy_line[4]
                    chr_pso = chr_pso.split(":")
                    chr_dic[easy_line[1]] = chr_pso[2]

                    seq_dic[easy_line[1]] = easy_line[2]

            return seq_dic, pos_dic, chr_dic
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []


def read_sim_answer_txt_to_make_vcf_v1_0(file_path):
    def parse_custom_format(s):
        # 移除括号和空格
        s = s.strip().strip("()").replace(" ", "")
        # 按逗号分隔
        items = s.split(",")
        # 去除末尾 + 或 -
        cleaned = [str(item) for item in items]
        return tuple(cleaned)  # 返回元组

    sim_ans_dic = {}
    ref_ans_dic = {}
    with open(file_path, "r") as file:
        lin_count = 0
        for line in file:
            li = line.strip().split("\t")
            ref_ans_dic[parse_custom_format(li[0])] = parse_custom_format(li[1])
            sim_ans_dic[parse_custom_format(li[0])] = parse_custom_format(li[2])
            lin_count += 1
    return ref_ans_dic, sim_ans_dic, lin_count


def reverse_complement(sequence):
    # 创建一个字典映射每个碱基到它的互补碱基
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    # 生成互补链，并反转序列
    reverse_complement_seq = "".join(complement[base] for base in reversed(sequence))

    return reverse_complement_seq


import argparse


def main_true_vcf(sim_answer, sim_vcf, sample_name, source_data):

    # sim_sv_vcf
    # 读入txt文件，分别存储为参考和变异
    ref_ans_dic, sim_ans_dic, sim_var_num = read_sim_answer_txt_to_make_vcf_v1_0(
        sim_answer
    )
    # 提取节点序列
    seq_dic, pos_dic, chr_dic = S_position_build_v1_0(source_data)
    # 新增记录节点染色体号
    ref_seq = ""
    sim_seq = ""
    with open(sim_vcf, "w") as file1:
        # header = ["##fileformat=VCFv4.2","##fileDate=" + datetime.datetime.now().strftime("%Y%m%d"),"##source=ExampleVCFGenerator","##contig=<ID=chr22,length=50818468>",'##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">','##ALT=<ID=DEL,Description="Deletion relative to the reference">','##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">','##ALT=<ID=INV,Description="Inversion of reference sequence">','##ALT=<ID=BND,Description="Breakend of translocation">','##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">','##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">','##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">','##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">','##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">','##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">','##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">','##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends">','##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">','##INFO=<ID=STRAND,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">','##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names of SVs (comma separated)">','##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">','##FILTER=<ID=q5,Description="Quality below 5">','##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">','##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# High-quality reference reads">','##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">','##FORMAT=<ID=PL,Number=G,Type=Integer,Description="# Phred-scaled genotype likelihoods rounded to the closest integer">','##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality">',"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]

        header = [
            "##fileformat=VCFv4.3",
            "##fileDate=" + datetime.datetime.now().strftime("%Y%m%d"),
            "##source=syri",
            "##contig=<ID=chr1,length=248387328>",
            "##contig=<ID=chr2,length=242696752>",
            "##contig=<ID=chr3,length=201105948>",
            "##contig=<ID=chr4,length=193574945>",
            "##contig=<ID=chr5,length=182045439>",
            "##contig=<ID=chr6,length=172126628>",
            "##contig=<ID=chr7,length=160567428>",
            "##contig=<ID=chr8,length=146259331>",
            "##contig=<ID=chr9,length=150617247>",
            "##contig=<ID=chr10,length=134758134>",
            "##contig=<ID=chr11,length=135127769>",
            "##contig=<ID=chr12,length=133324548>",
            "##contig=<ID=chr13,length=113566686>",
            "##contig=<ID=chr14,length=101161492>",
            "##contig=<ID=chr15,length=99753195>",
            "##contig=<ID=chr16,length=96330374>",
            "##contig=<ID=chr17,length=84276897>",
            "##contig=<ID=chr18,length=80542538>",
            "##contig=<ID=chr19,length=61707364>",
            "##contig=<ID=chr20,length=66210255>",
            "##contig=<ID=chr21,length=45090682>",
            "##contig=<ID=chr22,length=51324926>",
            "##contig=<ID=chrX,length=154259566>",
            "##contig=<ID=chrY,length=62460029>",
            '##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">',
            '##ALT=<ID=DEL,Description="Deletion relative to the reference">',
            '##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">',
            '##ALT=<ID=INV,Description="Inversion of reference sequence">',
            '##ALT=<ID=BND,Description="Breakend of translocation">',
            '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
            '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
            '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
            '##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends">',
            '##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">',
            '##INFO=<ID=STRAND,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">',
            '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names of SVs (comma separated)">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">',
            '##FILTER=<ID=q5,Description="Quality below 5">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# High-quality reference reads">',
            '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">',
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="# Phred-scaled genotype likelihoods rounded to the closest integer">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality"> ',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name,
        ]

        for line in header:
            file1.write(line + "\n")
        for key, values in sim_ans_dic.items():
            # pos_id = key[1]
            pos = int(pos_dic[key[0]]) + len(seq_dic[key[0]])
            pos = int(pos)
            # pos = pos_dic[key[1]]

            # for value in values:

            ref_seq = ""
            sim_seq = ""
            if len(ref_ans_dic[key]) == 2:
                ref_seq = ""
            else:
                for i in range(1, len(ref_ans_dic[key]) - 1):
                    ref_seq = ref_seq + seq_dic[ref_ans_dic[key][i][:-1]]

            if len(values) == 2:
                sim_seq = ""
            else:
                for i in range(1, len(values) - 1):
                    if values[i].endswith("-"):
                        sim_seq = sim_seq + reverse_complement(seq_dic[values[i][:-1]])
                    else:
                        sim_seq = sim_seq + seq_dic[values[i][:-1]]

            var_len = len(sim_seq) - len(ref_seq)

            if ref_seq == "":
                ref_seq = sim_seq[0]
            if sim_seq == "":
                sim_seq = ref_seq[0]

            if var_len > 0:
                var_type = "INS"

            else:
                var_type = "DEL"
            file1.write(
                f"{chr_dic[key[0]]}\t{str(pos)}\tSimPG.{var_type}\t{ref_seq}\t{sim_seq}\t.\tPASS\tPRECISE;SVTYPE={var_type};SVLEN={abs(var_len)};END={int(pos) + abs(var_len)}\tGT\t1\n"
            )


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rvcf_in", required=True)
    parser.add_argument("--vcf_out", required=True)
    parser.add_argument("--source_GFA", required=True)
    parser.add_argument("--sample_name", default="my_sim_answer")
    args = parser.parse_args()
    sim_answer = args.rvcf_in
    sim_vcf = args.vcf_out
    source_data = args.source_GFA
    sample_name = args.sample_name
    main_true_vcf(sim_answer, sim_vcf, sample_name, source_data)


if __name__ == "__main__":
    cli()
