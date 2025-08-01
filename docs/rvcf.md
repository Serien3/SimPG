## The Reference VCF（rvcf）Format

In order to facilitate the recording, analysis and verification of the simulated genome variation, We use a non-standard intermediate output file called `rvcf`. It uses segment ID in the pangenome as the index output. Although this makes `rvcf` completely dependent on the pangenome `GFA` file, it is very fast to generate and parse it by SimPG. 

Specifically, it has only three columns per row and is separated by `\t`. 

- The first column gives the interval of the variant region on the reference genome.
- The second column gives the segment path of the reference genome within the variant region.
- The third column gives the segment path for detecting genomic variation regions.

For example, 

```
(s144,s145)	(s144+,s145+)	(s144+,s1305104-,s145+)
```

Obviously, this indicates an insertional structural variation. Among them, "-" represents the reverse complementary sequence of the connected segment sequence。

If you want to convert the `rvcf` file into a standard `vcf` format data file that records sequence variations, we provide a [script](../scripts/rvcf_to_vcf.py) to help with this conversion. To use it, enter the following command in the command line :

```bash 
python3 rvcf_to_vcf.py --rvcf_in <rvcf_file_path>\
> --vcf_out <vcf_file_path>\
> --source_GFA <GFA_file_path>\
> [--sample_name <sample_name>]
```

**Note: <sample_name> defaults to "my_sim_answer"**