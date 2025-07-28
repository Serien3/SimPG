## Prepare Materials

Since the original purpose of developing this tool was to develop a set of scripts that can process [Minigraph pangenome graphs for HPRC samples](https://zenodo.org/records/15252892), this tool has certain requirements for the format of the received files.

Basically, please prepare the following two files：

1. The pangenome graph in the [GFA format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), or preferrably the [rGFA format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md). 
2. Genome annotation  `BED` files. It is best obtained by calling structural variants using the [Minigraph](https://github.com/lh3/minigraph) tool.
3. (Optional) Sample name file for simulated population.

If you don't have a graph, you can generate a graph from multiple samples by [Minigraph](https://github.com/lh3/minigraph).

Based on the standard file format, the tool also has certain requirements for the specific format of the file.

**By the way, much of the following content is excerpted and referenced from Minigraph's documentation, and we would like to express our gratitude for that.**

---

### rGFA data format

GFA primarily consists of two types of lines: segment lines (S-lines) and link lines (L-lines). An S-line represents a sequence; an L-line represents a connection between two oriented segments. The core of GFA can be loosly described by the following grammar:

```
<gfa>    <- <line>+
<line>   <- <S-line> | <L-line>
<S-line> <- 'S' <segId> <seq> <tag>*
<L-line> <- 'L' <segId> <strand> <segId> <strand> <cigar> <tag>*
<strand> <- '+' | '-'
```



rGFA is a strict subset of GFA. It disallows overlaps between segments and requires three additional tags on each segment. These tags trace the origin of the segment:

| Tag  | Type | Description                                               |
| ---- | ---- | --------------------------------------------------------- |
| `SN` | `Z`  | Name of stable sequence from which the segment is derived |
| `SO` | `i`  | Offset on the stable sequence                             |
| `SR` | `i`  | Rank. `0` if on a linear reference genome; `>0` otherwise |

For  files accepted for this tool, we require that `SN` and `SR` tags must be present and meaningful. 

The `SN` tag must be represented as:

```php
<SN-tag>   <- "SN:Z:" <PanSN>
<PanSN>        <- <SampleID> "#" <HaplotypeID> "#" [<ContigID>]
<SampleID>     <- <Alnum>+ [ "." <Digit>+ ]
<HaplotypeID>  <- <Digit>+ 
<ContigID>     <- <Alnum> <AlnumDotUnder>* [ "." <Digit>+ ]

<Alnum>        <- "A" … "Z" | "a" … "z" | "0" … "9"
<Digit>        <- "0" … "9"
<AlnumDotUnder><- <Alnum> | "." | "_" 
```

**Notes:**

- `+` means "one or more";
- `*` means "zero or more";
- `[...]` indicates optional
- Literals are enclosed in double or single quotes.

Actually, in the HPRC pangenome graphs, the **`SN:Z:`** field uses the Pangenome Sequence Naming (PanSN) convention to pack `sampleID`, `haplotype`, and `contig identifiers` into one string. For example, `SN:Z:NA21110.2#2#JBHIJJ010000042.1` is meaningful. We understand it this way: 

- `NA21110.2`

  ​	The sample identifier (here, individual NA21110, assembly version 2).

- `#2`

  ​	Haplotype 2 (e.g. paternal vs. maternal phase group).

- `#JBHIJJ010000042.1`

  ​	The original contig or scaffold name from that assembly.

Putting it all together, that tag tells you “this segment comes from contig `JBHIJJ010000042.1` of haplotype `2` of sample `NA21110 (v2)`.”

It should be noted that our tool processes `rGFA` files based on two important assumptions:

1. Assembly version in sample identifier and Haplotype are equivalent. So you can think `SN:Z:NA21110.2#2#JBHIJJ010000042.1` and `SN:Z:NA21110.2##` are equivalent.  (**Notice: This is only for the case where SR>0**)
2. The sample identifier and SR are one-to-one corresponding. for example, all segments derived from `NA21110.2` have a `SR` of `199`. In other words, when you use `Minigraph` to build a pangenome graph, you need to ensure that the `SampleID` part in the header of the `fasta` file of any two assembly is not repeated.

For example, for `SN:Z:NA21110#2#    SR:i:199` and `SN:Z:NA21110#1#    SR:i:198`, SimPG will assume that they are the genomes of the same sample chromosome, though their `SR` are different. 

On the other hand, we will notice that there is some linear reference genome segments with `SR=0` in the pangenome graph. For their SN tags, `<ContigID>` must exist and is used to represent the chromosome name where the segment is located. For example, `CHM13#0#chr1` is acceptable.

---

### BED data format

We recommend using the `Calling structural variations` function provided by `Minigraph` to obtain the BED file directly. The data format may be as follows: 

- (1 - 3) The first three columns give the position of a bubble/variation
- (4) # GFA segments in the bubble including the source and the sink of the bubble
- (5) # all possible paths through the bubble (not all paths present in input samples)
- (6) 1 if the bubble involves an inversion; 0 otherwise
- (7) length of the shortest path (i.e. allele) through the bubble
- (8) length of the longest path/allele through the bubble
- (9-11) please ignore
- (12) list of segments in the bubble; first for the source and last for the sink
- (13) sequence of the shortest path (`*` if zero length)
- (14) sequence of the longest path (NB: it may not be present in the input samples)

SimPG requires that the user's BED file contain at least the (1)(6)(12) lines. The remaining columns can be replaced by `*`.

BTW, this tool requires that (1) is consistent with the `SN` tag content format in the GFA file. For example, (1) need like `CHM13#0#chr1`.

---

### SampleID file for simulated population

When calling the `run_SimPG` or `simulate_population_every_walk` function, the Arg `population` will be a list of sample names, or a text file with only one sample name per line if you prepare it. `SampleID file for simulated population` file's data format is similar to the following:

```
HG04184.1
HG04184.2
HG03831.1
HG03831.2
HG03927.1
HG03927.2
HG02738.1
...
```

Every line is `sampleI identifier`(Including assembly version). However, most of the time, you only know the name of the sample you want to extract, not the complete `SampleID`. In order to avoid users manually copying and entering the version number, we provide a [script](./scripts/turn_sampleID_standard.py) to help you convert.

If you have a file like this :

```
HG04184
HG03831
HG03927
HG02738
```

Each sample is diploid, and each chromosome of the sample is involved in the construction of the pangenome. Execute the script via command 

```bash 
$ python3 <filePath_in> <filePath_out> <sample_chromosome_ploidy>(Default is 2)
```

You will get a new file like: 

```
HG04184.1
HG04184.2
HG03831.1
HG03831.2
HG03927.1
HG03927.2
HG02738.1
HG02738.2
```

