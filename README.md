# SimPG

![tip](https://badgen.net/badge/python/3.9/green?icon=github)
![tip](https://badgen.net/github/license/Serien3/SimPG)
![Last Update](https://img.shields.io/github/last-commit/Serien3/SimPG.svg?label=Last%20Update)

A population-specific haplotype genome simulation tool developed based on pangenome data.

## Table of contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Prepare Materials](#prepare-materials)
- [Getting Start](#getting-start)
- [Usage](#usage)
  - [API  Reference](#api--reference)
  - [CLI  Reference](#cli--reference)
- [Datasets generated from SimPG](#datasets-generated-from-simpg)
- [License](#license)
- [Contact](#contact)

## Introduction

The rapid advancement of high-throughput genome sequencing has enabled large-scale reconstruction of genome sequences at both individual and population levels. Existing linear genome simulation tools, typically based on a single reference genome, offer limited biological realism and fail to capture the genomic diversity and structural complexity needed for comprehensive evaluation. SimPG is a novel simulation tool that generates individual genomes with population-level characteristics by leveraging the rich variant and structural information embedded in pangenomes. SimPG produces realistic, high-quality simulated genomes that support diverse applications such as structural variant detection and population genetics research. 

## Installation

```bash
$ pip install git+https://github.com/Serien3/SimPG.git
-- OR --
$ git clone https://github.com/Serien3/SimPG.git && cd SimPG/ && pip install .
```

## Dependencies

```
This package temporarily only depends on the third-party library networkx at runtime.So it's very lightweight
We recommend using networkx==3.5.
```

## Prepare Materials

You should prepare at least two types of files :

1. The pangenome graph in the [GFA format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), or preferrably the [rGFA format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md). 
2. Genome annotation  `BED` files. It is best obtained by calling structural variants using the [Minigraph](https://github.com/lh3/minigraph) tool.
3. (Optional) Sample name file for simulated population.

Based on the standard file format, the tool also has certain requirements for the specific format of the file. So before using this tool, please refer to  [data formats](./docs/data-formats.md) to prepare the files accepted and then run the program.

## Getting Start

Before using this tool, please refer to the [data format](./docs/data-formats.md) to prepare the input file accepted. If you are sure your documents are acceptable and correctly installed `SimPG`, you can use the following python script to quickly start a simulation work.

```python
# start_quick.py
from SimPG import run_SimPG

run_SimPG("./Pangenome.gfa","./Pangenome.bed","./region_sample.txt")
```

**Notice: Before running, please replace the string parameters involving the file path with your own file path**.

After the program is finished running, you will see three more folders in the working directory where you ran the script. Among them, the `tmp` folder contains the `my_walks.pl` file, the `my_simulate_fa` folder contains the simulated `fasta` file, and the `my_simulate_rvcf` folder contains the simulated `rvcf` file. 

Among them, `rvcf` is a non-standard output file we created with the help of pan-genome to record mutations. For more information and how to standardize this file, please see [rvcf format](./docs/rvcf.md).

## Usage

### API  Reference

SimPG is not only a standardized process tools, but also a programming library. SimPG provides and maintains some python APIs. Full API reference documentation is available at [api reference](./docs/api.md) . SimPG aims to keep APIs in [SimPG.py](./src/SimPG/SimPG.py) and will ensure the stability of these APIs in the current version. The [SPexpe.py](./src//SimPG//SPexpe.py) file contains some experimental APIs, which may change frequently in the current version, but their general functions will not change.

The File [example.py](./scripts/example.py) demonstrates typical uses of python APIs. In fact, the effect of [example.py](./scripts/example.py) is the same as the script shown in [getting start](#getting-start).

### CLI  Reference

With the package, a command line tool called `SimPG` is also installed. It currently allows users to quickly perform a full-pipeline simulation of SimPG through the command line. Call it with `-h` or `--help` for help.

However, for some reasons, we cannot guarantee the stability and timely updates of the command line tools, even though we will always keep the most basic functions normal. Therefore, if you have more personalized needs or want stability, it is recommended to use APIs to write scripts to run.

## Datasets generated from SimPG

In [simulation](./simulation/), you can see the results and the whole process of simulation using the Minigraph pangenome graphs for HPRC samples (v4.0) as input data.

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for more information.

## Contact

For advising, bug reporting and requiring help, please post on [Github Issue](https://github.com/Serien3/SimPG/issues) or contact whhit1825@outlook.com





