# API Reference

This section describes in detail the public interfaces (classes/functions) of each part under `SimPG`.

---

## Data Structures

```python
from SimPG import Minigfa,Minibed
```

---

### 1. Class:  Minigfa

```python
class Minigfa:
	def __init__(self, file_path:Optional[str]=None) -> None:
        """
         The path of the GFA file that is preferably passed in when constructing the object.
         If you don't do this, you will just get an empty object. Please call the build_Minigfa method to construct
        """
		 ...
```

- **Description**
  Composite data storing GFA file information.
  Notice: Only lines S and L can be processed, lines starting with other letters are discarded

- **Methods**

  | Method                                                                       | Description                                                                                                                        |
  | ---------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
  | `build_Minigfa(self, file_path: str) -> None`                                | Constructor. If the file path is not passed in when creating the object, this method should be called.                             |
  | `get_linear_reference(self) -> str`                                          | Selector. Return the name of the pan-genome linear reference genome.                                                               |
  | `get_seq(self, segID: str) -> str`                                           | Selector. Return the sequence corresponding to segment ID.                                                                         |
  | `get_source_sample(self, segID: str) -> str`                                 | Selector. Return the name of stable sequence sample name from which the segment is derived corresponding to segment ID.            |
  | `get_SRank(self, *segI: str) -> int`                                         | Selector. Return SR corresponding to segment ID.                                                                                   |
  | `get_all_segID(self) -> Generator[str, Any, None]`                           | Provide a generator for iteration. Return a segment ID each time.                                                                  |
  | `get_all_Link(self) -> Generator[tuple[str, str, str, str, int], Any, None]` | Provide a generator for iteration. Return a five-tuple, fromID, fromOrient, toID,toOrient, SRank in order from a `Link` each time. |

- **Example**

  ```python
  from SimPG import Minigfa
  
  example_GFA = Minigfa()
  example_GFA.build_Minigfa("./pangenome.gfa")
  print(example_GFA.get_linear_reference())	# "CHM13"
  print(example_GFA.get_SRank("s15623"))	# 3
  for fromID, fromOrient, _, _, SR in example_GFA.get_All_Link():
      ...
  ```

---

### 2. Class:  Minibed

```python
class Minibed:
    def __init__(self, file_path: Optional[str] = None) -> None:
        """
        The path of the BED file that is preferably passed in when constructing the object.
        If you don't do this, you will just get an empty object. Please call the build_Minigfa method to construct.
        """
        ...
```

- **Description**
  Construct a composite data storing BED file information, and this class is iterable.
- **Methods**

| Methods                                                                         | Description                                                                                                                                                                                                                                             |
| ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `build_Minibed(self, file_path: str) -> None`                                   | Constructor. If the file path is not passed in when creating the object, this method should be called.                                                                                                                                                  |
| `__iter__(self) -> Generator[tuple[str, bool, int, int, list[str]], Any, None]` | Provide a generator for iteration. Return a five-tuple, chr_num, is_invered, segs_num, possible_path_num,list_of_segments in order from one line in BED file each time                                                                                  |
| `get_linear_sources_and_sinks(self) -> tuple[dict[str, str], dict[str, str]]`   | Selector. It returns the start and end nodes of each chromosome on the linear reference genome.The first dictionary of the tuple is all the starting node, and the second dictionary is all the ending node, expressed in the form of: {chr: segmentID} |

- Example
  ```python
  from SimPG import Minibed
  
  example_BED = Minibed("./pangenome.bed")
  for chr_num, is_inversed, segs_num, possible_path_num, list_of_segments in example_BED :
      ...
  ```

---

## Algorithm Functions

```python
# This is core pipeline function
from SimPG import run_SimPG
```

### Core Pipeline Function : run_SimPG

 ```python
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
     ...
 ```

- **Description**

  ​	Complete the entire simulation pipeline.

- **Args**

  ​	`GFA_file_path` (`str`) : Location of the GFA file to be processed.

  ​	`BED_file_path` (`str`) : Location of the GFA file to be processed.

  ​	`population` (`list[str] | str`) : Input a list of sample names, or a text file with only one sample name per line.

  ​	`enable_to_save_temporary_folder` (`bool`, optional) : Whether to keep intermediate files as `pickle`. Defaults to False.

  ​	`every_sample_Whole_Genome_Sequencing_filepath` (`Optional[str]`, optional) : The file location of the walking route of each sample. The default is the `my_walks.pl` file in the `/tmp` folder of the working directory. 

  ​	`whether_to_output_graph_information_in_terminal` (`bool`, optional) : Whether to output the key parameters of the graph to `stdout`. Defaults to False.

  ​	`sim_file_out_folder` (`Optional[str]`, optional) : Result output location. By default, the output is in the working directory.

  ​	`is_human` (`bool`, optional) : Is the pan-genome a human pan-genome? Defaults to `False`.

  ​	`population_name` (`Optional[str]`, optional) : Give your simulated crowd a name, which will also be used as the prefix for the output files. Defaults to "My".

  ​	`sim_num` (`int`, optional) : Number of simulations. Defaults to 1.

  ​	`logging_verbose` (`bool`, optional) : Whether to set the log output information level to at least `INFO` level .Default to `False`, set to `Warning` level . 


---

```python
# These are the steps functions that make up the core pipeline
from SimPG import turn_GFA_to_DiGraph, simulate_population_every_walk, simulate_Population_Pangenome, get_coreSeg_in_Pangenome, simulate_Whole_Genome_Sequencing_for_population
import networkx as nx
```

### 1. Function:  turn_GFA_to_DiGraph

```python
def turn_GFA_to_DiGraph(
    gfa_message: Minigfa,
    bed_message: Optional[Minibed] = None,
    is_output_inspection_results_in_graph: bool = False,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> nx.DiGraph:
    ...
```

- **Description**
  
  ​	Convert the GFA file information and Bed file information (optional, if default, you may find some loops or paths that should not exist in your graph) into a directed graph.

  ​	You can choose whether to output the key parameter information of the graph and whether to save it. 

  ​	Notice: It will generate a tmp folder storage under your working folder

- **Args**

  ​	`gfa_maessage` (`Minigfa`) : Composite data storing GFA file information

  ​	`bed_message` (`Minibed | None`, optional) : Composite data storing Bed file information. Defaults to `None`.

  ​	`is_output_inspection_results` (`bool`, optional) : Whether to output the key parameters of the graph to `stdout`. Defaults to `False`.

  ​	`is_saved_as_pickle` (`bool`, optional) : Whether to save as a pickle file for reuse. Defaults to `False`.

  ​	`file_path` (`str | None`, optional) : If you choose to save as a pickle file,the graph will be saved in `file_path`. By default, the file name will be `myMinigraph.pl` in folder /tmp under your working folder.

  

- **Returns**
  
  `DiGraph`: Directed graph representing the pan-genome. We think direction is 5' end to 3' end as you follow the diagram,and the ID of each node is a tuple, the first element is the ID of the segment in GFA, and the second element is the symbol "+" or "-". "-" represents the reverse complementary sequence of the connected segment sequence

---

### 2. Function:  simulate_population_every_walk

```python
def simulate_population_every_walk(
    gfa_message: Minigfa,
    bed_message: Minibed,
    G_full: nx.DiGraph,
    population: list[str] | str,
    saved_file_path: Optional[str] = None,
) -> None:
    ...
```

- **Description**

  ​	Without the need for original individual genome sequence information involved in building a pan-genome, this function can extract the path of individual genome sequences mapped in the graph.

  ​	Notice: This function does not return anything. It will save the walking route of each sample in `saved_file_path` file.

- **Args**

  ​	`gfa_message` (`Minigfa`) : Composite data storing GFA file information.

  ​	`bed_message` (`Minibed`) : Composite data storing Bed file information.

  ​	`G_full` (`nx.DiGraph`) : Pan-genome graph

  ​	`population` (`list[str] | str`) : Input a list of sample names, or a text file with only one sample name per line

  ​	`saved_file_path` (`str`) : Save file location.By default, it is saved in `my_walks.pl` in the `/tmp` folder of the working directory.

---

### 3. Function:  simulate_population_Pangenome

```python
def simulate_Population_Pangenome(
    bed_message: Minibed,
    every_sample_Whole_Genome_Sequencing_filepath: Optional[str] = None,
    is_added_linear_reference_genome: bool = False,
    is_output_inspection_results_in_graph: bool = False,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> nx.DiGraph:
```

- **Description**

  ​	Simulate the pan-genome of a specific population.

- **Args**

  ​	`bed_message` (`Minibed`) : Composite data storing Bed file information.

  ​	`every_sample_Whole_Genome_Sequencing_filepath` (`str | None`, optional) : The file location of the walking route of each sample. The default is the `my_walks.pl` file in the `/tmp` folder of the working directory

  ​	`is_added_linear_reference_genome` (`bool`, optional) : Whether to add a linear reference genome in new pan-genome graph. Defaults to `False`.

  ​	`is_output_inspection_results` (`bool`, optional) : Whether to output the key parameters of the graph to `stdout`. Defaults to `False`.

  ​	`is_saved_as_pickle` (`bool`, optional) : Whether to save as a pickle file for reuse. Defaults to `False`. 

  ​	`file_path` (`str | None`, optional) : If you choose to save as a pickle file, the graph will be saved in `file_path`. By default, the file name will be `myPangenome.pl` in folder `/tmp` under your working folder.

- **Returns**

  ​	`nx.DiGraph` : New pan-genome graph

---

### 4. Function:  get_coreSeg_in_Pangenome

```python
def get_coreSeg_in_Pangenome(
    gfa_message: Minigfa,
    every_sample_Whole_Genome_Sequencing_filepath: Optional[str] = None,
    is_saved_as_pickle: bool = False,
    file_path: Optional[str] = None,
) -> set[tuple[str, str]]:
```

- **Description**

  ​	Get the core sequence nodes  belonging to a certain group of people.

- **Args**

  ​	`gfa_message` (`Minigfa`) : Composite data storing GFA file information.

  ​	`every_sample_Whole_Genome_Sequencing_filepath` (`str | None`, optional) : The file location of the walking route of each sample. The default is the my_walks.pl file in the `/tmp` folder of the working directory.

  ​	`is_saved_as_pickle` (`bool`, optional) : Whether to save as a pickle file for reuse. Defaults to `False`.

  ​	`file_path` (`str | None`, optional) : If you choose to save as a pickle file, the graph will be saved in `file_path`. By default, the file name will be `myCoreseg.pl` in folder `/tmp` under your working folder.

- **Returns**

  ​	`set[tuple[str, str]]` : Save the collection of core sequence nodes.

---

### 5. Function:  simulate_Whole_Genome_Sequencing_for_population

```python
def simulate_Whole_Genome_Sequencing_for_population(
    Pangenome_graph: nx.DiGraph,
    gfa_message: Minigfa,
    coreSeg: set[tuple[str, str]],
    file_out_folder: Optional[str] = None,
    is_human: bool = False,
    population_name: Optional[str] = None,
    sim_num: int = 1,
) -> None:
```

- **Description**

  ​	Without using the method of pre-setting and saving the weight matrix, directly randomly select the next node to walk.

  ​	Notice: This function will generate two folders under `file_out`, namely `{population_name}_simulate_fasta` and `{population_name}_simulate_rvcf`. The `fasta` files and `rvcf` files will be saved in the following folders respectively.

- **Args**

  ​	`Pangenome_graph` (`nx.DiGraph`) : Pan-genome graph.

  ​	`gfa_message` (`Minigfa`) : Composite data storing GFA file information.

  ​	`coreSeg` (`set[tuple[str, str]]`) : A collection of core sequence nodes.

  ​	`file_out_folder` (`str | None`, optional) : Result output location. By default, the output is in the working directory.

  ​	`is_human` (`bool`, optional) : Is the pan-genome a human pan-genome? Defaults to `False`.

  ​	`population_name` (`str | None`, optional) : Give your simulated crowd a name, which will also be used as the prefix for the output files. Defaults to "My".

  ​	`sim_num` (`int`, optional) : Number of simulations. Defaults to `1`.

- **Raises**

  ​	`ValueError` : The start or end node is not in the graph.

  ​	`networkx.NetworkXNoPath` : There is no path from the start node to the end node.

---

## Additional callable helper functions

```python
from SimPG import sim_part, sim_part_for_num, set_default_logging
```

---

### 1. Function:  sim_part

```python
def sim_part(
    in_rvcf: str,
    out_fasta: str,
    out_rvcf: str,
    bed_message: Minibed,
    gfa_message: Minigfa,
    fraction: float,
    is_human: bool = False,
) -> None:
    ...
```

- **Description**

  ​	Keep the variation of fraction ratio in `in_vcf` file, and output the corresponding new `fasta` file & `rvcf` (unconverted vcf file) file.

- **Args**

  ​	`in_rvcf` (`str`) : Input `rvcf` (unconverted `vcf` file) file location.

  ​	`out_fasta` (`str`) : Output `fasta` file location.

  ​	`out_rvcf` (`str`) : Output `rvcf` (unconverted `vcf` file) file location.

  ​	`bed_message` (`Minibed`) : Composite data storing BED file information.

  ​	`gfa_message` (`Minigfa`) : Composite data storing GFA file information.

  ​	`fraction` (`float`) : The proportion of the number of retained variants.

  ​	`is_human` (`bool`, optional) : Is the pan-genome a human pan-genome? Defaults to `False`.

- **Raises**
  
    `ValueError` : The scale must be between 0.0 and 1.0. But the input parameters do not meet the requirements.

---

### 2. Function:  sim_part_for_num

```python
def sim_part_for_num(
    in_rvcf_folder:str,
    out_folder:str,
    bed_message: Minibed,
    gfa_message: Minigfa,
    population_name: str,
    fraction: float,
    is_human: bool = False,
    sim_num=10,
) -> None:
    ...
```

- **Description**

  ​	Perform `sim_part` processing on the `num` rvcf files in the `in_rvcf_folder` folder (the folder contains the rvcf of a population).

- **Args**

  ​	`in_rvcf_folder` (`str`) : Input `rvcf`(unconverted vcf file) file folder location.

  ​	`out_folder` (`str`) : Output file location(will generate two folder under this folder: `fasta` folder and `rvcf` folder).

  ​	`bed_message` (`Minibed`) : Composite data storing BED file information.

  ​	`gfa_message` (`Minigfa`) : Composite data storing GFA file information.

  ​	`population_name` (`str`) : `in_vcf_folder` 's Population name.

  ​	`fraction` (`float`) : The proportion of the number of retained variants.

  ​	`is_human` (`bool`, optional) : Is the pan-genome a human pan-genome? Defaults to `False`.

  ​	`num` (`int`, optional) : Number of `rvcf` files. The default is ten times.

- **Raises**

  ​	`ValueErro` r: The scale must be between 0.0 and 1.0. But the input parameters do not meet the requirements.

---

### 3. Function:  set_default_logging

```python
def set_default_logging(verbose: bool = False) -> None:
    ...
```

- **Description**

  ​	Set up logging for users.

- **Args**

  ​	`verbose` ( bool , optional) : Whether to set the log output information level to at least `INFO` level .Default to `False`, set to `Warning` level . 