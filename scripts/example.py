"""example.py: A simple example script for an API call"""

from SimPG import *
from SimPG.SPexpe import *

def run_example():
    set_default_logging()  # Configuring Logging
    # Parse GFA files and BED files
    myGFA = Minigfa("./Pangenome.gfa")
    myBED = Minibed("./Pangenome.bed")
    # Construction of pangenome Digraph
    myDigraph = turn_GFA_to_DiGraph(myGFA, myBED)
    # Extract the walking route of individual genome sequence mapping in the graph
    simulate_population_every_walk(myGFA, myBED, myDigraph, "./region_sample.txt")
    # Get the core sequence nodes of the crowd
    core_seg = get_coreSeg_in_Pangenome(myGFA)
    # Assemble a pangenome graph for a specific population
    pangenome_graph = simulate_Population_Pangenome(myBED)
    # Simulate the genomes of individuals in a population
    simulate_Whole_Genome_Sequencing_for_population(pangenome_graph, myGFA, core_seg)


if __name__ == "__main__":
    run_example()
