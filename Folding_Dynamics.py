import shlex, subprocess
from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import Domains
from collections import defaultdict
import nupack_functions
import argparse, math, random, gzip, pickle, types

# Change following routines for other environments:
L_init = 20  # Initiation unit
dL = 5  # elongation unit (also means CG unit)
dt = 1  # Folding time for each elongation step
population_size_limit = 25  # maximum type of strands in the pool

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    clargs = parser.parse_args()
    with open(clargs.sequence, 'r') as sequence_file:
        full_sequence = sequence_file.readline()


    #NOTE: Initiation [create active population]
    all_domains = Domains.DomainsCollection()
    all_foldons = Domains.FoldonCollection()
    all_pathways = Domains.Pathways()
    active_species_pool = Domains.SpeciesPool(all_pathways)

    # NOTE: Initiate transcription
    init_segment = full_sequence[:L_init]
    init_foldon = all_foldons.new_foldon(init_segment, 0, L_init, all_domains)
    active_species_pool.add_species(init_foldon, population=1.0)
    sequence_length = len(full_sequence)
    current_length = L_init

    while sequence_length < current_length:
        if current_length+dL > sequence_length:
            L_step = sequence_length-current_length
        else:
            L_step = dL
        current_length += L_step
        # TODO: structure_generation(single strain, elongation segment) [to be called in pool.map()]
        #Compute all IFR segments; link sequences; update IFRs

        # TODO: active pool update

        # TODO: population_selection (need a active population, fitness function)

        active_species_pool.evolution(all_pathways, dt)
        active_species_pool.selection(population_size_limit)
        # NOTE: compute_foldon(i,j)

        # TODO: pickle & outputs

exit()