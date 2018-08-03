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
from multiprocessing import Pool

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

        old_species_pool = active_species_pool
        active_species_pool.clear()

        if current_length+dL > sequence_length:
            L_step = sequence_length-current_length
        else:
            L_step = dL
        current_length += L_step

        # Generate all new foldons
        l_bounds = np.arange(0, current_length, dL)
        multi_pool = Pool()
        multi_pool.map(lambda l_bound: all_foldons.new_foldon(
            full_sequence[l_bound:current_length], l_bound, current_length, all_domains), l_bounds)

        old_species_list = old_species_pool.species_list()

        # NOTE: population is inherited without updating its IFR!! No new domain instance should be created.

        # NOTE: structure_generation(single strain, elongation segment) [to be called in pool.map()]
        # Compute all IFR segments; link sequences; update IFRs
        for old_species in old_species_list:    #  TODO: Need parallel
            for terminal_foldon in all_foldons.find_foldons(old_species.r_bound, current_length):
                active_species_pool.add_species(old_species.elongate(terminal_foldon))
            for rearrange_point in reversed(old_species.get_IFR()):
                for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                    shortersegment = ... #TODO
                    active_species_pool.add_species(old_species.elongate(overlapping_foldon))
                active_species_pool.add_species()

        # TODO: active pool update

        # TODO: population_selection (need a active population, fitness function)

        active_species_pool.evolution(all_pathways, dt)
        active_species_pool.selection(population_size_limit)
        # NOTE: compute_foldon(i,j)

        # TODO: pickle & outputs



exit()
