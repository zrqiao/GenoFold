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
import copy

def recombination(strand, all_foldons, all_domains, old_species_pool, active_species_pool):
			for terminal_foldon in all_foldons.find_foldons(strand.r_bound, current_length):
                active_species_pool.add_species(strand.elongate(terminal_foldon), 
				population = old_species_pool.get_population(strand)/ len(all_foldons.find_foldons(strand.r_bound, current_length)))
            for rearrange_point in reversed(strand.get_IFR()[:-1]):
				if rearrange_point == 0:#Global rearrangement
					active_species_pool.add_species(all_foldons.find_foldons(rearrange_point, current_length))
				else:
					for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
						unrearranged_domain = all_domains.get_domain(full_sequence[0:rearrange_point], strand.get_structure()[0:rearrange_point], 0, rearrange_point)
						active_species_pool.add_species(unrearranged_domain.elongate(overlapping_foldon))

# Change following routines for other environments:
L_init = 20  # Initiation unit
dL = 5  # elongation unit (also means CG unit)
dt = 1  # Folding time for each elongation step
population_size_limit = 25  # maximum type of strands in the pool

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline()

	checkpoint_pool = open(clargs.sequence + '_pool.p.tgz', 'w')
	structure_output = open(clargs.sequence + '.dat', 'w')
	log = open(clargs.sequence + '.log', 'w')
	
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
	pickle.dump(active_species_pool, checkpoint_pool)
	step = 0
	log.write('Step: %3d \n'%step)
	structure_output.write('#Time %g'%(dt*step))
	for domain in active_species_pool.species_list:
		structure_output.write('%s    %g'%(domain, active_species_pool.get_population(domain)))
	
    while sequence_length < current_length:
		step +=1
		log.write('Step: %3d \n'%step)
		
        old_species_pool = copy.deepcopy(active_species_pool)
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
			
		multi_pool.map(lambda strand: recombination(strand, all_foldons, 
			all_domains, old_species_pool, active_species_pool), old_species_list)	#  Parallelize

        # NOTE: population dynamics (master equation)

        active_species_pool.evolution(all_pathways, dt)
        active_species_pool.selection(population_size_limit)

        # pickle & outputs
		pickle.dump(active_species_pool, checkpoint_pool)
		structure_output.write('#Time %g'%(dt*step))
		for domain in active_species_pool.species_list:
			structure_output.write('%s    %g'%(domain, active_species_pool.get_population(domain)))
		
	with open(clargs.sequence + '_domains.p.tgz', 'w') as checkpoint_domains:
		pickle.dump(all_domains, checkpoint_domains)

	checkpoint_pool.close()
	log.close()
	structure_output.close()

	exit()
