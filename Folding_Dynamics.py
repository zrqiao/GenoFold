import shlex, subprocess
from difflib import SequenceMatcher
import numpy as np
import Domains
from collections import defaultdict
import nupack_functions
import argparse, math, random, gzip, pickle, types
from multiprocessing import Pool
import copy

# Change following routines for other environments:
L_init = 5  # Initiation unit
dL = 5  # elongation unit (also means CG unit)
dt = 100  # Folding time for each elongation step
population_size_limit = 25  # maximum type of strands in the pool
MULTI_PROCESS = 32

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline()

    #   NOTE: Initiation [create active population]
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
    step = 0
    # print('Population size: '+str(active_species_pool.size))

    # Start IO
    checkpoint_pool = gzip.open(clargs.sequence + '_pool.p.gz', 'w')
    pickle.dump(active_species_pool, checkpoint_pool)
    structure_output = open(clargs.sequence + '.dat', 'w')
    structure_output.write("#Time %g\n" % (dt * step))
    for domain in active_species_pool.species_list():
        structure_output.write('%s    %g\n' % (domain, active_species_pool.get_population(domain)))
    log = open(clargs.sequence + '.log', 'w')
    log.write('Step: %3d \n' % step)

    while sequence_length > current_length:
        step += 1
        log.write('Step: %3d \n'%step)
        
        old_species_list = copy.deepcopy(active_species_pool.species_list())
        old_species_pool = copy.deepcopy(active_species_pool)
        active_species_pool.clear()

        if current_length+dL > sequence_length:
            L_step = sequence_length-current_length
        else:
            L_step = dL
        current_length += L_step

        # Generate all new foldons
        l_bounds = np.arange(0, current_length, dL)
        # multi_pool = Pool(MULTI_PROCESS)
        for l_bound in l_bounds:
            all_foldons.new_foldon(full_sequence[l_bound:current_length], l_bound, current_length, all_domains)

               # print(old_species_list)

        # NOTE: population is inherited without updating its IFR!! No new domain instance should be created.

        # NOTE: structure_generation(single strain, elongation segment) [to be called in pool.map()]
        # Compute all IFR segments; link sequences; update IFRs
        for strand in old_species_list:
            Domains.recombination(
                strand, current_length, all_foldons,
                all_domains, old_species_pool, active_species_pool)  # parallelization

        # NOTE: population dynamics (master equation)

        print('Population size: '+str(active_species_pool.size))
        active_species_pool.evolution(all_pathways, dt)
        active_species_pool.selection(population_size_limit)

        # pickle & outputs
        # with gzip.open(clargs.sequence + '_pool.p.gz', 'a') as checkpoint_pool:
        pickle.dump(active_species_pool, checkpoint_pool)
        # with open(clargs.sequence + '.dat', 'a') as structure_output:
        structure_output.write("#Time %g\n" % (dt * step))
        for domain in active_species_pool.species_list():
            structure_output.write('%s    %g\n' % (domain, active_species_pool.get_population(domain)))
        # with open(clargs.sequence + '.log', 'a') as log:
        log.write('Step: %3d \n' % step)

        # pickle.dump(active_species_pool, checkpoint_pool)
        # structure_output.write('#Time %g\n'%(dt*step))
        # for domain in active_species_pool.species_list():
        #     structure_output.write('%s    %g\n'%(domain, active_species_pool.get_population(domain)))

    with gzip.open(clargs.sequence + '_domains.p.gz', 'w') as checkpoint_domains:
        pickle.dump(all_domains, checkpoint_domains)

    checkpoint_pool.close()
    log.close()
    structure_output.close()

exit()
