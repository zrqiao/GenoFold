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
L_init = 50  # Initiation unit
dL = 25  # elongation unit (also means CG unit)
dt = 0.1 * dL  # Folding time for each elongation step (0.1 s/nt)
population_size_limit = 10  # maximum type of strands in the pool
MULTI_PROCESS = 32

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    parser.add_argument('--k', type=np.float, default=1., \
                        help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '_k' + '%e' % clargs.k + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    #   NOTE: Initiation [create active population]
    all_domains = Domains.DomainsCollection()
    all_foldons = Domains.FoldonCollection()
    all_pathways = Domains.Pathways(clargs.k)
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
    checkpoint_pool = gzip.open(clargs.sequence + '_k' + '%e' % clargs.k + '_pool.p.gz', 'w')
    pickle.dump(active_species_pool, checkpoint_pool)
    structure_output = open(clargs.sequence + '_k' + '%e' % clargs.k + '.dat', 'w+')
    structure_output.write("#Time %g\n" % (dt * step))
    for domain in active_species_pool.species_list():
        structure_output.write('%s    %g\n' % (domain, active_species_pool.get_population(domain)))
    log = open(clargs.sequence + '_k' + '%e' % clargs.k + '.log', 'w+')
    log.write('Step: %3d \n' % step)
    log.flush()
    while sequence_length > current_length:
        step += 1
        # print('Step: %3d \n'%step)
        log.write('Step: %3d \n'%step)
        old_species_pool = copy.deepcopy(active_species_pool)
        old_species_list = old_species_pool.species_list()
        active_species_pool.clear()
        log.flush()

        if current_length+dL > sequence_length:
            L_step = sequence_length-current_length
        else:
            L_step = dL
        current_length += L_step
        log.write('Current length: %d \n'%current_length)
 
        # Generate all new foldons
        l_bounds = np.arange(0, current_length, dL)
        # multi_pool = Pool(MULTI_PROCESS)
        log.write('Calculate new foldons...')
        log.flush()
        for l_bound in l_bounds:
            all_foldons.new_foldon(full_sequence[l_bound:current_length], l_bound, current_length, all_domains)

               # print(old_species_list)
        log.write(' finished\n')
        log.flush()
        # NOTE: population is inherited without updating its IFR!! No new domain instance should be created.

        # NOTE: structure_generation(single strain, elongation segment) [to be called in pool.map()]
        # Compute all IFR segments; link sequences; update IFRs
        log.write('Generating new secondary structures...')
        log.flush()
        for strand in old_species_list:
            Domains.recombination(
                strand, current_length, all_foldons,
                all_domains, old_species_pool, active_species_pool)  # parallelization
        log.write(' finished\n')        

        log.write('active space: ')
        log.write(str(active_species_pool.species_list())+'\n')
        # NOTE: population dynamics (master equation)

        log.write('Population evolution... \n')
        log.write('Population size before selection: '+str(active_species_pool.size)+'\n')
        active_species_pool.evolution(all_pathways, dt)
        active_species_pool.selection(population_size_limit)
        log.flush()
        log.write('Time: %d \n'%active_species_pool.timestamp )
        log.write('Population size after selection: '+str(active_species_pool.size)+'\n')
        log.write('Selection finished \n')
        # pickle & outputs
        # with gzip.open(clargs.sequence + '_pool.p.gz', 'a') as checkpoint_pool:
        pickle.dump(active_species_pool, checkpoint_pool)
        # with open(clargs.sequence + '.dat', 'a') as structure_output:
        structure_output.write("#Time %g\n" % (dt * step))
        for domain in active_species_pool.species_list():
            structure_output.write('%s    %g\n' % (domain, active_species_pool.get_population(domain)))
        # with open(clargs.sequence + '.log', 'a') as log:
        # log.write('Step: %3d \n' % step)

        log.flush()
        structure_output.flush()
        # pickle.dump(active_species_pool, checkpoint_pool)
        # structure_output.write('#Time %g\n'%(dt*step))
        # for domain in active_species_pool.species_list():
        #     structure_output.write('%s    %g\n'%(domain, active_species_pool.get_population(domain)))

    with gzip.open(clargs.sequence + '_k' + '%e' % clargs.k + '_domains.p.gz', 'w') as checkpoint_domains:
        pickle.dump(all_domains, checkpoint_domains)

    checkpoint_pool.close()
    log.close()
    structure_output.close()

exit()
