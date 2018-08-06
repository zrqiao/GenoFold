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
Temperature = 37

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
    structure_output = open(clargs.sequence + '.dat', 'w+')
    structure_output.write("#Time %g\n" % (dt * step))
    for domain in active_species_pool.species_list():
        structure_output.write('%s    %g\n' % (domain, active_species_pool.get_population(domain)))
    log = open(clargs.sequence + '.log', 'w+')
    log.write('Step: %3d \n' % step)
    log.flush()
    def recombination(strand, current_length, all_foldons, all_domains, old_species_pool):
        result = []
        for terminal_foldon in all_foldons.find_foldons(strand.r_bound, current_length):
            # print(terminal_foldon)
            # print(terminal_foldon.get_IFR())
            result.append( (strand.elongate(terminal_foldon),
                                        old_species_pool.get_population(strand) /
                                        len(all_foldons.find_foldons(strand.r_bound, current_length))))
        for rearrange_point in reversed(strand.get_IFR()[:-1]):
            if rearrange_point == 0:  # Global rearrangement
                for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                    result.append( (overlapping_foldon, 0) )
            else:
                for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                    unrearranged_domain = all_domains.get_domain(strand.get_sequence()[0:rearrange_point],
                                                             strand.get_structure()[0:rearrange_point], 0,
                                                             rearrange_point)
                    result.append( (unrearranged_domain.elongate(overlapping_foldon), 0))
        return result

    def new_foldon_ss(sequence):  # foldon_collection is a DomainCollection
        mfe = nupack_functions.nupack_mfe(sequence, Temperature)
       # print(new_foldon.get_IFR())
        return mfe  # A very tricky solution

    while sequence_length > current_length:
        step += 1
        # print('Step: %3d \n'%step)
        log.write('Step: %3d \n'%step)
        log.write('Current length: %d \n'%current_length)
        old_species_pool = copy.deepcopy(active_species_pool)
        old_species_list = old_species_pool.species_list()
        active_species_pool.clear()
        log.flush()

        if current_length+dL > sequence_length:
            L_step = sequence_length-current_length
        else:
            L_step = dL
        current_length += L_step

        # Generate all new foldons
        l_bounds = np.arange(0, current_length, dL)
        multi_pool = Pool()
        log.write('Calculate new foldons...')
        log.flush()
        
        # new_foldon_sss = multi_pool.map(new_foldon_ss, [full_sequence[l_bounds[i]:current_length] for i in range(len(l_bounds))])
        # for i in range(len(l_bounds)):
        #     ss = new_foldon_ss(full_sequence[l_bounds[i]:current_length])
        #     sequence, l_bound, r_bound = full_sequence[l_bounds[i]:current_length], l_bounds[i], current_length
        #     new_foldon = all_domains.get_domain(sequence, ss, l_bound, r_bound)  # TODO: degeneracy
        #     new_foldon.foldonize()
        #     if new_foldon not in all_foldons.collection[new_foldon.l_bound, new_foldon.r_bound]:
        #         all_foldons.add_foldon(new_foldon)
 
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
        
        new_species = list(multi_pool.starmap(recombination, [(old_species_list[i], current_length, all_foldons,
                all_domains, old_species_pool) for i in range(len(old_species_list))]))  # parallelization
        print(new_species)
        multi_pool.close()
        multi_pool.join()
        for species_set in new_species:
            for species in species_set:
                active_species_pool.add_species(species[0],species[1])


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

    with gzip.open(clargs.sequence + '_domains.p.gz', 'w') as checkpoint_domains:
        pickle.dump(all_domains, checkpoint_domains)

    checkpoint_pool.close()
    log.close()
    structure_output.close()

exit()
