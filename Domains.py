import nupack_functions
from difflib import SequenceMatcher
import numpy as np
from collections import defaultdict
import operator
from multiprocessing import Pool
import re
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm, expm_multiply

#Change following routines for other environments:
Temperature = 37
k = 1
k0 = 1

##


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def rate(dG):
    return k0*np.exp(-(k*dG))


class DomainsCollection(object):
    def __init__(self):
        self.collection = dict()

    def __contains__(self, adomain):
        if (adomain.sequence, adomain.structure, adomain.l_bound, adomain.r_bound) in self.collection:
            return True
        else:
            return False

    def get_domain(self, sequence, structure, l_bound, r_bound):#Get another domain from collection or create a new one
        key = (sequence, structure, l_bound, r_bound)
        if key in self.collection:
            return self.collection[key]
        else:
            new_domain = Domain(sequence, structure, l_bound, r_bound, self)
            self.add_domain(new_domain)  # Add to repository when created
            return new_domain

    def add_domain(self, adomain):
        self.collection[(adomain.sequence, adomain.structure, adomain.l_bound, adomain.r_bound)] = adomain  # Add to repository when created


class FoldonCollection(object):
    def __init__(self):
        self.collection = defaultdict(lambda: defaultdict(set))

    def add_foldon(self, adomain):
        if adomain.is_foldon:
            self.collection[adomain.l_bound][adomain.r_bound].add(adomain)  # Add to repository when created
            # TODO: degenerate case
            return self
        else:
            return False

    def find_foldons(self, l_bound, r_bound):  # Get possible foldon configurations from collection, cannot create new instances
        try:
            return self.collection[l_bound][r_bound]
        except IndexError:
            print("Error: no such foldon")
            return False

    def new_foldon(self, sequence, l_bound, r_bound, domain_collection):#foldon_collection is a DomainCollection
        mfe = nupack_functions.nupack_mfe(sequence, Temperature)
        new_foldon = domain_collection.get_domain(sequence, mfe, l_bound, r_bound)  #	TODO: degeneracy
        new_foldon.foldonize()
        if new_foldon not in self.collection[l_bound][r_bound]:
            self.add_foldon(new_foldon)
        print(new_foldon.get_IFR())
        return new_foldon

class Domain(object):
    def __init__(self, sequence, structure, l_bound, r_bound, collection):
        # [ , ) rbound[i]==lbound[i+1](stopping base +1)
        self.sequence = sequence
        self.structure = structure
        self.l_bound = l_bound
        self.r_bound = r_bound
        self.length = r_bound - l_bound
        self.reducible = True
        self.G = None
        self.elements = set()
        self.is_foldon = False
        self.IFR = None #Irreducible foldons representation NOTE: IFR is automatically generated when initiate foldons
        self.collection = collection #NOTE: collection is a domains_collection

    def __eq__(self, other):
        return self.structure == other.structure and self.l_bound == other.lbound and self.r_bound == other.rbound

    def __hash__(self):
        return hash((self.sequence, self.structure, self.l_bound, self.r_bound))

    def __str__(self):
        return self.structure

    def __repr__(self):
        return self.structure

    def get_domain(self, sequence, structure, lbound, rbound):#Get another domain from collection or create a new one
        #key = (sequence, structure, lbound, rbound)
        #if key in self.collection.keys():
        #    return self.collection[key]
        #else:
        return self.collection.get_domain(sequence, structure, lbound, rbound)

    def get_sequence(self):
        return self.sequence

    def get_structure(self):
        return self.structure

    def get_IFR(self):
        return self.IFR

    def foldonize(self):#If self is a foldon
        self.is_foldon = True
        self.IFR = [self.l_bound, self.r_bound]
        return True

    def get_elements(self): #find closed domains. NOTE: No peudo knots
        if not self.elements:
            unpaired = []
            elements = []
            sub_l_bound = self.l_bound
            if self.reducible:
                for index in range(self.length):
                    symb = self.structure[index]
                    base = self.sequence[index]
                    if symb == '.':
                        if not unpaired:
                            elements.append(self.get_domain(base, symb, self.l_bound + index, self.r_bound + index + 1))
                        else:
                            pass
                    elif symb == '(':
                        sub_l_bound= self.l_bound + index
                        unpaired.append(symb)
                    elif symb == ')':
                        try:
                            unpaired.pop()
                            if not unpaired:
                                sub_r_bound = self.l_bound + index + 1
                                new_element = self.get_domain(self.sequence[sub_l_bound:sub_r_bound],
                                                       self.structure[sub_l_bound:sub_r_bound],
                                                              sub_l_bound, sub_r_bound)
                                new_element.reducible = False
                                elements.append(new_element)
                        except IndexError: #If ')' appears alone
                            print('Error: Invalid secondary structure')
                    else:
                        print('Error: Invalid symbol')
                if unpaired:
                    print('Error: Unfinished structure')
            if len(elements)==1:
                self.reducible = False
            self.elements = set(elements)
        return self.elements

    def get_G(self):#NOTE: G is not initialized, has to be called explicitly before k calculation
        if not self.G:
            self.G = nupack_functions.nupack_ss_free_energy(sequence=self.sequence, ss=self.structure, T=Temperature)
        return self.G

    def dissociate_to_loop(self):#Input is an element
        #nonempty = re.search(r"^\s*(\S.*?)\s*$", self.sequence)
        #cent, stid, endid = nonempty.group(1), nonempty.start(1), nonempty.end(1)
        #Can be a '.' or a helix
        if self.reducible:
            print('need decompose')
            return False
        else:
            if self.structure != '.':#Is a helix?
                #helix_core = sub_ss.structure[1:-1] #Check here!
                loop_ss = '('+'.'*(self.length-2)+')'
            else: #Is a '.'?
                loop_ss = self.structure
        loop_state = self.get_domain(self.sequence, ''.join(loop_ss), self.l_bound, self.r_bound)
        return loop_state

    def dissociate_energy(self):
        #NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        #if self.structure == '.':
        #    return 0
        #else:
        G_loop = self.dissociate_to_loop().calc_G()
        #G_empty=.0
        #dG_forward= G_loop1 - self.calc_G() + G_loop2 - G_empty
        #dG_backward = G_loop2 - other.calc_G() + G_loop1 - G_empty
        return G_loop-self.get_G()

    def loop_formation_energy(self):
        #NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        #if self.structure == '.':
        #    return 0
        G_empty = .0
        #else:
        G_loop = self.dissociate_to_loop().calc_G()
        return G_loop-G_empty

    def elongate(self, additional_domain):#Primary structure have a IFR; elongation structure is a Irreducible foldon;
        # generate another valid domain from one
        if self.r_bound != additional_domain.r_bound:
            print('Error: illegal elongation')
            return False
        else:
            longer_domain = self.get_domain(self.sequence + additional_domain.sequence,
                                            self.structure + additional_domain.structure,
                                            self.l_bound, additional_domain.r_bound)
            if not longer_domain.IFR: #update IFR
                longer_domain.IFR = self.IFR
                longer_domain.IFR.append(longer_domain.rbound)
            return longer_domain

    def check_availability(self, other): #Forward rearrangement
        IFRa = set(self.IFR)
        IFRb = set(other.IFR)
        if IFRa == IFRb:
            return False
        elif IFRa & IFRb == IFRb:
            diff = IFRa ^ IFRb
            for i in range(len(other.IFR)-1):
                if max(diff) < other.IFR[i+1] and min(diff) > other.IFR[i]:
                    return other.IFR[i], other.IFR[i+1] #Exact indices of rearrangement site
            return False
        else:
            return False
        #Add to collection

    ###def pathway_rate(self, other, pathway_collection):#Input

#class foldon(domain):
#    def __init__(self, adomain, foldon_collection):
#        domain.__init__(self, adomain.sequence, adomain.structure, adomain.lbound, adomain.rbound, foldon_collection)
#        self.isfoldon = True


class Pathways(object):     #  Indices of a pathway should be two domain(for robustness, can be improved by using IFR for indices)
    def __init__(self):
        self.collection = defaultdict(dict)

    def has_path(self, source, sink):
        if sink in self.collection[source]:
            return True
        else:
            return False

    def get_rate(self, source, sink):
        if not self.has_path(source, sink):
            self.pathway_link(source, sink)
        return self.collection[source][sink]

    def add(self, source, sink, rate):
        self.collection[source][sink] = rate

    def pathway_link(self, domain1, domain2):    #If no path exists, call this
        # NOTE: finding the minimum rearrangement site
        trial_forward = domain1.check_availability(domain1)
        if trial_forward:
            trial = trial_forward
            source, sink = domain1, domain2
        else:
            trial_backward = domain2.check_availability(domain2)
            if trial_backward:
                trial = trial_backward
                source, sink = domain2, domain1
            else:
                forward_rate = backward_rate = 0
                self.add(domain1, domain2, forward_rate)
                self.add(domain2, domain1, backward_rate)
                return forward_rate, backward_rate
        sub_l_bound, sub_r_bound = trial
        if sub_l_bound == domain1.lbound and sub_r_bound == domain2.rbound:#This is the minimum rearrangement site
            my_site, other_site= source, sink
        else:
            my_site = source.get_domain(source.sequence[sub_l_bound:sub_r_bound],
                                                       source.structure[sub_l_bound:sub_r_bound],
                                                       sub_l_bound, sub_r_bound)
            other_site = sink.get_domain(sink.sequence[sub_l_bound:sub_r_bound],
                                     sink.structure[sub_l_bound:sub_r_bound],
                                     sub_l_bound, sub_r_bound)
        # NOTE: rate calculation interface
        # if already exists
        if self.has_path(my_site, other_site):# If calculated
            forward_rate = self.get_rate(my_site, other_site)
            backward_rate = self.get_rate(other_site, my_site)
        else:
            #Rate calculation
            no_contribution_elements = my_site.get_elements() ^ other_site.get_elements()
            my_contributions = my_site.get_elements() - no_contribution_elements
            their_contributions = other_site.get_elements() - no_contribution_elements
            forward_energy = sum(map(lambda x: x.dissociate_energy(), my_contributions)) + \
                            sum(map(lambda x: x.loop_formation_energy(), their_contributions))
            backward_energy = sum(map(lambda x: x.dissociate_energy(), their_contributions)) + \
                             sum(map(lambda x: x.loop_formation_energy(), my_contributions))
            forward_rate = rate(forward_energy)
            backward_rate = rate(backward_energy)
            # newforwardpathway = pathway(my_site,other_site, forward_rate)
            # newbackwardpathway = pathway(other_site,my_site, backward_rate)
            self.add(my_site, other_site, forward_rate)
            self.add(other_site, my_site, backward_rate)
        self.add(source, sink, forward_rate)
        self.add(sink, source, backward_rate)
        return forward_rate, backward_rate


class SpeciesPool(object):
    def __init__(self, pathways):
        self.species = defaultdict(float)
        self.pathways = pathways
        self.size = 0
        self.timestamp = 0

    def add_species(self, domain, population=0):
        if domain.IFR:
            self.species[domain] += population  # NOTE: duplication means more
            self.size += 1

    def clear(self):
        self.species.clear()
        self.size = 0
        return self

    def species_list(self):
        return list(self.species.keys())

    def get_population(self, domain):
        return self.species[domain]

    def update_population(self, domain, population):
        self.species[domain] = population
        return self

    def evolution(self, pathways, time):
        print(self.size)
        rate_matrix = np.zeros((self.size, self.size))
        species_list = list(self.species.items())
        population_array = np.zeros(self.size)
        # TODO: parallel optimization
        for i in range(self.size):
            population_array[i] = species_list[i][1]
            for j in range(self.size):
                rate_matrix[i][j] = pathways.get_rate(species_list[i][0], species_list[j][0])
            rate_matrix[i][i] = -np.sum(rate_matrix[i])
	
        print(list(population_array))
        # Master Equation
        population_array = population_array.dot(expm(time*rate_matrix))
        self.timestamp += time

        print(list(population_array))

        # Remapping
        for i in range(self.size):
            self.update_population(species_list[i][0], population_array[i])

        return self

    def selection(self, size_limit):
        if self.size > size_limit:
            ordered_species = list(sorted(self.species.items(), key=operator.itemgetter(1), reverse=True))[:size_limit]
            remaining_population = sum(list(zip(ordered_species))[1])
            self.clear()
            for species in ordered_species:
                self.add_species(species[0], population=species[1]/remaining_population)
        return self


def recombination(strand, current_length, all_foldons, all_domains, old_species_pool, active_species_pool):
    for terminal_foldon in all_foldons.find_foldons(strand.r_bound, current_length):
        active_species_pool.add_species(strand.elongate(terminal_foldon),
                                        population=old_species_pool.get_population(strand) /
                                        len(all_foldons.find_foldons(strand.r_bound, current_length)))
    for rearrange_point in reversed(strand.get_IFR()[:-1]):
        if rearrange_point == 0:  # Global rearrangement
            active_species_pool.add_species(all_foldons.find_foldons(rearrange_point, current_length))
        else:
            for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                unrearranged_domain = all_domains.get_domain(strand.get_sequence()[0:rearrange_point],
                                                             strand.get_structure()[0:rearrange_point], 0,
                                                             rearrange_point)
                active_species_pool.add_species(unrearranged_domain.elongate(overlapping_foldon))
