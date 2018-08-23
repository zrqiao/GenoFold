import nupack_functions
from difflib import SequenceMatcher
import numpy as np
from collections import defaultdict
import operator
from decimal import *
import mpmath as mp
from scipy.linalg import expm
import scipy.linalg as linalg
import copy
import time
from sklearn import preprocessing
from sympy import exp, N, S
from sympy.matrices import Matrix
# import matlab
# import matlab.engine

#Change following routines for other environments:
Temperature = 37
k0 = 1
R = 1.9858775e-3  # G in kcal/mol
rate_cutoff = 1e-20  # minimum allowed rate constant
subopt_gap=0.99
mp.mp.prec = 333
mp.mp.dps = 60
mp.mp.pretty = True
##


def eigen(M):
    # print(M)
    # M = mp.matrix(M.tolist())
    # eigenValues, eigenVectors = mp.eig(M)
    # eigenValues = np.diag(np.array(eigenMatrix).astype('float128'))
    E, EL, ER = mp.eig(M, left = True, right = True)
    E, EL, ER = mp.eig_sort(E, EL, ER)
    # print(ER*mp.diag(E)*EL)
    eigenVectors = ER
    eigenValues = E
    # eigenVectors = np.array(ER.apply(mp.re).tolist(), dtype=float)
    # eigenValues = np.array([mp.re(x) for x in E], dtype=float)
    # idx = eigenValues.argsort()[::-1]
    # eigenValues = eigenValues
    if len(eigenVectors.shape) == 1:
        eigenVectors = [eigenVectors]
    # print(eigenValues)
    # print(eigenVectors)
    return eigenValues, eigenVectors


'''
def Propagate(M, p, time):

    e, U = eigen(M)
    # print(e)
    # the eigenvalues are distinct -- possibly complex, but
    # E will always be real
    Uinv = linalg.inv(U)
    # print(np.exp(time*np.diag(e)))
    E = np.real(np.dot(np.dot(U, np.diag(np.exp(time*e))), Uinv))
    # print(E)
    p1 = np.dot(p, E)
    return p1
'''


def Propagate_stationary(M, p, dt, ddt=1):

    E, EL, ER = mp.eig(M.T, left = True, right = True)
    E, EL, ER = mp.eig_sort(E, EL, ER)
    times = range(ddt, dt+ddt, ddt)
    UR = ER**-1
    R = [0]*(len(E)-1)
    R.append(1)
    # print(M.T)
    # print(E)
    intermediate_populations = []
    # print(np.exp(time*np.diag(e)))
    # E = np.real(np.dot(np.dot(U, np.diag(R)), Uinv))
    for i in range(len(times)):
        # print(mp.diag(R[i]))
        A = ER*mp.diag(R)*UR*p
        intermediate_populations.append(np.array((A.T).apply(mp.re).tolist()[0], dtype=float))
    print(intermediate_populations)
    # intermediate_populations = [np.array(((ER*mp.diag(R)*EL*p).T).apply(mp.re).tolist()[0], dtype=float) for t in times]
    # print(intermediate_populations)
    return intermediate_populations


def Propagate_trunc2(M, p, dt, ddt=1):

    E, EL, ER = mp.eig(M.T, left = True, right = True)
    E, EL, ER = mp.eig_sort(E, EL, ER)
    times = np.arange(0, dt, ddt)+dt
    N = int(dt/ddt)
    if len(p) == 1:
        intermediate_populations = [p for i in range(N)]
        return intermediate_populations
    else:
        R = [[0 for i in range(len(p))] for j in range(N)]
        for i in range(N) : 
            R[i][-2] = mp.exp(E[-2]*times[i])
            R[i][-1] = 1
        # print(M)
        # print(ER)
        print(E)
        # print(R)
        # print(EL)
        # print(p)
        # print(ER*mp.diag(R[0])*EL*p)
        intermediate_populations = []
        for i in range(N):
            # print(mp.diag(R[i]))
            A = ER*mp.diag(R[i])*EL*p
            intermediate_populations.append(np.array((A.T).apply(mp.re).tolist()[0], dtype=float))
        print(intermediate_populations)
        return intermediate_populations


def Propagate(M, p, dt, ddt=1): 
    E, EL, ER = mp.eig(M.T, left = True, right = True)
    UR = ER**-1
    # E, EL, ER = mp.eig_sort(E, EL, ER)   # time_series = np.arange(0, dt, ddt) + dt
    # print(mp.nstr(EL*ER, n=3))
    times = range(ddt, dt+ddt, ddt)
    # print(E)
    intermediate_populations = []
    # print(np.exp(time*np.diag(e)))
    # E = np.real(np.dot(np.dot(U, np.diag(R)), Uinv))
    for i in range(len(times)):
        R = [mp.exp(E[j]*times[i]) for j in range(len(E))]
        # R = [mp.exp(E[j]*0) for j in range(len(E))]
        A = ER*mp.diag(R)*UR*p
        # print(R)
        intermediate_populations.append(np.array((A.T).apply(mp.re).tolist()[0], dtype=float))
    # print(p)
    # print(intermediate_populations[-1])
    # intermediate_populations = [np.array(((ER*mp.diag(R)*EL*p).T).apply(mp.re).tolist()[0], dtype=float) for t in times]
    # print(intermediate_populations)
    return intermediate_populations



def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def rate(dG, k):
    rate = k*mp.exp(-(dG/(R*(273.15+Temperature))))
    # print(str(rate))
    return rate
    # return k


def boltzmann_factor(dG):
    return rate(dG, 1)


def disso(x): return x.dissociate_energy()


def loopf(x): return x.loop_formation_energy()


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
        if key in self.collection.keys():
            return self.collection[key]
        else:
            new_domain = Domain(sequence, structure, l_bound, r_bound, self)
            self.add_domain(new_domain)  # Add to repository when created
            return new_domain

    def add_domain(self, adomain):
        self.collection[(adomain.sequence, adomain.structure, adomain.l_bound, adomain.r_bound)] = adomain  # Add to repository when created


class FoldonCollection(object):
    def __init__(self):
        self.collection = defaultdict(set)

    def add_foldon(self, adomain):
        if adomain.is_foldon:
            self.collection[adomain.l_bound, adomain.r_bound].add(adomain)  # Add to repository when created
            # TODO: degenerate case
            return self
        else:
            return False

    def find_foldons(self, l_bound, r_bound):
        # Get possible foldon configurations from collection, cannot create new instances
        try:
            return self.collection[l_bound, r_bound]
        except IndexError:
            print("Error: no such foldon")
            return False

    def new_foldon(self, sequence, l_bound, r_bound, domain_collection, ss=None):  # foldon_collection is a DomainCollection
        if not ss:
            sss = nupack_functions.nupack_subopt(sequence, Temperature, subopt_gap)  # NOTE: subopt
            for ss in sss:
                self.new_foldon(sequence, l_bound, r_bound, domain_collection, ss=ss)
            return True
        newfoldon = domain_collection.get_domain(sequence, ss, l_bound, r_bound)
        newfoldon.foldonize()
        if newfoldon not in self.collection[l_bound, r_bound]:
            self.add_foldon(newfoldon)  # Don't forget to modify
        # print(new_foldon.get_IFR())
        return True


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
        self.IFR = np.array([])  # Irreducible foldons representation
        #  NOTE: IFR is automatically generated when initiate foldons
        self.collection = collection  # NOTE: collection is a domains_collection
        self.loop_state = None

    def __eq__(self, other):
        return self.structure == other.structure and self.l_bound == other.l_bound and self.r_bound == other.r_bound

    def __hash__(self):
        return hash((self.sequence, self.structure, self.l_bound, self.r_bound))

    def __str__(self):
        return self.structure

    def __repr__(self):
        return '^'.join(map(str, self.IFR))  # all information is encoded by IFR

    def get_domain(self, sequence, structure, l_bound, r_bound):#Get another domain from collection or create a new one
        # key = (sequence, structure, lbound, rbound)
        # if key in self.collection.keys():
        #    return self.collection[key]
        # else:
        return self.collection.get_domain(sequence, structure, l_bound, r_bound)

    def get_sequence(self):
        return self.sequence

    def get_structure(self):
        return self.structure

    def get_IFR(self):
        return self.IFR

    def foldonize(self):  # If self is a foldon
        self.is_foldon = True
        self.IFR = np.array([self.l_bound, self.r_bound])  # overwrite
        return True

    def get_elements(self):  # find closed domains. NOTE: No peudo knots
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
                            new_element = self.get_domain(base, symb, self.l_bound + index, self.l_bound + index + 1)
                            new_element.reducible = False
                            new_element.elements = set([new_element])
                            elements.append(new_element)
                        else:
                            pass
                    elif symb == '(':
                        if not unpaired:
                            sub_l_bound= self.l_bound + index
                        unpaired.append(symb)
                    elif symb == ')':
                        try:
                            unpaired.pop()
                            if not unpaired:
                                sub_r_bound = self.l_bound + index + 1
                                new_element = self.get_domain(self.sequence[sub_l_bound-self.l_bound:sub_r_bound-self.l_bound],
                                                       self.structure[sub_l_bound-self.l_bound:sub_r_bound-self.l_bound],
                                                              sub_l_bound, sub_r_bound)
                                new_element.reducible = False
                                new_element.elements = set([new_element])
                                elements.append(new_element)
                        except IndexError:  # If ')' appears alone
                            print('Error: Invalid secondary structure')
                    else:
                        print('Error: Invalid symbol')
                if unpaired:
                    print('Error: Unfinished structure')
            if len(elements) == 1:
                self.reducible = False
            # print([e.structure for e in elements])
            self.elements = set(elements)
        return self.elements

    def get_G(self):  # NOTE: G is not initialized, has to be called explicitly before k calculation
        if self.G is None:  # NOTE:0.00=False!
            if self.l_bound !=0:  # Redirect to domain with the same ss
                self.G = self.get_domain(self.sequence, self.structure, 0, self.length).get_G()
            else:
                # time1=time.time()
                self.G = nupack_functions.nupack_ss_free_energy(sequence=self.sequence, ss=self.structure, T=Temperature)
                time2=time.time()
                # print(time2-time1)
        # print(self.G)
        return self.G

    def dissociate_to_loop(self):#Input is an element
        # nonempty = re.search(r"^\s*(\S.*?)\s*$", self.sequence)
        # cent, stid, endid = nonempty.group(1), nonempty.start(1), nonempty.end(1)
        # Can be a '.' or a helix
        # time_1 = time.time()
        if self.loop_state: return self.loop_state
        if self.reducible:
            print('need decompose')
            return False
        else:
            if self.structure != '.':  # Is a helix?
                # helix_core = sub_ss.structure[1:-1] #Check here!
                loop_ss = '('+'.'*(self.length-2)+')'
            else:  # Is a '.'?
                loop_ss = self.structure
        loop_state = self.get_domain(self.sequence, ''.join(loop_ss), self.l_bound, self.r_bound)
        # time_2 = time.time()
        # print(time_2-time_1)
        self.loop_state = loop_state
        # print(loop_state.sequence)
        return loop_state

    def dissociate_energy(self):
        # time_1=time.time()
        # NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        # if self.structure == '.':
        #    return 0
        # else:
        G_loop = self.dissociate_to_loop().get_G()
        # G_empty=.0
        # dG_forward= G_loop1 - self.calc_G() + G_loop2 - G_empty
        # dG_backward = G_loop2 - other.calc_G() + G_loop1 - G_empty
        # print(G_loop)
        # time_2=time.time()
        # print(time_2-time_1)
        return G_loop-self.get_G()

    def loop_formation_energy(self):
        # NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        # if self.structure == '.':
        #    return 0
        G_empty = .0
        # else:
        G_loop = self.dissociate_to_loop().get_G()
        return G_loop-G_empty

    def elongate(self, additional_domain):
        # Primary structure have a IFR; elongation structure is a Irreducible foldon;
        # generate another valid domain from one
        if self.r_bound != additional_domain.l_bound:
            print('Error: illegal elongation')
            return False
        else:
            longer_domain = self.get_domain(self.sequence + additional_domain.sequence,
                                            self.structure + additional_domain.structure,
                                            self.l_bound, additional_domain.r_bound)
            # update IFR
            if longer_domain.IFR.size:
                if longer_domain.IFR[-2] <= self.r_bound:  # Is a better IFR
                    return longer_domain
            else:
                longer_domain.IFR = np.append(self.IFR, additional_domain.r_bound)  # Direct deepcopy can be very slow
            # print('segment after elongation: ')
            # print(longer_domain.get_IFR())
            return longer_domain

    def check_availability(self, other): #Forward rearrangement
        IFRa = set(self.IFR)
        IFRb = set(other.IFR)
        
        if IFRa >= IFRb:
            diff = IFRa - IFRb
            if not diff:  # NOTE: cannot discern identical domains here
                return other.IFR[0], other.IFR[-1]  # Global rearrangement
            for i in range(len(other.IFR)-1):
                if max(diff) < other.IFR[i+1] and min(diff) > other.IFR[i]:
                    return other.IFR[i], other.IFR[i+1]  # Exact indices of rearrangement site
            return False
        else:
            return False
        # Add to collection

# def pathway_rate(self, other, pathway_collection):#Input


# class foldon(domain):
#    def __init__(self, adomain, foldon_collection):
#        domain.__init__(self, adomain.sequence, adomain.structure, adomain.lbound, adomain.rbound, foldon_collection)
#        self.isfoldon = True


class Pathways(object):  # Indices of a pathway should be two domain(for robustness, can be improved by using IFR for indices)
    def __init__(self, k):  # Need pre-exponential factor
        self.collection = defaultdict(dict)
        self.k = k

    def has_path(self, source, sink):
        if sink in self.collection[source]:
            return True
        else:
            return False

    def get_rate(self, source, sink):
        if not self.has_path(source, sink):
            self.pathway_link(source, sink, self.k)
        return self.collection[source][sink]

    def add(self, source, sink, rate):
        self.collection[source][sink] = rate

    def pathway_link(self, domain1, domain2, k):    # If no path exists, call this
        # NOTE: finding the minimum rearrangement site
        trial_forward = domain1.check_availability(domain2)
        if trial_forward:
            trial = trial_forward
            source, sink = domain1, domain2
        else:
            trial_backward = domain2.check_availability(domain1)
            if trial_backward:
                trial = trial_backward
                source, sink = domain2, domain1
            else:
                forward_rate = backward_rate = 0
                self.add(domain1, domain2, forward_rate)
                self.add(domain2, domain1, backward_rate)
                return forward_rate, backward_rate
        sub_l_bound, sub_r_bound = trial
        if sub_l_bound == domain1.l_bound and sub_r_bound == domain2.r_bound:  # This is the minimum rearrangement site
            my_site, other_site= source, sink
        else:
            my_site = source.get_domain(source.sequence[sub_l_bound-source.l_bound:sub_r_bound-source.l_bound],
                                            source.structure[sub_l_bound-source.l_bound:sub_r_bound-source.l_bound],
                                                       sub_l_bound, sub_r_bound)
            other_site = sink.get_domain(sink.sequence[sub_l_bound-sink.l_bound:sub_r_bound-sink.l_bound],
                                            sink.structure[sub_l_bound-sink.l_bound:sub_r_bound-sink.l_bound],
                                     sub_l_bound, sub_r_bound)
        # NOTE: rate calculation interface
        # if already exists
        if self.has_path(my_site, other_site):  # If calculated
            forward_rate = self.get_rate(my_site, other_site)
            backward_rate = self.get_rate(other_site, my_site)
        else:
            # Rate calculation
            # time1 = time.time()
            no_contribution_elements = my_site.get_elements() & other_site.get_elements()
            my_contributions = my_site.get_elements() - no_contribution_elements
            # print('mysite: ')
            # print(sorted([[e.structure, e.l_bound, e.r_bound] for e in my_contributions],key=operator.itemgetter(2)))
             
            their_contributions = other_site.get_elements() - no_contribution_elements
            # print('othersite: ')
            # print(sorted([[e.structure, e.l_bound, e.r_bound] for e in their_contributions],key=operator.itemgetter(2)))
            # time1 = time.time()
            # print(time2-time1)
 
            # print([[e.structure, e.l_bound, e.r_bound] for e in other_site.get_elements()])

            # NOTE: This is the actual computational bottleneck 
            # multi_pool = Pool()
            # print(len(my_contributions))
            # for x in my_contributions:
                # print(x, x.sequence)
                # multi_pool.apply_async(disso, x)
            forward_energy = sum(map(disso, my_contributions)) + \
                            sum(map(loopf, their_contributions))
            # time2 = time.time()
            # print(time2-time1)

            backward_energy = sum(map(disso, their_contributions)) + \
                             sum(map(loopf, my_contributions))
            # multi_pool.close()
            # multi_pool.join()
            # time2 = time.time()

            forward_rate = rate(forward_energy, k)
            backward_rate = rate(backward_energy, k)

            # print(time2-time1)
            # print(backward_energy)
            # newforwardpathway = pathway(my_site,other_site, forward_rate)
            # newbackwardpathway = pathway(other_site,my_site, backward_rate)
            self.add(my_site, other_site, forward_rate)
            self.add(other_site, my_site, backward_rate)
            time2 = time.time()
            # print(time2-time1)
        self.add(source, sink, forward_rate)
        self.add(sink, source, backward_rate)
        return forward_rate, backward_rate


class SpeciesPool(object):
    def __init__(self, pathways):
        self.species = defaultdict(float)
        self.pathways = pathways
        self.size = 0
        self.timestamp = 0

    def add_species(self, domain, population=0.):
        if domain.get_IFR().size:
            if domain not in self.species : self.size += 1
            self.species[domain] += population  # NOTE: duplication means more

    def clear(self):
        self.species = defaultdict(float)
        self.size = 0
        return self

    def __deepcopy__(self, memo):
        memo[id(self)] = newself = self.__class__(self.pathways)  # NOTE: pathway collection is shallow copied!
        newself.size = copy.deepcopy(self.size, memo)
        newself.timestamp = copy.deepcopy(self.timestamp, memo)
        newself.pathways = self.pathways
        for myspecies in self.species:
            newself.add_species(myspecies, population=self.get_population(myspecies))
        return newself

    def species_list(self):
        return list(self.species.keys())

    def get_population(self, domain):
        return self.species[domain]

    def update_population(self, domain, population):
        self.species[domain] = population
        return self

    def evolution(self, pathways, dt, ddt, stationary=False):
        # print(self.size)
        # eng = matlab.engine.start_matlab()
        rate_matrix = mp.matrix(self.size, self.size)
        species_list = list(self.species.items())
        population_array = mp.matrix([0 for i in range(self.size)])
        # TODO: parallel optimization

        self.timestamp += dt
        time_array = np.arange(0, dt, ddt) + self.timestamp + ddt
        for i in range(self.size):
            diag_i = 0
            population_array[i] = species_list[i][1]
            for j in range(self.size):
                if j == i: continue
                rate = pathways.get_rate(species_list[i][0], species_list[j][0])
                diag_i += rate
                rate_matrix[i, j] = rate
            rate_matrix[i, i] = -diag_i

        # rate_matrix = matlab.double(rate_matrix)
        # population_array = matlab.double(population_array)

        if stationary:
            '''
            intermediate_population_arrays = \
                preprocessing.normalize([[rate(species[0].get_G(), 1) for species in species_list]
                                        for t in time_array], norm='l1', axis=1)
            '''
            intermediate_population_arrays = Propagate_stationary(rate_matrix, population_array, dt, ddt=ddt)
        else:
            '''
            k_fastest = np.max(rate_matrix)
            for i in range(self.size):  # Make it a REAL sparse matrix
                for j in range(self.size):
                    if rate_matrix[i][j] < k_fastest * rate_cutoff:
                        rate_matrix[i][j] = 0
                rate_matrix[i][i] = -np.sum(rate_matrix[i])
            '''
            # Master Equation
            intermediate_population_arrays = Propagate(rate_matrix, population_array, dt, ddt=ddt)
            # TODO: This is only for perturbation!
        population_array = intermediate_population_arrays[-1]
        # Remapping
        for i in range(self.size):
            self.update_population(species_list[i][0], population_array[i])

        return species_list, intermediate_population_arrays, time_array

    def selection(self, size_limit):
        if self.size > size_limit:
            ordered_species = list(sorted(self.species.items(), key=operator.itemgetter(1), reverse=True))[:size_limit]
            # print(ordered_species)
            remaining_population = sum([a[1] for a in ordered_species])
            self.clear()
            for species in ordered_species:
                self.add_species(species[0], population=species[1]/remaining_population)
        return self


def recombination(strand, current_length, all_foldons, all_domains, old_species_pool, active_species_pool):
    # test_time_1 = time.time()
    for rearrange_point in strand.get_IFR():
        if rearrange_point == 0:  # Global rearrangement
            for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                active_species_pool.add_species(overlapping_foldon)
        elif rearrange_point == strand.r_bound:
            terminal_set = all_foldons.find_foldons(strand.r_bound, current_length)
            terminal_pfunc = np.sum([boltzmann_factor(fd.get_G()) for fd in terminal_set])
            for terminal_foldon in terminal_set:
                active_species_pool.add_species(strand.elongate(terminal_foldon),
                                                population=old_species_pool.get_population(strand) *
                                                boltzmann_factor(terminal_foldon.get_G()) / terminal_pfunc)
        else:
            for overlapping_foldon in all_foldons.find_foldons(rearrange_point, current_length):
                unrearranged_domain = all_domains.get_domain(strand.get_sequence()[0:rearrange_point],
                                                             strand.get_structure()[0:rearrange_point], 0,
                                                             rearrange_point)
                active_species_pool.add_species(unrearranged_domain.elongate(overlapping_foldon))
    # test_time_2 = time.time()
    # print(test_time_2-test_time_1)
