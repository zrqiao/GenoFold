import numpy as np
import scipy as sp
import nupack_functions
from difflib import SequenceMatcher
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import re
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm, expm_multiply

#Change following routines for other environments:
global Temperature
global k, k0
Temperature = 37
k = 1
k0 = 1

foldons = defaultdict(dict)
domains = dict()
##


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def rate(dG):
    return k0*np.exp(-(k*dG))

class domains_collection():

    def __init__(self):
        self.collection = defaultdict(dict)

    def get_domain(self, sequence, structure, lbound, rbound):#Get another domain from collection or create a new one
        key = (sequence, structure, lbound, rbound)
        if key in self.collection.keys():
            return self.collection[key]
        else:
            return domain(sequence, structure, lbound, rbound, self.collection)

    def add_domain(self, adomain):
        self.collection[(adomain.sequence, adomain.structure, adomain.lbound, adomain.rbound)] = adomain  # Add to repository when created


def get_foldon(sequence, lbound, rbound, domain_collection, foldon_collection):#foldon_collection[lbound][rbound]=object(foldon)
    mfe = nupack_functions.nupack_mfe(sequence, Temperature)
    newfoldon = domain_collection.get_domain(sequence, mfe, lbound, rbound)
    newfoldon.foldonize()
    if newfoldon not in foldon_collection[lbound].values():
        foldon_collection[lbound][rbound] = newfoldon
    return newfoldon

class domain(object):
    def __init__(self, sequence, structure, lbound, rbound, collection, population=0):
        # [ , ) rbound[i]==lbound[i+1](stopping base +1)
        self.sequence = sequence
        self.structure = structure
        self.lbound = lbound
        self.rbound = rbound
        self.length = rbound-lbound
        self.reducible = True
        self.G = None
        self.elements = set()
        self.isfoldon = False
        self.population = population
        self.IFR = None #Irreducible foldons representation NOTE: IFR is automatically generated when initiate foldons
        self.collection = collection #NOTE: collection is a domains_collection
        self.collection.add_domain(self) #Add to repository when created

    def __eq__(self, other):
        return self.structure == other.structure and self.lbound == other.lbound and self.rbound == other.rbound

    def __hash__(self):
        return hash((self.sequence, self.structure, self.lbound, self.rbound))

    def __str__(self):
        return self.structure

    def __repr__(self):
        return self.lbound, self.rbound, self.structure

    def get_domain(self, sequence, structure, lbound, rbound):#Get another domain from collection or create a new one
        #key = (sequence, structure, lbound, rbound)
        #if key in self.collection.keys():
        #    return self.collection[key]
        #else:
        return self.collection.get_domain(sequence, structure, lbound, rbound)

    def get_structure(self):
        return self.structure

    def foldonize(self):#If self is a foldon
        self.isfoldon = True
        self.IFR = [self.lbound, self.rbound]
        return True

    def change_population(self,population):
        self.population = population

    def get_elements(self): #find closed domains. NOTE: No peudo knots
        if not self.elements:
            unpaired = []
            elements = []
            sub_lbound = self.lbound
            if self.reducible:
                for index in range(self.length):
                    symb = self.structure[index]
                    base = self.sequence[index]
                    if symb == '.':
                        if not unpaired:
                            elements.append(self.get_domain(base, symb, self.lbound+index, self.rbound+index+1))
                        else:
                            pass
                    elif symb == '(':
                        sub_lbound=self.lbound+index
                        unpaired.append(symb)
                    elif symb == ')':
                        try:
                            unpaired.pop()
                            if not unpaired:
                                sub_rbound = self.lbound+index+1
                                newelement = self.get_domain(self.sequence[sub_lbound:sub_rbound],
                                                       self.structure[sub_lbound:sub_rbound],
                                                       sub_lbound, sub_rbound)
                                newelement.reducible = False
                                elements.append(newelement)
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

    def calc_G(self):#NOTE: G is not initialized, has to be called explicitly before k calculation
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
                loopss = '('+'.'*(self.length-2)+')'
            else: #Is a '.'?
                loopss = self.structure
        loopstate = self.get_domain(self.sequence, ''.join(loopss), self.lbound, self.rbound)
        return loopstate

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
        return G_loop-self.calc_G()

    def loop_formation_energy(self):
        #NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        #if self.structure == '.':
        #    return 0
        G_empty = .0
        #else:
        G_loop = self.dissociate_to_loop().calc_G()
        return G_loop-G_empty

    def elongate(self, d_seq, d_ss, dl):#Primary structure have a IFR; elongation structure is a Irreducible Foldon
        longerdomain= self.get_domain(self.sequence+d_seq, self.structure+d_ss, self.lbound, self.rbound+dl)
        if not longerdomain.IFR: #update IFR
            longerdomain.IFR = self.IFR
            longerdomain.IFR.append(longerdomain.rbound)
        return longerdomain

    def check_availability(self, other): #Forward rearrangement
        IFRa = set(self.IFR)
        IFRb = set(other.IFR)
        if IFRa & IFRb == IFRb:
            diff = IFRa ^ IFRb
            for i in range(len(other.IFR)-1):
                if max(diff) < other.IFR[i+1] and min(diff) > other.IFR[i]:
                    return other.IFR[i], other.IFR[i+1] #Exact indices of rearrangement site
            return False
        else:
            return False

    def pathway_link(self, other, pathways):
        trial = self.check_availability(other)
        if trial:
            sub_lbound, sub_rbound = trial
            mysite = self.get_domain(self.sequence[sub_lbound:sub_rbound],
                                                       self.structure[sub_lbound:sub_rbound],
                                                       sub_lbound, sub_rbound)
            othersite = other.get_domain(other.sequence[sub_lbound:sub_rbound],
                                     other.structure[sub_lbound:sub_rbound],
                                     sub_lbound, sub_rbound)
            forwardrate = 0; backwardrate = 0
            #NOTE: rate calculation interface
            #if already exists
            forwardrate = pathways.get_rate(mysite, othersite)  # TODO: implementation
            backwardrate = pathways.get_rate(othersite, mysite)
            if not forwardrate and backwardrate:
                nocontributionelements = mysite.get_elements() ^ othersite.get_elements()
                mycontributions = mysite.get_elements() - nocontributionelements
                theircontributions = othersite.get_elements() - nocontributionelements
                forwardenergy = sum(map(lambda x:x.dissociate_energy(), mycontributions)) +\
                                sum(map(lambda x:x.loop_formation_energy(), theircontributions))
                backwardenergy = sum(map(lambda x: x.dissociate_energy(), theircontributions)) + \
                                sum(map(lambda x: x.loop_formation_energy(), mycontributions))
                forwardrate = rate(forwardenergy)
                backwardrate = rate(backwardenergy)
                #newforwardpathway = pathway(mysite,othersite, forwardrate)
                #newbackwardpathway = pathway(othersite,mysite, backwardrate)
                pathways.add(mysite, othersite, forwardrate)
                pathways.add(othersite, mysite, backwardrate)
            pathways.add(self, other, forwardrate)
            pathways.add(other, self, backwardrate)

            return forwardrate, backwardrate
        #Add to collection

    ###def pathway_rate(self, other, pathway_collection):#Input

#class foldon(domain):
#    def __init__(self, adomain, foldon_collection):
#        domain.__init__(self, adomain.sequence, adomain.structure, adomain.lbound, adomain.rbound, foldon_collection)
#        self.isfoldon = True

class pathways(object):#Indices of a pathway should be two domain(for robustness, can be improved by using IFR for indices)
    def __init__(self):
        self.collection = defaultdict(defaultdict(np.float))

    def get_rate(self, source, sink):
        return self.collection[source][sink]

    def add(self, source, sink, rate):
        self.collection[source][sink] = rate

pathway_collection = defaultdict(dict)

class species(object):
    def __init__(self, pathways):
        self.domains = set()
        self.pathways = pathways
