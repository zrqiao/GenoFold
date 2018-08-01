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
Temperature=37
##


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

class domain(object):
    def __init__(self, sequence, structure, lbound, rbound, collection, reducible=True, population=0):
        # [ , ) rbound[i]==lbound[i+1](stopping base +1)
        self.sequence = sequence
        self.structure = structure
        self.lbound = lbound
        self.rbound = rbound
        self.length = rbound-lbound
        self.population = population
        self.reducible = reducible
        self.G = None
        self.elements = set()
        self.IFR = [] #Irreducible foldons representation
        self.collection = collection #collection is a dictionary
        collection[(self.sequence, self.structure, self.lbound, self.rbound)] = self

    def __eq__(self, other):
        return self.structure == other.structure and self.lbound == other.lbound and self.rbound == other.rbound

    def __hash__(self):
        return hash((self.sequence, self.structure, self.lbound, self.rbound))

    def __str__(self):
        return self.structure

    def get_domain(self, sequence, structure, lbound, rbound):#Get another domain from collection or create a new one
        key = (sequence, structure, lbound, rbound)
        if key in self.collection.keys():
            return self.collection[key]
        else:
            return domain(sequence, structure, lbound, rbound, self.collection)

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

    def rate_nonoverlap(self, other):#Forward rate constant i->j for non-overlap domains
        #NOTE: This version is a very primary estimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        G_loop1=self.dissociate_to_loop().calc_G()
        G_loop2=other.dissociate_to_loop().calc_G()
        G_empty=.0
        dG_forward= G_loop1 - self.calc_G() + G_loop2 - G_empty
        dG_backward = G_loop2 - other.calc_G() + G_loop1 - G_empty
        return (np.exp(-(dG_forward)), np.exp(-(dG_backward)))

    def elongate(self, d_seq, d_ss, dl):
        return self.get_domain(self.sequence+d_seq, self.structure+d_ss, self.lbound, self.rbound+dl)

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

    def pathway_link(self, other):
        trial = self.check_availability(other)
        if trial:
            sub_lbound, sub_rbound = trial
            mysite = self.get_domain(self.sequence[sub_lbound:sub_rbound],
                                                       self.structure[sub_lbound:sub_rbound],
                                                       sub_lbound, sub_rbound)
            othersite = other.get_domain(other.sequence[sub_lbound:sub_rbound],
                                     other.structure[sub_lbound:sub_rbound],
                                     sub_lbound, sub_rbound)
            nocontributionelements = mysite.get_elements() ^ othersite.get_elements()
            mycontribution = mysite.get_elements() - nocontributionelements
            othercontribution = othersite.get_elements() - nocontributionelements
            #TODO: rate calculation interface
        #Add to collection

class foldon(domain):
    def __init__(self, adomain, foldon_collection):
        domain.__init__(self, adomain.sequence, adomain.structure, adomain.lbound, adomain.rbound, foldon_collection)

class pathway(object):
    def __init__(self, source, sink, collection, rate=0):
        self.source = source
        self.sink = sink
        self.rate = rate
        self.collection = collection
        collection[source][sink] = rate

pathway_collection = defaultdict(defaultdict(pathway))
