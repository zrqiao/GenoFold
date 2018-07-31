import numpy as np
import scipy as sp
import nupack_functions
import re
#Change following routines for other environments:

global Temperature
Temperature=37

class domain(object):
    def __init__(self, sequence, structure, lbound, rbound, reducible=True):# [ , ) rbound[i]==lbound[i+1](stopping base +1)
        self.sequence = sequence
        self.structure = structure
        self.lbound = lbound
        self.rbound = rbound
        self.length = rbound-lbound
        self.reducible = reducible
        self.G = None

    def decompose(self):#NOTE: No peudo knots
        unpaired = []
        elements = []
        sub_lbound = self.lbound
        if self.reducible:
            for index in range(self.length):
                symb = self.structure[index]
                base = self.sequence[index]
                if symb == '.':
                    if not unpaired:
                        elements.append(domain(base, symb, self.lbound+index, self.rbound+index+1))
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
                            elements.append(domain(self.sequence[sub_lbound:sub_rbound],
                                                   self.structure[sub_lbound:sub_rbound],
                                                   sub_lbound, sub_rbound, reducible=False))
                    except IndexError: #If ')' appears alone
                        print('Error: Invalid secondary structure')
                else:
                    print('Error: Invalid symbol')
            if unpaired:
                print('Error: Unfinished structure')
        return elements

    def calc_G(self):#NOTE: G is not initialized, has to be called explicitly before k calculation
        if not self.G:
            self.G = nupack_functions.nupack_ss_free_energy(sequence=self.sequence, ss=self.structure, T=Temperature)
        return self.G

    def dissociate_to_loop(self):
        #nonempty = re.search(r"^\s*(\S.*?)\s*$", self.sequence)
        #cent, stid, endid = nonempty.group(1), nonempty.start(1), nonempty.end(1)
        loopss=[]
        for sub_ss in self.decompose():#Can be a '.' or a helix
            if sub_ss.structure != '.':#Is a helix?
                #helix_core = sub_ss.structure[1:-1] #Check here!
                loopss.append('('+'.'*(sub_ss.length-2)+')')
            else: #Is a '.'?
                loopss.append(sub_ss.structure)
        loopstate = domain(self.sequence, ''.join(loopss), self.lbound, self.rbound)
        return loopstate

    def rate(self, other):#Forward rate constant i->j
        #NOTE: This version is a very primary esimation for unfolding-folding free energy that did not discriminate
        #  stem/ hairpin/ interior/ multi loops. This part is left for optimization.
        G_loop1=self.dissociate_to_loop().calc_G()
        G_loop2=other.dissociate_to_loop().calc_G()
        G_empty=.0
        dG_forward= G_loop1 - self.calc_G() + G_loop2 - G_empty
        dG_backward = G_loop2 - other.calc_G() + G_loop1 - G_empty
        return (np.exp(-(dG_forward)), np.exp(-(dG_backward)))


