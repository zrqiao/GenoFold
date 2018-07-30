import numpy as np
import scipy as sp
import nupack_functions
import re
#Change following routines for other environments:

global Temperature
Temperature=37

class domain(object):
    def __init__(self, sequence, structure, lbound, rbound,reducible=True):# [ , ) rbound[i]==lbound[i+1](stopping base +1)
        self.sequence = sequence
        self.structure = structure
        self.lbound = lbound
        self.rbound = rbound
        self.length = rbound-lbound
        self.reducible = reducible
        self.G = None

    def decompose(self):
        unpaired = []
        elements = []
        sub_lbound = self.lbound
        if self.reducible:
            for index in range(self.length):
                symb = self.structure[index]
                base = self.sequence[index]
                if symb == '.':
                    if not unpaired:
                        elements.append(domain(base,symb,self.lbound+index,self.rbound+index+1))
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
                                                   ''.join(unpaired), sub_lbound, sub_rbound, reducible=False))
                    except:
                        print('Error: Invalid secondary structure')
                else: print('Error: Invalid symbol')
            if unpaired:
                print('Error: Unfinished structure')
        return elements

    def calc_G(self):#NOTE: G is not initialized, has to be called explicitly before k calculation
        self.G=nupack_functions.nupack_ss_free_energy(sequence=self.sequence, ss=self.structure, T=Temperature)
        return self.G

    def dissociate_to_loop(self):
        #nonempty = re.search(r"^\s*(\S.*?)\s*$", self.sequence)
        #cent, stid, endid = nonempty.group(1), nonempty.start(1), nonempty.end(1)
        if not self.sequence.strip:
            loopss = self.sequence
        else:
            cent, stid, endid=self.sequence.strip('.'), self.length-len(self.sequence.lstrip('.')),\
                    len(self.sequence.rstrip('.'))
            loopss = '.'*stid+'('+'.'*(endid-stid-2)+')'+'.'*endid
        loopstate = domain(self.sequence, '(' + loopss, self.lbound, self.rbound)
        return loopstate

    def rate(self,other):#Forward rate constant i->j
        G_loop1=self.dissociate_to_loop().calc_G()
        G_loop2=other.dissociate_to_loop().calc_G()
        G_empty=.0
        dG_forward= G_loop1 - self.G + G_loop2 - G_empty
        dG_backward = G_loop2 - other.G + G_loop1 - G_empty
        return (np.exp(-(dG_forward)), np.exp(-(dG_backward)))


