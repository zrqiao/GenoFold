import shlex, subprocess
from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import Domains
import nupack_functions

#Change following routines for other environments:
L_init = 20 #Initiation unit
dL = 5 #elongation unit (also CG unit)
dt = 1 #Folding time for each elongation step

#TODO: Initiation [create active population]

#TODO: structure_generation(single strain, elgation segment) [to be called in pool.map()]
#Compute all IFR segments; link sequences; update IFRs

#TODO: population_selection (need a active population, fittness function)

#TODO: compute_foldon(i,j)

#TODO: need a collection for all foldons

#TODO: A active map for all pathways