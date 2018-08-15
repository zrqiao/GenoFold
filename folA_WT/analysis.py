from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

plt.style.use('ggplot')
fig = plt.figure(figsize=(10, 10))

colors = [plt.cm.jet(lt) for lt in range(0,8)]
fig.add_axes()

#mpl.rcParams['axes.color_cycle'] = colors
mpl.rcParams['axes.titlesize']=20
mpl.rcParams['axes.titleweight'] = 15
TPTDistribution=np.array([])


fig = plt.figure(figsize=(10, 7.5))
fig.add_axes()
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 2.5])

ax_energy = fig.add_subplot(gs[0,0])
ax_energy.set_title('Free Energy')
#ax_energy.set_xlabel('Subsequence Length', fontsize=12.5)
ax_energy.set_ylabel('Free Energy', fontsize=12.5)

ax_structure = fig.add_subplot(gs[1,0])
ax_structure.set_title('First N nt local structure predominance among subsequences')
ax_structure.set_xlabel('Subsequence Length', fontsize=12.5)
ax_structure.set_ylabel('Similarity', fontsize=12.5)

with open('summary_energy.dat','r') as energy_dat:
    energy=np.loadtxt(energy_dat)
    seqcount=np.arange(1,len(energy)+1)
    ax_energy.plot(seqcount,energy,ls='solid')

with open('summary_structure.dat','r') as structure_dat:
    f=open('Similarity.dat','w')
    structure=np.array(list(map(lambda x: x.rstrip('\n'),structure_dat.readlines())))
    seqcount=np.arange(1,len(energy)+1)
    for N in (np.arange(5,50,5)):
        print('N=' + str(N) + ' nt')
        f.write('#N=' + str(N) + ' nt')
        sim=np.zeros(len(structure))
        first_sequences=np.array(list(map(lambda x:x[:N],structure[N:])))
        for i in range(len(first_sequences)):
            temp_cmp=[]
            for cmpseq in first_sequences:
                temp_cmp.append(similar(cmpseq,first_sequences[i]))
            sim[i+N]=np.average(np.array(temp_cmp))
            #print(sim)
        ax_structure.plot(seqcount[N:],sim[N:],ls='solid',label='N='+str(N)+' nt')
        np.savetxt(f,sim)
    ax_structure.legend(loc='best')
fig.tight_layout()
plt.show()

fig.savefig('G&Similarity.eps')