import os, sys


def rna_seq(sequence):
    # Convert a DNA sequence to a RNA sequence
    return ''.join(c if c != 'T' else 'U' for c in sequence)

UTR = 'ATACCCGTTTTTTGGGCTAACAGGAGGAATTACAT'
input = sys.argv[1]

with open(input, 'r+') as f:
    for line in f.readlines():
        if line.startswith('>'):
            name=sys.argv[2]+'_'+line[1:].strip('\n')
            if not os.path.exists(name):
                os.makedirs(name)
            with open('mutants.in', 'a') as f1:
                f1.write(name+'\n')
        elif line.rstrip():
            dnaseq=UTR+line.strip('\n')
            rnaseq=rna_seq(dnaseq)
            with open('sequences/'+name+'.in', 'w+') as rnaf:
                rnaf.write(rnaseq)
    
        

       
