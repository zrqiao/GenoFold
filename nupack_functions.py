import subprocess, random, os, sys

# Modify the following lines according to your NUPACK installation:
nupack_path = os.environ['HOME'] + '/.local/source/nupack3.0.6/bin'
nupack_env = {'NUPACKHOME' : os.environ['HOME'] + '/.local/source/nupack3.0.6/'}

def rna_seq(sequence):
    # Convert a DNA sequence to a RNA sequence
    return ''.join(c if c != 'T' else 'U' for c in sequence)

def nupack_free_energy(sequence, T):
    # Use NUPACK to calculate the -log partition function of a sequence
    # NOTE: This free energy sums over all possible secondary structures
    seq = rna_seq(sequence)
    with subprocess.Popen([nupack_path + '/pfunc', '-T', str(T)],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                          stderr=subprocess.PIPE, env=nupack_env) as proc:
        nupack_input = bytes("%s\n" % seq, 'ascii')
        nupack_output, nupack_err = proc.communicate(nupack_input)
        dG = float(nupack_output.decode('ascii').split('\n')[-3])
    return dG

def nupack_mfe(sequence, T):
    # Use NUPACK to calculate the minimum-free-energy secondary structure of a sequence
    # NOTE: Returns a secondary structure in the (((.))) notation
    seq = rna_seq(sequence)
    rint = int(random.random() * 1.e9)
    tmp = '/tmp/%d' % rint
    with open(tmp + '.in', 'w') as f:
        f.write("%s\n" % seq)
    subprocess.call([nupack_path + '/mfe', '-T', str(T), tmp], env=nupack_env)
    with open(tmp + '.mfe', 'r') as f:
        flag = False
        for line in f:
            if len(line) > 1 and all(c in '(.)' for c in line.strip()):
                ss = line.strip()
    os.remove(tmp + '.in')
    os.remove(tmp + '.mfe')
    return ss

def nupack_ss_free_energy(sequence, ss, T):
    # Use NUPACK to calculate the free energy of a specific secondary structure
    # NOTE: 'ss' should use the (((.))) notation
    seq = rna_seq(sequence)
    with subprocess.Popen([nupack_path + '/energy', '-T', str(T)],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                          stderr=subprocess.PIPE, env=nupack_env) as proc:
        nupack_input = bytes("%s\n%s\n" % (seq, ss), 'ascii')
        nupack_output, nupack_err = proc.communicate(nupack_input)
        dG = float(nupack_output.decode('ascii').split('\n')[-2])
    return dG
