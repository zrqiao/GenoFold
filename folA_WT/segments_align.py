import numpy as np

with open("../folA_WT.in", 'r') as m:
    sequence=m.readline()

print('RNA='+sequence)

energy=open('summary_energy.dat','w')
structure=open('summary_structure.dat','w')
summary_pairs=open('summary_pairs.dat','w')

for i in range(1,len(sequence)+1):
    with open('segments/sample_'+str(i)+'.in', 'w') as w:
        w.write(sequence[:i]+'\n')
    with open('segments/sample_'+str(i)+'.mfe', 'r') as result:
        q=result.readlines()
        # w.write(q[15])
        energy.write(q[14])
        structure.write(q[15])

    with open('segments/sample_'+str(i)+'.ppairs', 'r+') as result_pairs:
        # summary_pairs.write(f'# {i}')
        raw_data = [list(map(float, line.split())) for line in result_pairs if not line.startswith('%')]
        data = np.zeros(i)
        p_unbound = np.zeros(i)+1
        for dat in raw_data:
            if len(dat) == 3:
                if dat[1] > i:  # p_unbound
                    data[int(dat[0]-1)] = dat[2]
        summary_pairs.write(' '.join(map(str, data))+'\n')
energy.close()
structure.close()
summary_pairs.close()
