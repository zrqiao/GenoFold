with open("RNA.in",'r') as m:
    sequence=m.readline()

print('RNA='+sequence)

energy=open('summary_energy.dat','w')
structure=open('summary_structure.dat','w')


for i in range(1,len(sequence)+1):
    with open('segments/sample_'+str(i)+'.in','w') as w:
        w.write(sequence[:i]+'\n')
    with open('segments/sample_'+str(i)+'.mfe','r') as result:
        q=result.readlines()
        energy.write(q[14])
        structure.write(q[15])

energy.close()
structure.close()
