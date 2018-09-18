import os, sys
mut=sys.argv[1]
with open('sequences/'+mut+'.in','r') as m:
    sequence=m.readline()

# print('RNA='+sequence)

#energy=open(mut+'summary_energy.dat','w')
#structure=open(mut+'summary_structure.dat','w')


for i in range(1,len(sequence)+1):
    with open(mut+'/segments/sample_'+str(i)+'.in','w') as w:
        w.write(sequence[:i]+'\n')
        # with open(mut+'/segments/sample_'+str(i)+'.mfe','r') as result:
            # q=result.readlines()
            # w.write(q[15])
            # energy.write(q[14])
            # structure.write(q[15])

#energy.close()
#structure.close()
