import sys
import numpy as np

path = sys.argv[1]
lmodel = sys.argv[2]
step = sys.argv[3]

output_string_t = "output/"+str(step)+"_"+str(lmodel)
output_string   = "output/"+str(lmodel)

# Posterior distributions
ola = np.loadtxt(path+output_string_t+'_postdist.txt')
part_a = ola[:,2:4]
part_b = ola[:,4:9]
part_c = ola[:,9:]
new = np.column_stack((ola[:,:2],part_b,part_a,part_c))
np.savetxt('tmp_postdist.txt',new) 


# Parameter names
old_lines = [line.rstrip('\n') for line in open(path+output_string+'_postdist.paramnames')]
new_lines = ['a' for i in range(0,len(old_lines)-1)]
new_lines[0:5] = old_lines[2:7]
new_lines[5:7] = old_lines[0:2]
new_lines[7:9] = old_lines[7:9]
new_lines[9:11] = old_lines[9:11]

with open('tmp_postdist.paramnames','w') as f:
    for line in new_lines:
        f.write(line+"\n")


# Parameter ranges
old_lines = [line.rstrip('\n') for line in open(path+output_string+'_postdist.ranges')]
new_lines = ['a' for i in range(0,len(old_lines)-1)]
new_lines[0:5] = old_lines[2:7]
new_lines[5:7] = old_lines[0:2]
new_lines[7:9] = old_lines[7:9]
new_lines[9:11] = old_lines[9:11]

with open('tmp_postdist.ranges','w') as f:
    for line in new_lines:
        f.write(line+"\n")
