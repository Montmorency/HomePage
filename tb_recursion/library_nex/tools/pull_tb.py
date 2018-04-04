import os
import sys
import pickle
import numpy as np

with open('tb.dat','r') as f:
  tb_file = f.read()

#for line in tb_file.split('\n'):
#  print line

tb_lines = tb_file.split('\n')

date_created = tb_lines[0]

cell_vec_1 = tb_lines[1]
cell_vec_2 = tb_lines[2]
cell_vec_3 = tb_lines[3]

cell_latt = np.zeros([3,3])
cell_vec_1 = np.array(map(float, cell_vec_1.split()))
cell_vec_2 = np.array(map(float, cell_vec_2.split()))
cell_vec_3 = np.array(map(float, cell_vec_3.split()))

cell_latt = np.array([cell_vec_1, cell_vec_2, cell_vec_3])

print 'cell_latt', cell_latt

matel = tb_lines[5:-1]
print matel[0], matel[1]

print len(matel), 
wan_blocks = len(matel)/326

atoms_dict = {}

for i in range(wan_blocks+1):
  start = i*326
  end = (i+1)*326
  chunk = matel[start:end]
  position = tuple(map(int, chunk[0].split()))
  print position
  H_wan = np.zeros([18, 18], dtype=complex)
  wan_ham_list = [x.split() for x in chunk[1:]]
  for line in wan_ham_list[:-1]:
    i,j = int(line[0])-1, int(line[1])-1
    H_wan[i][j] = complex(float(line[2]), float(line[3]))

  atoms_dict[position] = H_wan

with open('wan_tb.pckl','w') as f:
  pickle.dump(atoms_dict,f)


