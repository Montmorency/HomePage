import os
import re
import sys

#script to convert LATTICE output from CamRecLib to extended xyz format
#File should have format
#TITLE( LATTICE COORDINATES)
#Index1  x y z Index2 x y z Index3 x y z Index4 x y z

#


exyz_str = 'cutoff=-1.00000000 nneightol=1.20000000 pbc="F F F" Lattice="5.000000       0.00000000       0.00000000       0.00000000 5.000000       0.00000000       0.00000000       0.00000000      6.000000" Properties=species:S:1:pos:R:3:Z:I:1'

with open('lattice.txt','r') as f:
  latt_file = f.read()

latt_coords = latt_file.split('\n')[1:-1]
new_latt = []
for latt in latt_coords:
  latt = latt.split()
  for n in range(len(latt)/4): 
    a,b,c,d = latt[4*n:4*n+4]
    new_latt.append([int(a), float(b), float(c), float(d)])

print len(new_latt)
print exyz_str
for latt_vec in new_latt: 
  frmt_str = '{At}{ls}{x}{ss}{y}{ss}{z}{ss}{Z}'.format(At='Fe',x=latt_vec[1],y=latt_vec[2],z=latt_vec[3],Z=26,ss=6*' ', ls=14*' ')
  print frmt_str
    


