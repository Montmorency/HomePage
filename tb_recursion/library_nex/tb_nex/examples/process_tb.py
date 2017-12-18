import re
import os
import sys
import numpy as np

def pull_rec_coeffs():
  with open('fort.1', 'r') as f:
    rec_coeffs = f.read()
  #split fort file of coefficients.
  rec_coeffs = rec_coeffs.split('CONTINUED FRACTION COEFFICIENTS A(I),B2(I),I=1, 11\n')
  #exclude first line.
  N_orbitals = len(rec_coeffs[1:])

  rec_mat = np.zeros([5, 11, 2])
  for num_orbital, line in enumerate(rec_coeffs[1:]):
  #each orbital has N pairs (a,b) recursion coefficients.
    for num_coeff, ab in enumerate(line.split('\n')[:-1]):
  #each row of recursion coefficients a_{i}, b_{i}.
      rec_mat[num_orbital, num_coeff, :] = map(float, ab.split())
  print rec_mat


if __name__=='__main__':
  pull_rec_coeffs()


