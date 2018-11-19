import numpy as np

def g_const(x,a,b):
  if x.real >= 0:
    return ((x-a)/(2.0*b**2) - (0.25*(((x-a)/(b**2))**2) - 1.0/(b**2))**0.5)
  else:
    return ((x-a)/(2.0*b**2) + (0.25*(((x-a)/(b**2))**2) - 1.0/(b**2))**0.5)

def g_const_chain(x=0, a=0.0, b=0.4, emin=-3.0, emax=3.01, ediff=0.01):
  en_range = np.arange(emin, emax,0.01)
  g = np.array([0.0]*len(en_range), dtype=complex)
  for i, E in enumerate(en_range):
    x = E + 0.01j
    if x.real >= 0:
      g[i] = ((x-a)/(2.0*b**2) - (0.25*(((x-a)/(b**2))**2) - 1.0/(b**2))**0.5)
    else:
      g[i] = ((x-a)/(2.0*b**2) + (0.25*(((x-a)/(b**2))**2) - 1.0/(b**2))**0.5)
  return g, en_range

def g_defect(x=0, a_0=-1.2, b_0=0.4, a=0.0, b=0.4, emin=-3.0, emax=3.01, ediff=0.01):
  en_range = np.arange(emin, emax,0.01)
  g = np.array([0.0]*len(en_range), dtype=complex)
  for i, E in enumerate(en_range):
    x = E + 0.01j
    g[i] = (1.0/b_0)/((x-a_0)/b_0 - g_const(x,a,b))
  return g, en_range

g, en_range = g_const_chain()
with open('ge.dat','w') as f:
  for x, y in zip(en_range, g):
    print >>f, x, y.real, y.imag

a_0s =[-1.2, -1.2, -1.2]
b_0s =[0.04, 0.4, 1.2]
for inum, a_0, b_0 in zip(range(len(a_0s)), a_0s, b_0s):
  g, en_range = g_defect(a_0=a_0, b_0=b_0)
  with open('g_defect_{}.dat'.format(inum+1), 'w') as f:
    for x, y in zip(en_range, g):
      print >>f, x, y.real, y.imag

a_0s =[-0.1, -0.1, -0.1]
b_0s =[0.04, 0.4, 1.0]
for inum, a_0, b_0 in zip(range(len(a_0s)), a_0s, b_0s):
  g, en_range = g_defect(a_0=a_0, b_0=b_0)
  with open('g_inband_{}.dat'.format(inum+1), 'w') as f:
    for x, y in zip(en_range, g):
      print >>f, x, y.real, y.imag

a_0s = [12, 12, 12]
b_0s = [2., 4., 8.]
for inum, a_0, b_0 in zip(range(len(a_0s)), a_0s, b_0s):
  g, en_range = g_defect(a_0=a_0, b_0=b_0, a=0, b=0.1, emin=-5,emax=15, ediff=0.01)
  with open('g_plasmon_{}.dat'.format(inum+1), 'w') as f:
    for x, y in zip(en_range, g):
      print >>f, x, y.real, y.imag

gtmp = g.copy()
N = len(g)
gtmp[0] = g[N//2]
print en_range[N//2]
gtmp[1:N//2+1]=g[N//2+1:]
gtmp[N//2+1:]=g[:N//2]
print gtmp[0]
gt = np.fft.ifft(gtmp)
with open('gt.dat','w') as f:
  for  y in gt:
    print >>f, y.real, y.imag, y.real**2+y.imag**2

with open('dgt2.dat','w') as f:
  gt2 = [y.real**2+y.imag**2 for y in gt]
  for  y in np.diff(gt2):
    print >>f, y

