import sys

import pickle
import numpy  as np
import tensorflow as tf

from ase import units

def direction_cosine(vec):
  """
  The vectors are represented as 
  p*a*\hat{i}, q*a*j, r*a*k where p,q,r
  are all integers.

  l = p^{2}/(p^{2}+q^{2}+r^{2})^{1/2}
  m = q^{2}/(p^{2}+q^{2}+r^{2})^{1/2}
  n = r^{2}/(p^{2}+q^{2}+r^{2})^{1/2}
  """
  l = vec[0]/np.linalg.norm(vec)
  m = vec[1]/np.linalg.norm(vec)
  n = vec[2]/np.linalg.norm(vec)
  return np.array([l,m,n])

def fit_tb(features, labels, verbose=True):
  """
  The input data consists of the direction cosines connecting the two orbitals
  and the indices of the orbitals. The appropriate coefficients of the energy 
  integrals in the two center approximation are considered the features to be
  passed in. The weight they are accorded will correspond to the energy
  integrals (ss\sigma, dd\delta, etc.). The optimal weights as determined by
  gradient descent will constitute the optimal tight binding model.

  Args:
    n_matel(int): Number of matrix elements to fit

  E_{s,s}(l,m,n) = (ss\sigma)
  E_{xy,xy}(l,m,n) = (3l^{2}m^{2})*(dd\sigma)+ (l^{2}+m^{2}-4.0*l^{2}*m^{2})*(dd\pi) + (n^{2}+l^{2}*m^{2})*(dd\delta)
  """
#https://stackoverflow.com/questions/35298326/freeze-some-variables-scopes-in-tensorflow-stop-gradient-vs-passing-variables
# If we are interested in freezing the sub space on the ddsigma, dddelta, 
# and ddpi parameters we can do that.
# sp3d2:
# SP3D2 = [ss\sigma. sp\sigma, sd\sigma, pp\sigma, pp\pi, pd\pi, dd\sigma, dd\pi, \dd\delta]
# DD =[dd\sigma, dd\pi, dd\delta]
# Initialize with pettifor parameters
#  dd_sig = tf.Variable(-0.065, trainable=True, name='dd_sig')
#  dd_pi = tf.Variable(0.043, trainable=True, name='dd_pi')
#  dd_delta = tf.Variable(-0.010, trainable=True, name='dd_delta')

  dd_sig = tf.Variable(0.0, trainable=True, name='dd_sig')
  dd_pi = tf.Variable(0.0, trainable=True, name='dd_pi')
  dd_delta = tf.Variable(0.0, trainable=True, name='dd_delta')

  DD = tf.Variable([dd_sig, dd_pi, dd_delta])
  R  = tf.placeholder(shape=[len(features), 3], dtype=tf.float32, name='direction_cosines')
  #Matrix element Wannier
  M_wan  = tf.placeholder(shape=[len(labels)], dtype=tf.float32, name='matrix_elements')
  #Matrix element tb
  M_tb = tf.reduce_sum(tf.multiply(DD,R),axis=1)

  #define lost, optimizer, coeffs to minimize
  loss = tf.reduce_sum((M_wan-M_tb)**2)
  #loss = tf.reduce_sum(tf.abs(M_wan-M_tb))

  optimizer = tf.train.AdamOptimizer(0.0001)
  train_step = optimizer.minimize(loss, var_list=[DD])
  nsteps = 10000
  init_op = tf.initialize_all_variables()
  #np.random.shuffle(features)
  feed_dict={R:features, M_wan:labels}
  print '\n'
  print '\n'
  print 'Fitting TB Coefficients'
  with tf.Session() as sess:
    sess.run(init_op)
    for n in range(nsteps):
      sess.run(train_step, feed_dict={R:features, M_wan:labels})
      if (n+1) % 500 == 0:
        #print('iter %i, %f' % (n+1, objval))
        print 'DD-params:', sess.run(DD,feed_dict={R:features, M_wan:labels}), 'loss: ', sess.run(loss, feed_dict=feed_dict)
        dd_out =  sess.run(DD,feed_dict={R:features, M_wan:labels})
        loss_out = sess.run(loss, feed_dict=feed_dict)

# To freeze a subspace.
# train_W = optimizer.minimize(loss, var_list=[W])
# train_C = optimizer.minimize(loss, var_list=[C])
  print '\n'
  print 'TB coefficients', dd_out, 'loss:', loss_out
  print '\n'
  dd_out = np.array(dd_out)

  if verbose:
    print 'Direction Cosine', 'DD Parameters', 'TB Matrix element', 'Wannier Matrix Element'
    for feat, label in zip(features, labels):
      print (' '.join(['{:3.3f}'.format(round(x,4)) for x in feat]), ' '.join(['{:3.3f}'.format(round(x,4)) for x in dd_out]), 
             'SK:', np.round(dd_out.dot(feat),5), 'target:', np.round(label,5))
  return

def matel_coeffs(l,m,n,wan_i,wan_j):
  """
  Takes direction cosines and returns
  the appropriate combination of direction cosines for
  the required energy integral.
  Args:
    l,m,n: direction cosines (p/(p^{2}+q^{2}+r^{2})).
  The orbital indices are wan_i, wan_j  (Energy Integral)
  xy, yz, zx 
  """
  #E_xy,xy
  l2 = np.square(l)
  m2 = np.square(m)
  n2 = np.square(n)
  #Diagonal interactions (dd\sigma, dd\pi, dd\delta)
  if (wan_i == 'xy' and wan_j == 'xy'):
    return [3.0*l2*m2, l2+m2-4.0*l2*m2, n2+l2*m2]
  elif (wan_i == 'yz' and wan_j == 'yz'):
    return [3.0*m2*n2, m2+n2-4.0*m2*n2, l2+m2*n2]
  elif (wan_i == 'zx' and wan_j == 'zx'):
    return [3.0*n2*l2, n2+l2-4.0*n2*l2, m2+n2*l2]
  elif (wan_i == 'xy' and wan_j == 'zx'):
    return [-3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'xy' and wan_j == 'yz'): 
    return [3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'xy'): 
    return [-3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'yz'):
    return [3.0*m*n2*l, m*l*(1.0-4.0*n2), m*l*(n2-1.0)]
  elif (wan_i == 'yz' and wan_j == 'xy'):
    return [3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'yz' and wan_j == 'zx'): 
    return [3.0*m*n2*l, m*l*(1.0-4.0*n2), m*l*(n2-1.0)]
  #cyclic permutation vs. SK paper. check works!
  #elif wan_i == 'xy' and wan_j='zx':
  #  return [3.0*l2*m*n, m*n*(1.0-4.0*l2), m*n*l2-1.0)]
  else:
    sys.exit('Invalid Combo')

def load_matels(l,m,n,subspace_matrix):
  """
  Takes direction cosines, and 3x3 subspace matrix.
  """
  matel_dict={}
  #diagonal
  matel_dict[('xy', 'xy')] = (0, 0)
  matel_dict[('zx', 'zx')] = (1, 1)
  matel_dict[('yz', 'yz')] = (2, 2)
  #off diagonal
  matel_dict[('xy', 'zx')] = (0, 1)
  matel_dict[('xy', 'yz')] = (0, 2)
  matel_dict[('zx', 'xy')] = (1, 0)
  matel_dict[('zx', 'yz')] = (1, 2)
  matel_dict[('yz', 'xy')] = (2, 0)
  matel_dict[('yz', 'zx')] = (2, 1)
  pairs = [('xy', 'xy'), ('yz', 'yz'), ('zx', 'zx'), ('xy', 'zx'), ('xy', 'yz'),
           ('zx', 'xy'), ('zx', 'yz'), ('yz', 'xy'), ('yz', 'zx')]
  features = []
  matels = []
  for pair in pairs:
      i,j = matel_dict[pair]
      feature = matel_coeffs(l,m,n, pair[0], pair[1])
      matel = subspace_matrix[i][j]
      features.append(feature)
      matels.append(matel)
  return features, matels

def gen_input(Fe_nn, R_ji, spin_mask=None):
  """
  Args:
    atom_dict(dict): Dictionary keyed off by atomic coordinates, containing Hamiltonian
    matrix elements between Wannier functions on different centers.
    nn_list: List of coordinates (as tuples) for the desired neighbour shell.
    spin_mask: If present only spin up or spin down (alternating rows) are considered.
  """
  features, labels = [], []
  #for nn in nn_list:
  #Fe_nn = atom_dict[nn]/units.Ry
  if spin_mask != None:
    Fe_nn = Fe_nn[spin_mask,:]
    Fe_nn = Fe_nn[:, spin_mask]
  subspace_matrix = Fe_nn[-3:,-3:].real
  l, m, n = direction_cosine(np.array(R_ji))
  features_tmp, labels_tmp = load_matels(l, m, n, subspace_matrix)
  features.extend(features_tmp)
  labels.extend(labels_tmp)
  return features, labels
  
def pprint(H, mask=None, units=None):
  """
  Pretty print matrix.
  """
  rows, cols = np.shape(H)
  if units != None:
    H = 1./units * H
  for i in range(rows):
    if mask == None:
      frmt_str = [' {:3.4f}'.format(x) for x in H[i,:].real.round(4)]
    else:
      frmt_str = [' {:3.4f}'.format(x) for x in H[i,mask].real.round(4)]
    pretty_string = [] 
    for x in frmt_str:
      if '-' in x and abs(float(x)) >= 10.:
        pass
      elif '-' in x:
        x = ' '+x
      else:
        x = '  '+x
      pretty_string.append(x)
    print ' '.join(pretty_string)

with open('wan_tb.pckl') as f:
  atom_dict = pickle.load(f)

with open('tb.dat','r') as f:
  tb_file = f.read()

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
alat = cell_latt[0][0]

#On site
Fe_00 = atom_dict[(0,0,0)]
spinup_mask = [i%2!=0 for i in range(18)]
spindown_mask = [i%2==0 for i in range(18)]
print 'Full On Site Matrix:'
pprint(Fe_00, units=units.Ry)
#print 'Spin Up On Site Matrix (with spin/spin):'
#pprint(Fe_00[spinup_mask,:])
#print 'Spin Down On Site Matrix (with spin/spin):'
#pprint(Fe_00[spindown_mask,:])
print 'Spin Up On Site Matrix (without spin/spin):'
pprint(Fe_00[spinup_mask,:], mask=spinup_mask, units=units.Ry)
#print 'Spin Down On Site Matrix (without spin/spin):'
#pprint(Fe_00[spindown_mask,:], mask=spindown_mask)

NN_1s = [(-1,0,0),(1,0,0),(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(1,-1,1),(-1,1,-1)] 
pairs = [('xy', 'xy'), ('yz', 'yz'), ('zx', 'zx'), ('xy', 'zx'), ('xy', 'yz'),
         ('zx', 'xy'), ('zx', 'yz'), ('yz', 'xy'), ('yz', 'zx')]

combined_feats = []
combined_labels = []
## Nearest Neighbours
for NN in NN_1s:
  print 'Spin Up 1NN Matrix (without spin/spin):'
  print NN
  cart = cell_latt.dot(np.array(NN))
  cart_int = cart/alat
  print 'Cartesian:', cart
  print 'Cartesian (Int):', cart_int
  print 'Direction Cosine:', direction_cosine(cart_int)
  print '\n'
  Fe_1NN = atom_dict[NN]
  pprint(Fe_1NN[spinup_mask,:], mask=spinup_mask, units=units.Ry)
  print '\n'
  spinup_subspace = Fe_1NN[spinup_mask,:]
  spinup_subspace = spinup_subspace[:, spinup_mask]
  pprint(spinup_subspace[-3:,-3:],units=units.Ry)
  features, labels =  gen_input(Fe_1NN/units.Ry, cart_int, spin_mask=spinup_mask)
  combined_feats.extend(features)
  combined_labels.extend(labels)

fit_tb(combined_feats, combined_labels)

#2nd NearestNeighbours
combined_feats = []
combined_labels = []
NN_2s = [(-1,1,0), (1,-1,0), (0,-1,1), (0,1,-1), (-1,0,-1), (1,0,1)] 
for NN in NN_2s:
  print 'Spin Up 2NN Matrix (without spin/spin):'
  print NN
  #cart = cell_latt.T.dot(np.array(NN))
  cart = cell_latt.dot(np.array(NN))
  cart_int = cart/alat
  print 'Cartesian:', cart
  print 'Cartesian (Int):', cart_int
  print 'Direction Cosine:', direction_cosine(cart_int)
  print '\n'
  Fe_2NN = atom_dict[NN]
  pprint(Fe_2NN[spinup_mask,:], mask=spinup_mask, units=units.Ry)
  features, labels =  gen_input(Fe_2NN/units.Ry, cart_int, spin_mask=spinup_mask)
  combined_feats.extend(features)
  combined_labels.extend(labels)

fit_tb(combined_feats, combined_labels)

