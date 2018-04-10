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
#########################################################
#  Universal Ordering of Slater-Koster parameters:      #
#       [ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma, #
#        pd\sigma, pd\pi,                               #
#        dd\sigma,                                      #
#        dd\pi,                                         #
#        dd\delta]                                      #
#########################################################
# DD =[dd\sigma, dd\pi, dd\delta]
# Initialize with pettifor parameters

  ss_sig = tf.Variable(0.0, trainable=True, name='ss_sig')
  sp_sig = tf.Variable(0.0, trainable=True, name='sp_sig')
  pp_sig = tf.Variable(0.0, trainable=True, name='pp_sig')
  pp_pi = tf.Variable(0.0, trainable=True, name='pp_pi')
  sd_sig = tf.Variable(0.0, trainable=True, name='sd_sig')
  pd_sig = tf.Variable(0.0, trainable=True, name='pd_sig')
  pd_pi = tf.Variable(0.0, trainable=True, name='pd_pi')
  #dd_sig = tf.Variable(0.0, trainable=True, name='dd_sig')
  #dd_pi = tf.Variable(0.0, trainable=True, name='dd_pi')
  #dd_delta = tf.Variable(0.0, trainable=True, name='dd_delta')
  #Fixing for debug?
  dd_sig = tf.Variable(0.021, trainable=False, name='dd_sig')
  dd_pi = tf.Variable(-0.027, trainable=False, name='dd_pi')
  dd_delta = tf.Variable(-0.027, trainable=False, name='dd_delta')


  SK = tf.stack([ss_sig, sp_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_delta])
  DD = [dd_sig, dd_pi, dd_delta]
  R  = tf.placeholder(shape=[len(features), 10], dtype=tf.float32, name='direction_cosines')
  #Matrix element Wannier
  M_wan  = tf.placeholder(shape=[len(labels)], dtype=tf.float32, name='matrix_elements')
  #Matrix element tb
  M_tb = tf.reduce_sum(tf.multiply(SK,R), axis=1)

  #loss = tf.reduce_sum((M_wan-M_tb)**2)
  loss = tf.reduce_sum(tf.abs(M_wan-M_tb))

  optimizer = tf.train.AdamOptimizer(0.0001)
  #uses re.match on the variable names
  #first_train_vars = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,"dd_")
  first_train_vars = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES)
  #optimize only subspace:
  #train_step = optimizer.minimize(loss, var_list=first_train_vars)
  #optimize full coefficient space:
  train_step = optimizer.minimize(loss, var_list=first_train_vars)
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
        print 'DD-params:', sess.run(SK,feed_dict={R:features, M_wan:labels}), 'loss: ', sess.run(loss, feed_dict=feed_dict)
        sk_out =  sess.run(SK, feed_dict={R:features, M_wan:labels})
        loss_out = sess.run(loss, feed_dict=feed_dict)

# To freeze a subspace.
# train_W = optimizer.minimize(loss, var_list=[W])
# train_C = optimizer.minimize(loss, var_list=[C])
  print '\n'
  print 'TB coefficients', sk_out, 'loss:', loss_out
  print '\n'
  sk_out = np.array(sk_out)

  if verbose:
    print 'lmn coeffs', 'SK Parameters', 'TB Matrix element', 'Wannier Matrix Element'
    for feat, label in zip(features, labels):
      print (' '.join(['{:3.3f}'.format(round(x,4)) for x in feat]), ' '.join(['{:3.3f}'.format(round(x,4)) for x in sk_out]), 
             'SK:', np.round(sk_out.dot(feat),5), 'target:', np.round(label,5))
  return

def matel_coeffs(l,m,n,wan_i,wan_j):
  """
  Takes direction cosines and returns
  the appropriate combination of direction cosines for
  the required energy integral.
  Args:
    l,m,n: direction cosines (p/(p^{2}+q^{2}+r^{2})).
  The orbital indices are wan_i, wan_j  (Energy Integral)
  Subspace:
    [xy, yz, zx]
  Full SK Weight Vector:
  [ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma, pd\sigma, pd\pi, dd\sigma, dd\pi, dd\delta]
  """
  #E_xy,xy
  l2 = np.square(l)
  m2 = np.square(m)
  n2 = np.square(n)
  #Diagonal interactions in {xy, xz,yz} subspace
  if (wan_i == 'xy' and wan_j == 'xy'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l2*m2, l2+m2-4.0*l2*m2, n2+l2*m2]
  elif (wan_i == 'yz' and wan_j == 'yz'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*m2*n2, m2+n2-4.0*m2*n2, l2+m2*n2]
  elif (wan_i == 'zx' and wan_j == 'zx'):
  #Off-Diagonal interactions in {xy, xz,yz} subspace
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*n2*l2, n2+l2-4.0*n2*l2, m2+n2*l2]
  elif (wan_i == 'xy' and wan_j == 'zx'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'xy' and wan_j == 'yz'): 
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'xy'): 
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'yz'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*m*n2*l, m*l*(1.0-4.0*n2), m*l*(n2-1.0)]
  elif (wan_i == 'yz' and wan_j == 'xy'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'yz' and wan_j == 'zx'): 
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*m*n2*l, m*l*(1.0-4.0*n2), m*l*(n2-1.0)]
#  cyclic permutation vs. SK paper. check works!
#  elif wan_i == 'xy' and wan_j='zx':
#  return [3.0*l2*m*n, m*n*(1.0-4.0*l2), m*n*l2-1.0)]
#########################################################
#  Universal Ordering of Slater-Koster parameters:      #
#       [ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma, #
#        pd\sigma, pd\pi,                               #
#        dd\sigma,                                      #
#        dd\pi,                                         #
#        dd\delta]                                      #
#########################################################
#  Off-Diagonal interactions between 6*{sp3d2} and the {xy,xz,yz} subspace.
  elif (wan_i == 'sp3d2-1' and wan_j == 'xy') or (wan_i=='xy' and wan_j=='sp3d2-2'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*l2*m, -m*(1.0-2.0*l2), #pd\sigma  pd\pi #oddparity
           -np.sqrt(3)*l*m*(n2-0.5*(l2+m2))+(3.0/2.0)*l*m*(l2-m2), #dd\sigma
            2.0*np.sqrt(3)*l*m*n2 + 2.0*l*m*(m2-l2), #dd\pi
           -np.sqrt(3.0)/(2.0)*l*m*(1.0+n2)+0.5*l*m*(l2-m2)] #dd\delta
  elif (wan_i=='xy', wan_j=='sp3d2_1') or (wan_i=='sp3d2_2' and wan_j=='xy'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, 
            np.sqrt(3.0)*l2*m, m*(1.0-2.0*l2), #pd\sigma  pd\pi  #oddparity
           -np.sqrt(3.0)*l*m*(n2-0.5*(l2+m2))+(3.0/2.0)*l*m*(l2-m2), #dd\sigma
            2.0*np.sqrt(3)*l*m*n2 + 2.0*l*m*(m2-l2), #dd\pi
           -np.sqrt(3.0)/(2.0)*l*m*(1+n2)+0.5*l*m*(l2-m2)] #dd\delta
  elif (wan_i == 'sp3d2_1' and wan_j == 'zx') or (wan_i=='zx' and wan_j=='sp3d2_2'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*l2*n, -n*(1.0-2.0*l2), #pd\sigma pd\pi #oddparity
            1.0*(np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2))-(3.0/2.0)*n*l(l2-m2)), #dd\sigma
            1.0*(np.sqrt(3.0)*l*n*(l2+m2-n2)-n*l*(1.0-2.0*(l2-m2))), #dd\pi
            1.0*(-0.5*np.sqrt(3.0)*l*n*(l2+m2)+n*l*(1.0-0.5*(l2-m2)))] #dd\delta
  elif (wan_i == 'zx' and wan_j=='sp3d2_1') or (wan_i == 'sp3d2_2' and wan_j=='zx'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*l2*n, n*(1.0-2.0*l2), #pd\sigma pd\pi #oddparity
            1.0*(np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2))-(3.0/2.0)*n*l(l2-m2)), #dd\sigma
            1.0*(np.sqrt(3.0)*l*n*(l2+m2-n2)-n*l*(1.0-2.0*(l2-m2))), #dd\pi
            1.0*(-0.5*np.sqrt(3.0)*l*n*(l2+m2)+n*l*(1.0-0.5*(l2-m2)))] #dd\delta
  elif (wan_i == 'sp3d2_1' and wan_j == 'yz') or (wan_i =='yz' and wan_j== 'sp3d2_2'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*l*m*n, -2.0*l*m*n, #pd\sigma  pd\pi #oddparity
           -np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)) + (3.0/2.0)*m*n*(l2-m2), #dd\sigma
           -np.sqrt(3.0)*m*n*(l2+m2-n2) - m*n*(1.0+2.0*(l2-m2)), #dd\pi
           -0.5*np.sqrt(3.0)*m*n*(l2+m2) + m*n*(1.0+0.5*(l2-m2))] #dd\delta
  elif (wan_i == 'yz' and wan_j == 'sp3d2_1') or ( wan_j=='sp3d2_2' and wan_i=='yz'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*l*m*n, 2.0*l*m*n, #pd\sigma  pd\pi #odd parity
           -np.sqrt(3)*m*n*(n2-0.5*(l2+m2)) + (3.0/2.0)*m*n*(l2-m2), #dd\sigma
           -np.sqrt(3)*m*n*(l2+m2-n2) - m*n*(1.0+2.0*(l2-m2)), #dd\pi
           -0.5*np.sqrt(3)*m*n*(l2+m2) + m*n*(1.0+0.5*(l2-m2))] #dd\delta
  elif (wan_i == 'sp3d2_3' and wan_j == 'xy') or (wan_i=='xy' and wan_j =='sp3d2_4'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*m2*l, -l*(1.0-2.0*m2), #pd\sigma  pd\pi #oddparity
           -np.sqrt(3.0)*l*m*(n2-0.5*(l2+m2)) - (3.0/2.0)*l*m*(l2-m2), #dd\sigma
            2.0*np.sqrt(3.0)*l*m*n2 - 2.0*l*m*(m2-l2), #dd\pi
           -np.sqrt(3.0)/(2.0)*l*m*(1.0+n2) - 0.5*l*m*(l2-m2)] #dd\delta
  elif (wan_i=='xy' and wan_j=='sp3d2_3') or (wan_i=='sp3d2_4' and wan_j=='xy'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, 
            np.sqrt(3.0)*m2*l, l*(1.0-2.0*m2), #pd\sigma  pd\pi #oddparity
           -np.sqrt(3)*l*m*(n2-0.5*(l2+m2)) - (3.0/2.0)*l*m*(l2-m2), #dd\sigma
            2.0*np.sqrt(3)*l*m*n2 - 2.0*l*m*(m2-l2), #dd\pi
           -np.sqrt(3.0)/(2.0)*l*m*(1.0+n2) - 0.5*l*m*(l2-m2)] #dd\delta
  elif (wan_i == 'sp3d2_3' and wan_j == 'zx') or (wan_i=='zx' and wan_j=='sp3d2_4'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*m*n*l, -2.0*m*n*l, #pd\sigma pd\pi #oddparity
            1.0*(np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2)) + (3.0/2.0)*n*l(l2-m2)), #dd\sigma
            1.0*(np.sqrt(3.0)*l*n*(l2+m2-n2) + n*l*(1.0-2.0*(l2-m2))), #dd\pi
            1.0*(-0.5*np.sqrt(3.0)*l*n*(l2+m2) - n*l*(1.0-0.5*(l2-m2)))] #dd\delta
  elif (wan_i == 'zx' and wan_j=='sp3d2_3') or (wan_i=='sp3d2_4' and wan_j=='zx'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            -np.sqrt(3.0)*m*n*l, 2.0*m*n*l, #pd\sigma pd\pi #oddparity
            1.0*(np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2)) + (3.0/2.0)*n*l(l2-m2)), #dd\sigma
            1.0*(np.sqrt(3.0)*l*n*(l2+m2-n2) + n*l*(1.0-2.0*(l2-m2))), #dd\pi
            1.0*(-0.5*np.sqrt(3.0)*l*n*(l2+m2) - n*l*(1.0-0.5*(l2-m2)))] #dd\delta
  elif (wan_i == 'sp3d2_3' and wan_j == 'yz') or (wan_i =='sp3d2_4' and wan_j=='yz'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*m2*n, -n*(1.0-2.0*m2), #pd\sigma  pd\pi #oddparity
           -np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)) - (3.0/2.0)*m*n*(l2-m2), #dd\sigma
           -np.sqrt(3.0)*m*n*(l2+m2-n2) + m*n*(1.0+2.0*(l2-m2)), #dd\pi
           -0.5*np.sqrt(3.0)*m*n*(l2+m2) - m*n*(1.0+0.5*(l2-m2))] #dd\delta
  elif (wan_i == 'yz' and wan_j == 'sp3d2_3') or (wan_i =='sp3d2_4'and wan_i =='yz'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*m2*n, n*(1.0-2.0*m2), #pd\sigma  pd\pi #oddparity
           -np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)) - (3.0/2.0)*m*n*(l2-m2), #dd\sigma
           -np.sqrt(3.0)*m*n*(l2+m2-n2) + m*n*(1.0+2.0*(l2-m2)), #dd\pi
           -0.5*np.sqrt(3.0)*m*n*(l2+m2) - m*n*(1.0+0.5*(l2-m2))] #dd\delta
  elif (wan_i == 'sp3d2_5' and wan_j == 'xy') or (wan_i =='xy' and wan_j=='sp3d2_6'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -np.sqrt(3.0)*l*m*n, 2.0*n*m*l, #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*l*m*(n2-0.5*(l2+m2)), #dd\sigma
           -2.0*np.sqrt(3.0)*l*m*n2, #dd\pi
            0.5*np.sqrt(3.0)*l*m*(1.0+n2)] #dd\delta
  elif (wan_i == 'sp3d2_6' and wan_j == 'xy') or (wan_i =='xy' and wan_j=='sp3d2_5'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*l*m, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*l*m*n, -2.0*n*m*l, #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*l*m*(n2-0.5*(l2+m2)), #dd\sigma
           -2.0*np.sqrt(3.0)*l*m*n2, #dd\pi
            0.5*np.sqrt(3.0)*l*m*(1.0+n2)] #dd\delta
  elif (wan_i == 'sp3d2_5' and wan_j == 'zx') or (wan_i =='zx' and wan_j=='sp3d2_6'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*n2*l, -l*(1.0-2.0*n2), #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2)), #dd\sigma
            np.sqrt(3.0)*l*n*(l2+m2-n2), #dd\pi
            -0.5*np.sqrt(3)*l*n*(l2+m2)] #dd\delta
  elif (wan_i == 'sp3d2_6' and wan_j == 'zx') or (wan_i =='zx' and wan_j=='sp3d2_5'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*n*l, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            -np.sqrt(3.0)*n2*l, l*(1.0-2.0*n2), #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2)), #dd\sigma
            np.sqrt(3.0)*l*n*(l2+m2-n2), #dd\pi
            -0.5*np.sqrt(3)*l*n*(l2+m2)] #dd\delta
  elif (wan_i == 'sp3d2_5' and wan_j == 'yz') or (wan_i =='yz' and wan_j=='sp3d2_6'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            -np.sqrt(3)*n2*m, -m*(1.0-2.0*n2), #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)), #dd\sigma
            np.sqrt(3.0)*m*n*(l2+m2-n2), #dd\pi
            -0.5*np.sqrt(3.0)*m*n*(l2+m2)] #dd\delta
  elif (wan_i == 'sp3d2_6' and wan_j == 'yz') or (wan_i =='yz' and wan_j=='sp3d2_5'):
    return [0.0, 0.0, 0.0, 0.0, np.sqrt(3.0)*m*n, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            np.sqrt(3.0)*n2*m, m*(1.0-2.0*n2), #pd\sigma  pd\pi #oddparity
            np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)), #dd\sigma
            np.sqrt(3.0)*m*n*(l2+m2-n2), #dd\pi
            -0.5*np.sqrt(3.0)*m*n*(l2+m2)] #dd\delta
  else:
    sys.exit('Invalid Combo')

def load_matels(l,m,n,subspace_matrix):
  """
  Takes direction cosines, and 3x3 subspace matrix.
  """
  matel_dict = {}
  dxy_subspace = [('xy', 'xy'), ('yz', 'yz'), ('zx', 'zx'), ('xy', 'zx'), ('xy', 'yz'),
                  ('zx', 'xy'), ('zx', 'yz'), ('yz', 'xy'), ('yz', 'zx')]
  #diagonal
  matel_dict[('xy', 'xy')] = (6, 6)
  matel_dict[('zx', 'zx')] = (7, 7)
  matel_dict[('yz', 'yz')] = (8, 8)
  #off diagonal
  matel_dict[('xy', 'zx')] = (6, 7)
  matel_dict[('xy', 'yz')] = (6, 8)
  matel_dict[('zx', 'xy')] = (7, 6)
  matel_dict[('zx', 'yz')] = (7, 8)
  matel_dict[('yz', 'xy')] = (8, 6)
  matel_dict[('yz', 'zx')] = (8, 7)

  sp3d2_dxy = [('sp3d2-1', 'xy'), ('xy', 'sp3d2-1'), 
               ('sp3d2-1', 'zx'), ('zx', 'sp3d2-1'),
               ('sp3d2-1', 'yz'), ('yz', 'sp3d2-1'),
               ('sp3d2-2', 'xy'), ('xy', 'sp3d2-2'), 
               ('sp3d2-2', 'zx'), ('zx', 'sp3d2-2'),
               ('sp3d2-2', 'yz'), ('yz', 'sp3d2-2'),
               ('sp3d2-3', 'xy'), ('xy', 'sp3d2-3'), 
               ('sp3d2-3', 'zx'), ('zx', 'sp3d2-3'),
               ('sp3d2-3', 'yz'), ('yz', 'sp3d2-3'),
               ('sp3d2-4', 'xy'), ('xy', 'sp3d2-4'), 
               ('sp3d2-4', 'zx'), ('zx', 'sp3d2-4'),
               ('sp3d2-4', 'yz'), ('yz', 'sp3d2-4'),
               ('sp3d2-5', 'xy'), ('xy', 'sp3d2-5'), 
               ('sp3d2-5', 'zx'), ('zx', 'sp3d2-5'),
               ('sp3d2-5', 'yz'), ('yz', 'sp3d2-5'),
               ('sp3d2-6', 'xy'), ('xy', 'sp3d2-6'), 
               ('sp3d2-6', 'zx'), ('zx', 'sp3d2-6'),
               ('sp3d2-6', 'yz'), ('yz', 'sp3d2-6')]

  matel_dict[('sp3d2-1', 'xy')], matel_dict[('xy', 'sp3d2-1')] = (0,6), (6,0)
  matel_dict[('sp3d2-1', 'zx')], matel_dict[('zx', 'sp3d2-1')] = (0,7), (7,0)
  matel_dict[('sp3d2-1', 'yz')], matel_dict[('yz', 'sp3d2-1')] = (0,8), (8,0)
  matel_dict[('sp3d2-2', 'xy')], matel_dict[('xy', 'sp3d2-2')] = (1,6), (6,1)
  matel_dict[('sp3d2-2', 'zx')], matel_dict[('zx', 'sp3d2-2')] = (1,7), (7,1)
  matel_dict[('sp3d2-2', 'yz')], matel_dict[('yz', 'sp3d2-2')] = (1,8), (8,1)
  matel_dict[('sp3d2-3', 'xy')], matel_dict[('xy', 'sp3d2-3')] = (2,6), (6,2)
  matel_dict[('sp3d2-3', 'zx')], matel_dict[('zx', 'sp3d2-3')] = (2,7), (7,2)
  matel_dict[('sp3d2-3', 'yz')], matel_dict[('yz', 'sp3d2-3')] = (2,8), (8,2)
  matel_dict[('sp3d2-4', 'xy')], matel_dict[('xy', 'sp3d2-4')] = (3,6), (6,3)
  matel_dict[('sp3d2-4', 'zx')], matel_dict[('zx', 'sp3d2-4')] = (3,7), (7,3)
  matel_dict[('sp3d2-4', 'yz')], matel_dict[('yz', 'sp3d2-4')] = (3,8), (8,3)
  matel_dict[('sp3d2-5', 'xy')], matel_dict[('xy', 'sp3d2-5')] = (4,6), (6,4)
  matel_dict[('sp3d2-5', 'zx')], matel_dict[('zx', 'sp3d2-5')] = (4,7), (7,4)
  matel_dict[('sp3d2-5', 'yz')], matel_dict[('yz', 'sp3d2-5')] = (4,8), (8,4)
  matel_dict[('sp3d2-6', 'xy')], matel_dict[('xy', 'sp3d2-6')] = (5,6), (6,5)
  matel_dict[('sp3d2-6', 'zx')], matel_dict[('zx', 'sp3d2-6')] = (5,7), (7,5)
  matel_dict[('sp3d2-6', 'yz')], matel_dict[('yz', 'sp3d2-6')] = (5,8), (8,5)

  features = []
  matels = []
  #for pair in dxy_subspace:
  for pair in sp3d2_dxy + dxy_subspace:
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
  #subspace_matrix = Fe_nn[-3:,-3:].real
  #use full matrix
  subspace_matrix = Fe_nn[:,:].real
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

