import sys

import pickle
import numpy  as np
import tensorflow as tf

from ase import units

from pandas import DataFrame

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


def print_predictions(pairs, predictions, matel_dict):
  """
  Pretty printing function. Takes pairs in orbital notation,
  writes them to the
  """
  pred_mat = np.zeros([len(pairs)])
  for pair, prediction in zip(pairs, predictions):
    pred_mat[matel_dict[pair]] = prediction
  print DataFrame(pred_mat)
  #print('\n'.join([''.join(['{:1.4f}'.format(item) for item in row]) for row in pred_mat]))

def fit_tb(features, labels, pairs, verbose=True):
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
  dd_sig = tf.Variable(0.0, trainable=True, name='dd_sig')
  dd_pi = tf.Variable(0.0, trainable=True, name='dd_pi')
  dd_delta = tf.Variable(0.0, trainable=True, name='dd_delta')
  #Fixing for debug?
  #dd_sig = tf.Variable(0.021, trainable=False, name='dd_sig')
  #dd_pi = tf.Variable(-0.027, trainable=False, name='dd_pi')
  #dd_delta = tf.Variable(-0.027, trainable=False, name='dd_delta')


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
    for pair, feat, label in zip(pairs, features, labels):
      print ' '.join(['{:3.3f}'.format(round(x,5)) for x in feat])
      print ' '.join(['{:3.3f}'.format(round(x,5)) for x in sk_out])
      print pair
      print 'SK:', np.round(sk_out.dot(feat),5), 'target:', np.round(label,5)
      print '\n'
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
  nsqrt2 = 1.0/np.sqrt(2.0)
  nsqrt3 = 1.0/np.sqrt(3.0)
  nsqrt6 = 1.0/np.sqrt(6.0)
  nsqrt12 = 1.0/np.sqrt(12.0)

#Empirical remapping: maybe states don't correspond to users guide? Should wok for 2nd NN.
#remap_dict = {'sp3d2-1':'sp3d2-5', 'sp3d2-2':'sp3d2-6','sp3d2-3':'sp3d2-1',
#              'sp3d2-4':'sp3d2-2','sp3d2-5':'sp3d2-3','sp3d2-6':'sp3d2-4'}
#remap_dict = {'yz': 'zx','zx':'yz'}
#if wan_i in remap_dict.keys():
#  wan_i = remap_dict[wan_i]
#  print 'remap'
#if wan_j in remap_dict.keys():
#  wan_j = remap_dict[wan_j]
#grouped together in cyclic permutation groups from SK tables
###########################################
  sdsig_s_xy = -nsqrt6*np.sqrt(3)*l*m
  sdsig_s_yz = nsqrt6*np.sqrt(3)*m*n
  sdsig_s_zx = -nsqrt6*np.sqrt(3)*n*l
###########################################

#trying to address sign problem when I permute the SK table
#similar to the xy,zx matrix elements derived in the table for dxy subspace
#still don't fully understand this. (swapping n/l?)
###########################################

  pdsig_x_yz = -nsqrt2*np.sqrt(3.0)*l*m*n
  pdsig_y_zx = -nsqrt2*np.sqrt(3.0)*m*n*l
  pdsig_z_xy = -nsqrt2*np.sqrt(3.0)*l*m*n

  pdpi_x_yz = nsqrt2*2.0*l*m*n
  pdpi_y_zx = nsqrt2*2.0*m*n*l
  pdpi_z_xy = nsqrt2*2.0*n*l*m

  pdsig_x_zx = nsqrt2*np.sqrt(3.0)*l2*n
###########################################
  pdsig_y_xy = -nsqrt2*np.sqrt(3.0)*m2*l
###########################################
  pdsig_z_yz = nsqrt2*np.sqrt(3)*n2*m

  pdpi_x_zx = nsqrt2*n*(1.0-2.0*l2)
###########################################
  pdpi_y_xy = -nsqrt2*l*(1.0-2.0*m2)
###########################################
  pdpi_z_yz = nsqrt2*m*(1.0-2.0*n2)

  pdsig_x_xy = nsqrt2*np.sqrt(3.0)*l2*m
  pdsig_y_yz = nsqrt2*np.sqrt(3.0)*m2*n 
###########################################
  pdsig_z_zx = -nsqrt2*np.sqrt(3.0)*l2*n
###########################################

  pdpi_x_xy = nsqrt2*m*(1.0-2.0*l2)
  pdpi_y_yz = nsqrt2*n*(1.0-2.0*m2)
###########################################
  pdpi_z_zx = -nsqrt2*l*(1.0-2.0*n2)
###########################################

  ddsig_xy_z2 = nsqrt12*np.sqrt(3.0)*l*m*(n2-0.5*(l2+m2))
  ddpi_xy_z2 = nsqrt12*2.0*np.sqrt(3.0)*l*m*n2 
  dddelta_xy_z2 = -nsqrt12*np.sqrt(3.0)/(2.0)*l*m*(1.0+n2) 

  ddsig_xy_x2my2 = 0.5*(3.0/2.0)*l*m*(l2-m2)
  ddpi_xy_x2my2 = 0.5*2.0*l*m*(m2-l2)
  ddelta_xy_x2my2 = 0.5*0.5*l*m*(l2-m2)

  ddsig_zx_z2 = nsqrt12*np.sqrt(3.0)*l*n*(n2-0.5*(l2+m2))
  ddpi_zx_z2 = nsqrt12*np.sqrt(3.0)*l*n*(l2+m2-n2) 
  dddelta_zx_z2 = -nsqrt12*0.5*np.sqrt(3.0)*l*n*(l2+m2)

  ddsig_zx_x2my2 = 0.5*(3.0/2.0)*n*l*(l2-m2)
  ddpi_zx_x2my2 = 0.5*n*l*(1.0-2.0*(l2-m2))
  dddelta_zx_x2my2 = -0.5*n*l*(1.0-0.5*(l2-m2))

  ddsig_yz_z2 = nsqrt12*np.sqrt(3.0)*m*n*(n2-0.5*(l2+m2)) 
  ddpi_yz_z2 = nsqrt12*np.sqrt(3.0)*m*n*(l2+m2-n2)
  dddelta_yz_z2 = -nsqrt12*0.5*np.sqrt(3.0)*m*n*(l2+m2)

  ddsig_yz_x2my2 = 0.5*(3.0/2.0)*m*n*(l2-m2)
  ddpi_yz_x2my2 = -0.5*m*n*(1.0+2.0*(l2-m2))
  dddelta_yz_x2my2 = 0.5*m*n*(1.0+0.5*(l2-m2))

  sssig = 1.0
  spsig_s_x = l
  spsig_s_y = m
  spsig_s_z = n
  sdsig_s_z2 = (n2 - 0.5*(l2+m2))
  sdsig_s_x2my2 = 0.5*np.sqrt(3.0)*(l2-m2)
  
  ppsig_x_x = l2
  ppsig_y_y = m2
  ppsig_z_z = n2
  
  pppi_x_x = (1.0-l2)
  pppi_y_y = (1.0-m2)
  pppi_z_z = (1.0-n2)
  
  ppsig_x_y = l*m
  ppsig_y_z = m*n
  ppsig_z_x = n*l

  pppi_x_y = -l*m
  pppi_y_z = -m*n
  pppi_z_x = -n*l

  ppsig_x_z = l*n
  ppsig_y_x = m*l
  ppsig_z_x = n*m

  pppi_x_z = -l*n
  pppi_y_x = -m*l
  pppi_z_x = -n*m

#Normalizations for diagonal elements of sp3d2 subspace
#From W90 user_guide.
  norm_s_s = (1.0/6.0)
  norm_s_p = nsqrt6*nsqrt2
  norm_s_dz2 = nsqrt12*nsqrt6
  norm_s_dx2y2 = nsqrt2*nsqrt6
  norm_p_p = (1.0/2.0)
  norm_spxyd2_dz2_dz2 = (1.0/12.0)
  norm_spxyd2_dx2y2_dx2y2 = (1.0/4.0)
  norm_spxyd2_dz2_dx2y2 = nsqrt12
#This normalization applied in the conditional part below.
  norm_spzd2_dz_dz = (1.0/3.0)
#Combined terms (1,2), (2,1)
  spxd2_spxd2_ddsigma = norm_spxyd2_dz2_dz2*(n2-(1.0/2.0)*(l2+m2)**2) + \
                       -1.0*norm_spxyd2_dz2_dx2y2*(np.sqrt(3.0)/2.0)*(l2-m2)*(n2-0.5*(l2+m2)) + \
                        norm_spxyd2_dx2y2_dx2y2*((3.0/4.0)*(l2-m2)**2) 

  spxd2_spxd2_ddpi = norm_spxyd2_dz_dz*(3.0*n2*(l2+m2)) + \
                -1.0*norm_spxyd2_dz_dx2y2*(np.sqrt(3.0)*n2*(m2-l2)) + \
                     norm_spxyd2_dx2y2_dx2y2*(l2+m2-(l2-m2)**2) 

  spxd2_spxd2_dddelta = norm_spxyd2_dz_dz*(3.0/4.0)*(l2+m2) + \
                   -1.0*norm_spxyd2_dz_dx2y2*((1.0/4.0)*np.sqrt(3.0)*(1.0+n2)*(l2-m2)) +\
                        norm_spxyd2_dx2y2_dx2y2*(n2+(1.0/4.0)**(l2-m2)**2)  
#Combined terms (3,4), (4,3) 
  spyd2_spyd2_ddsigma = norm_spxyd2_dz2_dz2*(n2-(1.0/2.0)*(l2+m2)**2) + \
                        norm_spxyd2_dz2_dx2y2*(np.sqrt(3.0)/2.0)*(l2-m2)*(n2-0.5*(l2+m2)) + \
                        norm_spxyd2_dx2y2_dx2y2*((3.0/4.0)*(l2-m2)**2) 

  spyd2_spyd2_ddpi = norm_spxyd2_dz_dz*(3.0*n2*(l2+m2)) + \
                     norm_spxyd2_dz_dx2y2*(np.sqrt(3.0)*n2*(m2-l2)) + \
                     norm_spxyd2_dx2y2_dx2y2*(l2+m2-(l2-m2)**2)

  spyd2_spyd2_dddelta = norm_spxyd2_dz_dz*(3.0/4.0)*(l2+m2) + \
                        norm_spxyd2_dz_dx2y2*((1.0/4.0)*np.sqrt(3.0)*(1.0+n2)*(l2-m2)) +\
                        norm_spxyd2_dx2y2_dx2y2*(n2+(1.0/4.0)**(l2-m2)**2)

#Combined terms (1,3), (3,1) and (1,4) (4,1) no (dz2|dx2y2) cross terms.
  spxyd2_spxyd2_ddsigma_od = norm_spxyd2_dz_dz*(n2-0.5*(l2+m2)**2) -\
                             norm_spxyd2_dx2y2_dx2y2*((3.0/4.0)*(l2-m2)**2)

  spxyd2_spxyd2_ddpi_od = norm_spxyd2_dz_dz*(3.0*n2*(l2+m2)) - \
                             norm_spxyd2_dx2y2_dx2y2*(l2+m2-(l2-m2)**2)

  spxyd2_spxyd2_dddelta_od = norm_spxyd2_dz_dz*((3.0/4.0)*(l2+m2)) - \
                             norm_spxyd2_dx2y2_dx2y2*(n2+(1.0/4.0)*(l2-m2)**2)

  ddsig_z2_z2 = n2-0.5*(l2+m2) 
  ddpi_z2_z2 = 3.0*n2*(l2+m2)  
  dddelta_z2_z2 = 0.75*(l2+m2) 

  norm_pxy_dz2 = nsqrt2*nsqrt12
  norm_pz_dz2 = nsqrt2*nsqrt3
  norm_pxy_dx2y2 = nsqrt2*(0.5)

  px_d2_pdsig = norm_pxy_dz2*(l*(n2-0.5*(l2+m2))) + norm_pxy_dx2y2*(0.5*np.sqrt(3.0)*l*(l2-m2))
  px_d2_pdpi = -norm_pxy_dz2*(np.sqrt(3.0)*l*n2) + norm_pxy_dx2y2*(l*(1.0-l2+m2))
  py_d2_pdsig = norm_pxy_dz2*(m*(n2-0.5*(l2+m2))) + norm_pxy_dx2y2*(0.5*np.sqrt(3.0)*m*(l2-m2))
  py_d2_pdpi = -norm_pxy_dz2*(np.sqrt(3.0)*m*n2) - norm_pxy_dx2y2*(m*(1.0+l2-m2))
  pz_d2_pdsig = norm_pz_dz2*(n*(n2-0.5*(l2+m2)))
  pz_d2_pdpi = -norm_pz_dz2*(np.sqrt(3.0)*m*n2)

#p/d terms occurring for off diagonal elements (1,3), (1,4) py have odd terms.
  px_d2_pdsig_od = norm_pxy_dz2*(l*(n2-0.5*(l2+m2))) + norm_pxy_dx2y2*((1.0/2.0)*np.sqrt(3.0)*l*(l2-m2))
  px_d2_pdpi_od = -norm_pxy_dz2*(np.sqrt(3.0)*l*n2) + norm_pxy_dx2y2*(l*(1.0-l2+m2))

  py_d2_pdsig_od = -norm_pxy_dz2*(m*(n2-0.5*(l2+m2))) + norm_pxy_dx2y2*(0.5*np.sqrt(3.0)*m*(l2-m2))
  py_d2_pdpi_od = norm_pxy_dz2*(np.sqrt(3.0)*m*n2) - norm_pxy_dx2y2*(m*(1.0+l2-m2))

#(1,5) and (1,6) pd terms.
#"cool" normalization factor:
  cool_norm = (np.sqrt(2.0)*(2.0-np.sqrt(12.0)))/(4.0*np.sqrt(12.0))

  pxz_dz2_pdsig = (-1.0)/(np.sqrt(2)*np.sqrt(3))*(l*(n2-0.5*(l2+m2))) + cool_norm*(n*(n2-0.5*l2+m2))
  pxz_dz2_pdpi = -0.5*l*n2 + cool_norm*np.sqrt(3) 
  
  sp3d2xz_d2_ddsig = -1.0/(nsqrt3*nsqrt12)*(n2-0.5*(l2+m2)**2) + \
                     (1.0/(2.0*nsqrt(3.0)))*((3.0/4.0)*(l2-m2)**2)

  sp3d2xz_d2_ddpi = -1.0/(nsqrt3*nsqrt12)*(3.0*n2*(l2+m2)) + \
                   (1.0/(2.0*nsqrt(3.0)))*(l2+m2-(l2-m2)**2)

  sp3d2xz_d2_dddelta = -1.0/(nsqrt3*nsqrt12)*((3.0/4.0)*(l2+m2)) + \
                     (1.0/(2.0*nsqrt(3.0)))*(n2+(1.0/4.0)*(l2-m2)**2)

  #Diagonal interactions in {xy, xz,yz} subspace
  if (wan_i == 'xy' and wan_j == 'xy'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l2*m2, l2+m2-4.0*l2*m2, n2+l2*m2]
  elif (wan_i == 'yz' and wan_j == 'yz'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*m2*n2, m2+n2-4.0*m2*n2, l2+m2*n2]
  elif (wan_i == 'zx' and wan_j == 'zx'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*n2*l2, n2+l2-4.0*n2*l2, m2+n2*l2]
  #Off-Diagonal interactions in {xy, xz,yz} subspace
  elif (wan_i == 'xy' and wan_j == 'zx'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'xy'): 
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.0*n*l2*m, -n*m*(1.0-4.0*l2), -n*m*(l2-1.0)]
  elif (wan_i == 'xy' and wan_j == 'yz'): 
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'yz' and wan_j == 'xy'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*l*m2*n, l*n*(1.0-4.0*m2), l*n*(m2-1.0)]
  elif (wan_i == 'zx' and wan_j == 'yz'):
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0*m*n2*l, m*l*(1.0-4.0*n2), m*l*(n2-1.0)]
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
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, 
           -pdsig_x_xy, -pdpi_x_xy,
           -ddsig_xy_z2 + ddsig_xy_x2my2, 
           -ddpi_xy_z2 + ddpi_xy_x2my2, 
           -dddelta_xy_z2 + ddelta_xy_x2my2]
  elif (wan_i=='xy' and wan_j=='sp3d2-1') or (wan_i=='sp3d2-2' and wan_j=='xy'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, 
            pdsig_x_xy, pdpi_x_xy, 
           -ddsig_xy_z2 + ddsig_xy_x2my2, 
           -ddpi_xy_z2 + ddpi_xy_x2my2, 
           -dddelta_xy_z2 + ddelta_xy_x2my2] 
  elif (wan_i == 'sp3d2-1' and wan_j == 'zx') or (wan_i=='zx' and wan_j=='sp3d2-2'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, 
           -pdsig_x_zx,  -pdpi_x_zx, 
           -ddsig_zx_z2 + ddsig_zx_x2my2, 
           -ddpi_zx_z2 + ddpi_zx_x2my2, 
           -dddelta_zx_z2 + dddelta_zx_x2my2]
  elif (wan_i == 'zx' and wan_j=='sp3d2-1') or (wan_i == 'sp3d2-2' and wan_j=='zx'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, 
            pdsig_x_zx, pdpi_x_zx, 
            -ddsig_zx_z2 + ddsig_zx_x2my2,
            -ddpi_zx_z2 + ddpi_zx_x2my2, 
            -dddelta_zx_z2 + dddelta_zx_x2my2] 
  elif (wan_i == 'sp3d2-1' and wan_j == 'yz') or (wan_i =='yz' and wan_j== 'sp3d2-2'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -pdsig_x_yz, -pdpi_x_yz, #pd\sigma  pd\pi #oddparity
           -ddsig_yz_z2 + ddsig_yz_x2my2, #dd\sigma
           -ddpi_yz_z2 + ddpi_yz_x2my2 , #dd\pi
           -dddelta_yz_z2 + dddelta_yz_x2my2] #dd\delta
  elif (wan_i == 'yz' and wan_j == 'sp3d2-1') or ( wan_i=='sp3d2-2' and wan_j=='yz'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, 
            pdsig_x_yz, pdpi_x_yz, 
           -ddsig_yz_z2 + ddsig_yz_x2my2, 
           -ddpi_yz_z2 + ddpi_yz_x2my2 , 
           -dddelta_yz_z2 + dddelta_yz_x2my2] 
  elif (wan_i == 'sp3d2-3' and wan_j == 'xy') or (wan_i=='xy' and wan_j =='sp3d2-4'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -pdsig_y_xy, -pdpi_y_xy, #pd\sigma  pd\pi #oddparity
           -ddsig_xy_z2 - ddsig_xy_x2my2, 
           -ddpi_xy_z2 - ddpi_xy_x2my2, 
           -dddelta_xy_z2 - ddelta_xy_x2my2]
  elif (wan_i=='xy' and wan_j=='sp3d2-3') or (wan_i=='sp3d2-4' and wan_j=='xy'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, 
            pdsig_y_xy, pdpi_y_xy, #pd\sigma  pd\pi #oddparity
           -ddsig_xy_z2 - ddsig_xy_x2my2, 
           -ddpi_xy_z2 - ddpi_xy_x2my2, 
           -dddelta_xy_z2 - ddelta_xy_x2my2]
  elif (wan_i == 'sp3d2-3' and wan_j == 'zx') or (wan_i=='zx' and wan_j=='sp3d2-4'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
           -pdsig_y_zx, -pdpi_y_zx, #pd\sigma pd\pi #oddparity
           -ddsig_zx_z2 - ddsig_zx_x2my2, 
           -ddpi_zx_z2 - ddpi_zx_x2my2, 
           -dddelta_zx_z2 - dddelta_zx_x2my2]
  elif (wan_i == 'zx' and wan_j=='sp3d2-3') or (wan_i=='sp3d2-4' and wan_j=='zx'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, #ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma,
            pdsig_y_zx, pdpi_y_zx, #pd\sigma pd\pi #oddparity
           -ddsig_zx_z2 - ddsig_zx_x2my2, 
           -ddpi_zx_z2 - ddpi_zx_x2my2, 
           -dddelta_zx_z2 - dddelta_zx_x2my2]
  elif (wan_i == 'sp3d2-3' and wan_j == 'yz') or (wan_i =='yz' and wan_j=='sp3d2-4'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, 
           -pdsig_y_yz, -pdpi_y_yz, 
           -ddsig_yz_z2 - ddsig_yz_x2my2, 
           -ddpi_yz_z2 - ddpi_yz_x2my2, 
           -dddelta_yz_z2 - dddelta_yz_x2my2] 
  elif (wan_i == 'yz' and wan_j == 'sp3d2-3') or (wan_i =='sp3d2-4' and wan_j=='yz'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, 
            pdsig_y_yz, pdpi_y_yz, 
           -ddsig_yz_z2 - ddsig_yz_x2my2, 
           -ddpi_yz_z2 - ddpi_yz_x2my2 , 
           -dddelta_yz_z2 - dddelta_yz_x2my2] 
  elif (wan_i == 'sp3d2-5' and wan_j == 'xy') or (wan_i =='xy' and wan_j=='sp3d2-6'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, 
           -pdsig_z_xy, -pdpi_z_xy, 
            2.0*ddsig_xy_z2,    
            2.0*ddpi_xy_z2,    
            2.0*dddelta_xy_z2]
  elif (wan_i =='xy' and wan_j=='sp3d2-5') or (wan_i == 'sp3d2-6' and wan_j == 'xy'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_xy, 
            pdsig_z_xy, pdpi_z_xy, 
            2.0*ddsig_xy_z2,    
            2.0*ddpi_xy_z2,    
            2.0*dddelta_xy_z2]
  elif (wan_i == 'sp3d2-5' and wan_j == 'zx') or (wan_i =='zx' and wan_j=='sp3d2-6'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, 
           -pdsig_z_zx, -pdpi_z_zx, 
            2.0*ddsig_zx_z2, 
            2.0*ddpi_zx_z2, 
            2.0*dddelta_zx_z2]
  elif (wan_i == 'sp3d2-6' and wan_j == 'zx') or (wan_i =='zx' and wan_j=='sp3d2-5'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_zx, 
            pdsig_z_zx, pdpi_z_zx, 
            2.0*ddsig_zx_z2, 
            2.0*ddpi_zx_z2, 
            2.0*dddelta_zx_z2]
  elif (wan_i == 'sp3d2-5' and wan_j == 'yz') or (wan_i =='yz' and wan_j=='sp3d2-6'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, 
           -pdsig_z_yz, -pdpi_z_yz, 
            2.0*ddsig_yz_z2,   
            2.0*ddpi_yz_z2,    
            2.0*dddelta_yz_z2] 
  elif (wan_i == 'sp3d2-6' and wan_j == 'yz') or (wan_i =='yz' and wan_j=='sp3d2-5'):
    return [0.0, 0.0, 0.0, 0.0, sdsig_s_yz, 
            pdsig_z_yz, pdpi_z_yz, 
            2.0*ddsig_yz_z2,   
            2.0*ddpi_yz_z2,    
            2.0*dddelta_yz_z2] 
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-1'):
    return [1.0*norm_s_s, 0.0, norm_p_p*ppsig_x_x, norm_p_p*pppi_x_x, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            0.0, 0.0,
            spxyd2_spxyd2_ddsigma,
            spxyd2_spxyd2_ddpi,
            spxyd2_spxyd2_dddelta]
  elif (wan_i == 'sp3d2-2' and wan_j == 'sp3d2-2'):
    return [1.0*norm_s_s, 0.0, norm_p_p*ppsig_x_x, norm_p_p*pppi_x_x, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            0.0,0.0,
            spxyd2_spxyd2_ddsigma,
            spxyd2_spxyd2_ddpi,
            spxyd2_spxyd2_dddelta]
  elif (wan_i == 'sp3d2-3' and wan_j == 'sp3d2-3'):
    return [norm_s_s, 0.0, norm_p_p*ppsig_y_y, norm_p_p*pppi_y_y, -2.0*(norm_s_dz2*sdsig_s_z2 - norm_s_dx2y2*sdsig_s_x2my2), 
            0.0,0.0,
            spxyd2_spxyd2_ddsigma,
            spxyd2_spxyd2_ddpi,
            spxyd2_spxyd2_dddelta]
  elif (wan_i == 'sp3d2-4' and wan_j == 'sp3d2-4'):
    return [1.0*norm_s_s, 0.0, norm_p_p*ppsig_y_y, norm_p_p*pppi_y_y, -2.0*(norm_s_dz2*sdsig_s_z2 - norm_s_dx2y2*sdsig_s_x2my2), 
            0.0,0.0,
            spxyd2_spxyd2_ddsigma,
            spxyd2_spxyd2_ddpi,
            spxyd2_spxyd2_dddelta]
  elif (wan_i == 'sp3d2-5' and wan_j == 'sp3d2-5'):
    return [1.0*norm_s_s, 0.0, norm_p_p*ppsig_z_z, norm_p_p*pppi_z_z, -2.0*(norm_s_dz2*sdsig_s_z2), 
            0.0,0.0,
            norm_spzd2_dz_dz*ddsig_z2_z2,
            norm_spzd2_dz_dz*ddpi_z2_z2,
            norm_spzd2_dz_dz*dddelta_z2_z2]
  elif (wan_i == 'sp3d2-6' and wan_j == 'sp3d2-6'):
    return [1.0*norm_s_s, 0.0, norm_p_p*ppsig_z_z, norm_p_p*pppi_z_z, -2.0*(norm_s_dz2*sdsig_s_z2), 
            0.0,0.0,
            norm_spzd2_dz_dz*ddsig_z2_z2,
            norm_spzd2_dz_dz*ddpi_z2_z2,
            norm_spzd2_dz_dz*dddelta_z2_z2]
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-2'):
    return [1.0*norm_s_s, 2.0*norm_s_p*spsig_s_x, -norm_p_p*ppsig_x_x, -norm_p_p*pppi_x_x, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            2.0*(px_d2_pdsig), 2.0*(px_d2_pdpi),
            spxd2_spxd2_ddsigma,
            spxd2_spxd2_ddpi,
            spxd2_spxd2_dddelta]
  elif (wan_i=='sp3d2-2' and wan_j =='sp3d2-1'):
    return [1.0*norm_s_s, -2.0*norm_s_p*spsig_s_x, -norm_p_p*ppsig_x_x, -norm_p_p*pppi_x_x, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            -2.0*(px_d2_pdsig), -2.0*(px_d2_pdpi),
            spxd2_spxd2_ddsigma,
            spxd2_spxd2_ddpi,
            spxd2_spxd2_dddelta]
  elif (wan_i == 'sp3d2-3' and wan_j == 'sp3d2-4'): 
    return [1.0*norm_s_s, 2.0*norm_s_p*spsig_s_y, -norm_p_p*ppsig_y_y, -norm_p_p*pppi_y_y, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            2.0*(py_d2_pdsig), 2.0*(py_d2_pdpi),
            spyd2_spyd2_ddsigma,
            spyd2_spyd2_ddpi,
            spyd2_spyd2_dddelta]
  elif (wan_i=='sp3d2-4' and wan_j =='sp3d2-3'):
    return [1.0*norm_s_s, -2.0*norm_s_p*spsig_s_y, -norm_p_p*ppsig_y_y, -norm_p_p*pppi_y_y, -2.0*(norm_s_dz2*sdsig_s_z2 + norm_s_dx2y2*sdsig_s_x2my2), 
            -2.0*(py_d2_pdsig), -2.0*(py_d2_pdpi),
            spyd2_spyd2_ddsigma,
            spyd2_spyd2_ddpi,
            spyd2_spyd2_dddelta]
  elif (wan_i == 'sp3d2-5' and wan_j == 'sp3d2-6'):
    return [1.0*norm_s_s, 2.0*norm_s_p*spsig_s_z, -norm_p_p*ppsig_z_z, -norm_p_p*pppi_z_z, -2.0*(nsqrt2*nsqrt3**sdsig_s_z2), 
            2.0*(pz_d2_pdsig), 2.0*(pz_d2_pdpi),
            norm_spzd2_dz_dz*ddsig_z2_z2,
            norm_spzd2_dz_dz*ddpi_z2_z2,
            norm_spzd2_dz_dz*dddelta_z2_z2]
  elif (wan_i=='sp3d2-6' and wan_j =='sp3d2-5'):
    return [norm_s_s, -2.0*norm_s_p*spsig_s_z, -norm_p_p*ppsig_z_z, -norm_p_p*pppi_z_z, -2.0*(nsqrt2*nsqrt3**sdsig_s_z2), 
           -2.0*(pz_d2_pdsig), -2.0*(pz_d2_pdpi),
            norm_spzd2_dz_dz*ddsig_z2_z2,
            norm_spzd2_dz_dz*ddpi_z2_z2,
            norm_spzd2_dz_dz*dddelta_z2_z2]
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-3'):
    return [norm_s_s, norm_s_p*(spsig_s_x-spsig_s_y), norm_p_p*ppsig_x_y, norm_p_p*pppi_x_y, -2.0*(norm_s_dz2*sdsig_s_z2),
            px_d2_pdsig_od + py_d2_pdsig_od,
            px_d2_pdpi_od + py_d2_pdpi_od,
            spxyd2_spxyd2_ddsigma_od,
            spxyd2_spxyd2_ddpi_od,
            spxyd2_spxyd2_dddelta_od]
  elif (wan_i == 'sp3d2-3' and wan_j == 'sp3d2-1'):
    return [norm_s_s, -norm_s_p*(spsig_s_x-spsig_s_y), norm_p_p*ppsig_x_y, norm_p_p*pppi_x_y, -2.0*(norm_s_dz2*sdsig_s_z2),
            -1.0*(px_d2_pdsig_od + py_d2_pdsig_od),
            -1.0*(px_d2_pdpi_od + py_d2_pdpi_od),
            spxyd2_spxyd2_ddsigma_od,
            spxyd2_spxyd2_ddpi_od,
            spxyd2_spxyd2_dddelta_od]
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-4'):
    return [norm_s_s, norm_s_p*(spsig_s_x +  spsig_s_y), -norm_p_p*ppsig_x_y, -norm_p_p*pppi_x_y, -2.0*(norm_s_dz2*sdsig_s_z2),
            px_d2_pdsig_od - py_d2_pdsig_od,
            px_d2_pdpi_od - py_d2_pdpi_od,
            spxyd2_spxyd2_ddsigma_od,
            spxyd2_spxyd2_ddpi_od,
            spxyd2_spxyd2_dddelta_od]
  elif (wan_i == 'sp3d2-4' and wan_j == 'sp3d2-1'):
    return [norm_s_s, -norm_s_p*(spsig_s_x +  spsig_s_y), -norm_p_p*ppsig_x_y, -norm_p_p*pppi_x_y, -2.0*(norm_s_dz2*sdsig_s_z2),
            1.0*(px_d2_pdsig_od - py_d2_pdsig_od),
            1.0*(px_d2_pdpi_od - py_d2_pdpi_od),
            spxyd2_spxyd2_ddsigma_od,
            spxyd2_spxyd2_ddpi_od,
            spxyd2_spxyd2_dddelta_od]
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-5'):
    return [norm_s_s, norm_s_p*(sspsig_s_x + sspsig_s_z), norm_p_p*ppsig_x_z, norm_p_p*pppi_x_z, 
           (nsqrt3*nsqrt6-nsqrt6*nsqrt12)*sdsig_s_z2 + nsqrt2*nsqrt6*(sdsig_s_x2my2),
            pxz_dz2_pdsig,
            pxz_dz2_pdpi,
            sp3d2xz_d2_ddsig,
            sp3d2xz_d2_ddpi, 
            sp3d2xz_d2_dddelta
            ]
  elif (wan_i == 'sp3d2-5' and wan_j == 'sp3d2-1'):
    return [norm_s_s, -norm_s_p*(sspsig_s_x + sspsig_s_z), norm_p_p*ppsig_x_z, norm_p_p*pppi_x_z, 
           (nsqrt3*nsqrt6-nsqrt6*nsqrt12)*sdsig_s_z2 + nsqrt2*nsqrt6*(sdsig_s_x2my2),
            -pxz_dz2_pdsig,
            -pxz_dz2_pdpi,
            sp3d2xz_d2_ddsig,
            sp3d2xz_d2_ddpi, 
            sp3d2xz_d2_dddelta]
  elif (wan_i == 'sp3d2-1' and wan_j == 'sp3d2-6'):
    return [norm_s_s, norm_s_p*(sspsig_s_x - sspsig_s_z), -norm_p_p*ppsig_x_z, -norm_p_p*pppi_x_z, 
           (nsqrt3*nsqrt6-nsqrt6*nsqrt12)*sdsig_s_z2 + nsqrt2*nsqrt6*(sdsig_s_x2my2),
            pxz_dz2_pdsig,
            pxz_dz2_pdpi,
            sp3d2xz_d2_ddsig,
            sp3d2xz_d2_ddpi, 
            sp3d2xz_d2_dddelta
            ]

  else:
    sys.exit('invalid combo')



#  norm_p_p = (1.0/2.0)
#  norm_spzd2_dz_dz = (1.0/3.0)
#  norm_spxyd2_dz_dz = (1.0/12.0)
#  norm_spxyd2_dz_dx2y2 = nsqrt12*nsqrt3
# [ss\sigma, sp\sigma, pp\sigma, pp\pi, sd\sigma, #
#  pd\sigma, pd\pi,                               #
#  dd\sigma,                                      #
#  dd\pi,                                         #
#  dd\delta]                                      #

def load_matels(l,m,n,subspace_matrix, verbose=True):
  """
  Takes direction cosines, and 3x3 subspace matrix.
  """
  matel_dict = {}
  dxy_subspace = [('xy', 'xy'), ('yz', 'yz'), ('zx', 'zx'), ('xy', 'zx'), ('xy', 'yz'),
                  ('zx', 'xy'), ('zx', 'yz'), ('yz', 'xy'), ('yz', 'zx')]

  sp3d2_dxy = [('sp3d2-1', 'xy'), ('xy', 'sp3d2-1'), 
               ('sp3d2-2', 'xy'), ('xy', 'sp3d2-2'), 
               ('sp3d2-1', 'zx'), ('zx', 'sp3d2-1'),
               ('sp3d2-2', 'zx'), ('zx', 'sp3d2-2'),
               ('sp3d2-1', 'yz'), ('yz', 'sp3d2-1'),
               ('sp3d2-2', 'yz'), ('yz', 'sp3d2-2'),
               ('sp3d2-3', 'xy'), ('xy', 'sp3d2-3'), 
               ('sp3d2-4', 'xy'), ('xy', 'sp3d2-4'), 
               ('sp3d2-3', 'zx'), ('zx', 'sp3d2-3'),
               ('sp3d2-4', 'zx'), ('zx', 'sp3d2-4'),
               ('sp3d2-3', 'yz'), ('yz', 'sp3d2-3'),
               ('sp3d2-4', 'yz'), ('yz', 'sp3d2-4'),
               ('sp3d2-5', 'xy'), ('xy', 'sp3d2-5'), 
               ('sp3d2-5', 'zx'), ('zx', 'sp3d2-5'),
               ('sp3d2-5', 'yz'), ('yz', 'sp3d2-5'),
               ('sp3d2-6', 'xy'), ('xy', 'sp3d2-6'), 
               ('sp3d2-6', 'zx'), ('zx', 'sp3d2-6'),
               ('sp3d2-6', 'yz'), ('yz', 'sp3d2-6')]

  sp3d2_sp3d2 = [('sp3d2-1', 'sp3d2-1'), 
                 ('sp3d2-1', 'sp3d2-2'), 
                 ('sp3d2-2', 'sp3d2-1'), 
                 ('sp3d2-2', 'sp3d2-2'), 
                 ('sp3d2-3', 'sp3d2-3'), 
                 ('sp3d2-4', 'sp3d2-4'), 
                 ('sp3d2-5', 'sp3d2-5'), 
                 ('sp3d2-6', 'sp3d2-6')]

# This should line up the second nearest neighbours i.e. pdpi bonds should be right?
  #diagonal
  matel_dict[('xy', 'xy')] = (6, 6)
  matel_dict[('zx', 'zx')] = (7, 7)
  matel_dict[('yz', 'yz')] = (8, 8)
  #off diagonal
  matel_dict[('xy', 'yz')] = (6, 8)
  matel_dict[('xy', 'zx')] = (6, 7)
  matel_dict[('zx', 'xy')] = (7, 6)
  matel_dict[('zx', 'yz')] = (7, 8)
  matel_dict[('yz', 'xy')] = (8, 6)
  matel_dict[('yz', 'zx')] = (8, 7)

#These chosen to make pdsigma work
  matel_dict[('sp3d2-3', 'xy')], matel_dict[('xy', 'sp3d2-3')] = (0,6), (6,0)
  matel_dict[('sp3d2-4', 'xy')], matel_dict[('xy', 'sp3d2-4')] = (1,6), (6,1)
  matel_dict[('sp3d2-5', 'xy')], matel_dict[('xy', 'sp3d2-5')] = (2,6), (6,2)
  matel_dict[('sp3d2-6', 'xy')], matel_dict[('xy', 'sp3d2-6')] = (3,6), (6,3)
  matel_dict[('sp3d2-1', 'xy')], matel_dict[('xy', 'sp3d2-1')] = (4,6), (6,4)
  matel_dict[('sp3d2-2', 'xy')], matel_dict[('xy', 'sp3d2-2')] = (5,6), (6,5)

  matel_dict[('sp3d2-3', 'zx')], matel_dict[('zx', 'sp3d2-3')] = (0,7), (7,0)
  matel_dict[('sp3d2-4', 'zx')], matel_dict[('zx', 'sp3d2-4')] = (1,7), (7,1)
  matel_dict[('sp3d2-5', 'zx')], matel_dict[('zx', 'sp3d2-5')] = (2,7), (7,2)
  matel_dict[('sp3d2-6', 'zx')], matel_dict[('zx', 'sp3d2-6')] = (3,7), (7,3)
  matel_dict[('sp3d2-1', 'zx')], matel_dict[('zx', 'sp3d2-1')] = (4,7), (7,4)
  matel_dict[('sp3d2-2', 'zx')], matel_dict[('zx', 'sp3d2-2')] = (5,7), (7,5)

  matel_dict[('sp3d2-3', 'yz')], matel_dict[('yz', 'sp3d2-3')] = (0,8), (8,0)
  matel_dict[('sp3d2-4', 'yz')], matel_dict[('yz', 'sp3d2-4')] = (1,8), (8,1)
  matel_dict[('sp3d2-5', 'yz')], matel_dict[('yz', 'sp3d2-5')] = (2,8), (8,2)
  matel_dict[('sp3d2-6', 'yz')], matel_dict[('yz', 'sp3d2-6')] = (3,8), (8,3)
  matel_dict[('sp3d2-1', 'yz')], matel_dict[('yz', 'sp3d2-1')] = (4,8), (8,4)
  matel_dict[('sp3d2-2', 'yz')], matel_dict[('yz', 'sp3d2-2')] = (5,8), (8,5)
#
  matel_dict[('sp3d2-3', 'sp3d2-3')] = (0,0)
  matel_dict[('sp3d2-4', 'sp3d2-4')] = (1,1)
  matel_dict[('sp3d2-5', 'sp3d2-5')] = (2,2)
  matel_dict[('sp3d2-6', 'sp3d2-6')] = (3,3)

  matel_dict[('sp3d2-1', 'sp3d2-1')] = (4,4)
  matel_dict[('sp3d2-2', 'sp3d2-2')] = (5,5)
  matel_dict[('sp3d2-1', 'sp3d2-2')] = (4,5)
  matel_dict[('sp3d2-2', 'sp3d2-1')] = (4,5)

  features = []
  matels = []
  pair_list = []

  print 'Direction Cosines', l,m,n
  #for pair in sp3d2_dxy:
  full_matrix = sp3d2_diag + dxy_subspace + sp3d2_dxy
  #for pair in full_matrix:
  for pair in sp3d2_sp3d2:
      i,j = matel_dict[pair]
      feature = matel_coeffs(l,m,n, pair[0], pair[1])
      matel = subspace_matrix[i][j]
      if verbose:
        print pair
        print feature, np.round(matel,4)
      features.append(feature)
      matels.append(matel)
      pair_list.append(pair)
  return features, matels, pair_list

def gen_input(Fe_nn, R_ji, spin_mask=None):
  """
  Args:
    atom_dict(dict): Dictionary keyed off by atomic coordinates, containing Hamiltonian
    matrix elements between Wannier functions on different centers.
    nn_list: List of coordinates (as tuples) for the desired neighbour shell.
    spin_mask: If present only spin up or spin down (alternating rows) are considered.
  """
  features, labels, pair_list= [], [], []
  #for nn in nn_list:
  #Fe_nn = atom_dict[nn]/units.Ry
  if spin_mask != None:
    Fe_nn = Fe_nn[spin_mask,:]
    Fe_nn = Fe_nn[:, spin_mask]
  #subspace_matrix = Fe_nn[-3:,-3:].real
  #use full matrix
  subspace_matrix = Fe_nn[:,:].real
  l, m, n = direction_cosine(np.array(R_ji))
  features_tmp, labels_tmp, pairs_tmp = load_matels(l, m, n, subspace_matrix)
  features.extend(features_tmp)
  labels.extend(labels_tmp)
  pair_list.extend(pairs_tmp)
  return features, labels, pair_list
  
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
combined_pairs = []
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
  features, labels, pair_list =  gen_input(Fe_1NN/units.Ry, cart_int, spin_mask=spinup_mask)
  combined_feats.extend(features)
  combined_labels.extend(labels)
  combined_pairs.extend(pair_list)

fit_tb(combined_feats, combined_labels, combined_pairs)

#2nd NearestNeighbours
combined_feats = []
combined_labels = []
combined_pairs = []
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
  features, labels, pair_list =  gen_input(Fe_2NN/units.Ry, cart_int, spin_mask=spinup_mask)
  combined_feats.extend(features)
  combined_labels.extend(labels)
  combined_pairs.extend(pair_list)

fit_tb(combined_feats, combined_labels, combined_pairs)

