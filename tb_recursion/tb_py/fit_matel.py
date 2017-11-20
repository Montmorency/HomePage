import json
import numpy as np
import tensorflow as tf
from pymc import deterministic, Normal


class Vxc(object):
  

class Evaluate(object):
  def integrate_wivxcwj(vxc, w_i, w_j):
#call fortran to evaluate and return matrix element with local density.
    vxc_matel(vxc, w_i, w_j)
    return vxc_matel

class FitVxcNN(object):
  def __init__(self, target_data=[], dim=[10]):
    self.targets = target_data
    self.N = 10

  def model(self):
    #the Pade matrix has powers of the density in it.
    W = tf.placeholder(tf.float64, [None, self.N])
    y = evaluate.calc_vxc(coeffs, self.target_data)
    cost = tf.reduce_mean(tf.square(self.targets-y))
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

  def train(self, N=100) 
    init = tf.flobal_variables_initializer()
    with tf.Session as sess():
      sess.run(init)
    self.vxc_model.train(steps=steps)
    return self.vxc_model

class FitVxcPade(object):
  """
  Fits coefficients of the Pade approximant using Bayesian regression.
  """
  def __init__(self):
    N_xc = 10
    coefficients = [Normal('c{}'.format(i), mu=1.0, tau=1.0) for i in range N_xc]

#load "observed" vxc matrix elements:
  with open('exc_elements.json') as f:
    exc_elements = json.load(f)

  observed_v_xc = Gaussian('observed_vxc', mu= value=exc_matels, observed=True)

  with Model() as fit_vxc: # model specifications in PyMC3 are wrapped in a with-statement
    start = find_MAP() # Find starting value by optimization
    step = NUTS(state=start) # Instantiate MCMC sampling algorithm
    trace = sample(2000, step, start=start, progressbar=False) 


if __name__ == "__main__":

rs = np.arange()






