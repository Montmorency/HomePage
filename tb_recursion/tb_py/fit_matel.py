import json
import numpy as np
#import tensorflow as tf
from pymc import deterministic, Normal
import matplotlib.pyplot as plt

class XC(object):
  """
  Methods that, given a local density return ec, vc, ex, vx.
  Correlation and exchange returned separately. Arbitrary functions 
  can be incorporated but pz, pw91, and pbe should
  be more than sufficient.
  """
  def __init__(self):
    pass

  def slater(self, rho):
    pi34 = 0.6203504908994

    if rho < 1.e-10: 
      return 0.0,0.0

    rs = (pi34/rho)**(1./3.)

    if rho < 0.0:
      sys.exit("Negative Density!")
    #f=-9/8*(3/2*pi)**(2/3)
    f = -0.687247939924714
    alpha = 2./3.
    ex = f*alpha/rs
    vx = 4.0/3.0*f*alpha/rs
    return ex, vx

  def pz(self, rho) :
    """
    Perdew-Zunger parameterization of the exchange correlation energy.
    Taken from functionals.f90 quantum-espresso org.
    """
    pi34 = 0.6203504908994

    a = 0.0311
    b = -0.048
    c = 0.002
    d = -0.0116
    b1 = 1.0529
    b2 = 0.3334
  
#if densities very small just get out:
    if rho < 1.e-10: 
      return 0.0,0.0
    else:
      rs = (pi34/rho)**(1./3.)
    if (rs < 1.0):#high density  limit
      ec = a*np.log(rs) + b + c*np.log(rs)*rs + d*rs
      vc = a*np.log(rs) + (b-a/3.0) + (2./3.)*c*rs*np.log(rs) + ((2.0*d-c)/3.0)*rs
    else:#interpolation formula
      ox = 1.0 + b1*np.sqrt(rs) + b2*rs 
      dox = 1.0 + (7.0/6.0)*b1*np.sqrt(rs) + (4.0/3.0)*b2*rs
      ec = (-0.1423)/ox
      vc = ec*dox/ox
    return [ec, vc]

class Evaluate(object):
  def integrate_wivxcwj(vxc, w_i, w_j):
#call fortran to evaluate and return matrix element with local density.
    vxc_matel(vxc, w_i, w_j)
    return vxc_matel

class VxcPerceptron(object):
  def __init__(self, target_data=[], dim=[10]):
    self.targets = target_data
    self.N = 10

  def model(self):
    #the Pade matrix has powers of the density in it.
    W = tf.placeholder(tf.float64, [None, self.N])
    y = evaluate.calc_vxc(coeffs, self.target_data)
    cost = tf.reduce_mean(tf.square(self.targets-y))
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

  def train(self):
    init = tf.flobal_variables_initializer()
    with tf.Session() as sess:
      sess.run(init)
    self.vxc_model.train(steps=steps)
    return self.vxc_model

#class VxcPade(object):
#  """
#  Fits coefficients of the Pade approximant using Bayesian regression,
#  original values of the coefficients chosen to be normally distributed.
#  """
#  def __init__(self):
#    N_xc = 10
#    coefficients = [Normal('c{}'.format(i), mu=1.0, tau=1.0) for i in range(N_xc)]
##load "observed" vxc matrix elements:
#  with open('exc_elements.json') as f:
#    exc_elements = json.load(f)
#  with Model() as fit_vxc: # model specifications in PyMC3 are wrapped in a with-statement
#    start = find_MAP() # Find starting value by optimization
#    step = NUTS(state=start) # Instantiate MCMC sampling algorithm
#    trace = sample(2000, step, start=start, progressbar=False) 

def plot_vxc_pz():
  densities = np.arange(0,3,0.1)
  xc = XC()
  ec = [xc.pz(rho)[0]  for rho in densities]
  vc = [xc.pz(rho)[1]  for rho in densities]
  ex = [xc.slater(rho)[0]  for rho in densities]
  vx = [xc.slater(rho)[1]  for rho in densities]
  fig = plt.figure()
  ax = fig.add_subplot(211)
  ax.plot(densities, ec)
  ax.plot(densities, vc)

  ax2 = fig.add_subplot(212)
  ax2.plot(densities, ex)
  ax2.plot(densities, vx)
  plt.show()

if __name__ == "__main__":
#plot densities
  plot_vxc_pz()


