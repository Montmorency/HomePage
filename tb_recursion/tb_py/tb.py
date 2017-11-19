import numpy as np

class DiamondElements():
  def __init__(k=np.array([0,0,0]), two_center=False, alat):
    """
    Returns Slater-Koster components of Hamiltonian for diamond structure.
    Variable names chosen to match as closely as possibly those of the paper.
    eta = alat*k_x
    """
    self.xi = alat*k[0]
    self.eta = alat*k[1]
    self.zet = alat*k[2]
#these are the parameters we need to fit.
    self.integrals = ['E_ss_000','E_ss_110','E_xx_000', 'E_xx_110', 'E_xx_011','E_ss_hhh','E_sx_hhh','E_sx_011','E_sx_110','E_xx_hhh','E_xy_hhh','E_xy_110', 'E_xy_011']

  def matrix_elements(self):
    """
    Returns the tight binding hamiltonian matrix elements for first and second neighbours for a given $\k$ point.
    """
    ss_11 = (self.integrals['E_ss_000'] + 4.0*self.integrals['E_ss_110']*(np.cos(self.xi)*np.cos(self.eta) +
                                                                          np.cos(self.eta)*np.cos(self.zeta) +
                                                                          np.cos(self.xi)*np.cos(self.zeta)))
    ss_22 = ss_11

    xx_11 = (self.integrals['E_xx_000'] + 4.0*self.integrals['E_xx_110']*(np.cos(self.xi)*np.cos(self.eta) + np.cos(self.xi)*np.cos(self.zeta))
                                        + 4.0*self.integrals['E_xx_011']*(np.cos(self.eta)*np.cos(self.zeta)))
    xx_22 = xx_11

    ss_12 = 4.0*self.integrals['E_ss_hhh']*((np.cos(0.5*self.xi)*np.cos(0.5*self.eta)*np.cos(0.5*self.zeta)) - 
                                             i*(np.sin(0.5*self.xi)*np.sin(0.5*self.eta)*np.sin(0.5*self.zeta)))

    sx_12 = 4.0*self.integrals['E_sx_hhh']*(i*(np.sin(0.5*self.xi)*np.cos(0.5*self.eta)*np.cos(0.5*self.zeta))-np.cos(0.5*self.xi)*np.sin(0.5*self.eta)*np.sin(0.5*self.zeta))
    sx_21 = -1.0*np.conjg(sx_12)

    sx_11 = -4.0*self.integrals['E_sx_011']*np.sin(self.xi)*np.sin(self.eta) + (4.0*i*self.integrals['E_sx_110']*(np.sin(self.xi)*np.cos(self.eta)
                                                                             +  np.sin(self.xi)*np.cos(self.zeta)))
    sx_22 = -1.0*np.conjg(sx_11)

    xx_12 = 4.0*self.integrals['E_xx_hhh']*((np.cos(0.5*self.xi)*np.cos(0.5*self.eta)*np.cos(0.5*self.zeta)) -
                                            i*(np.sin(0.5*self.xi)*np.sin(0.5*self.eta)*np.sin(0.5*self.zeta)))
    xx_21 = np.conjg(xx_12)
    xy_12 =  4.0*self.integrals['E_xy_hhh']*((i*np.cos(0.5*self.xi)*np.cos(0.5*self.eta)*np.sin(0.5*self.zeta)) -
                                                (np.sin(0.5*self.xi)*np.sin(0.5*self.eta)*np.cos(0.5*self.zeta)))
    xy_21 = np.conjg(xy_12)
    yx_12 = xy_12

    xy_11 = -4.0*self.integrals['E_xy_110']*np.sin(self.xi)*np.sin(self.eta) - (4.0*i*self.integrals['E_xy_011']*(np.sin(self.xi)*np.cos(self.zeta)
                                                                             -  np.sin(self.eta)*np.cos(self.zeta)))
    xy_22 = np.conjg(xy_11)
    return []

if __name__== "__main__":
  #take target data as a set of k points and eigenvalues.
  target_data = []
  d_el_k = DiamondElements()
  fit_integrals(d_el_k, DiamondElements())





