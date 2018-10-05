Library Overview
=================

Initialisation Routines
-------------------------

The routines in the library fall into two categories: those
which assist in the setting up of the problem, and those which
compute from the calculated coefficients the local density of
states and related functions. These are described separately in
the following two sections. There exists in the literature some slight confusion of notation,
but that used in the library is consistent within itself and with
the recursion formula:

.. math::

  b_{n+1}\psi_{n+1} = [H - a_{n}]\psi_{n} - b_{n}\psi_{n-1}, n=0,1,.
  
Where :math:`a_{n} = \langle\psi_{n}|H \psi_{n}\rangle` and :math:`b_{n+1}` is 
chosen such that :math:`\langle\psi(n+1)|\psi(n+1)\rangle = 1`, :math:`\psi_{-1} =0`, and 
:math:`\psi_{0}` are given. 

(Within the Fortran the subscripts are increased by 1 to conform to
Fortran IV convention, although all routines have also run in
Fortran 77 under the IBM vs compiler)

The system adopted for the setting up of the problem for the
application of the recursion method consists of a library of
routines divided into three sections, with only a subset being
used in any one problem, and for some problems the user supplying
their own routines to replace and supplement some of the library
ones. 

The three groups of routines are as follows: the first serve
to define a cluster of 'atom sites', the next set up a 'neighbour
map' (or matrix sparsity pattern) , and the last set up the
component submatrices of the Hamiltonian. The appropriate results
of these routines are loaded by the user into a common block which
should then contain all the information necessary to define the
matrix representation of the Hamiltonian. This is then used by the
subroutine at the heart of the recursion library ( :f:subr:`RECAL` ) to
operate  with the  Hamiltonian  and  produce  the recursion
coefficients. 

Subroutines for Generating a Cluster 
---------------------------------------------

:f:subr:`BCCLAT` and :f:subr:`FCCLAT` generate cuboid 
clusters in bcc and fcc lattices respectively. More exotic 
cluster shapes may be obtained by use  of the  following  
routines  but their  use needs justification. Onion will 
classify each atom according to the
number of 'hops' from a given group of 'central' atoms , but can
only be used after the atom 'neighbours' have been defined 
(e.g. by :f:subr:`NNCAL`); :f:subr:`PEEL` excises atoms from a cluster and
renumbers the remaining ones. The last two routines may be useful
to users with particular cluster-shape requirements, but in
general the cuboid clusters have proved satisfactory for almost
all applications to date. 

Setting-up of the Hamiltonian matrix. 
--------------------------------------

:f:subr:`NNCAL` is the central subroutine used to set up a 'neighbour map'
of the atom cluster. This uses another routine (:f:subr:`IBONDS`) defining
the neighbour criterion; this may be replaced by the user with
more particular requirements (i. e. first and second nearest
neighbours) or a different lattice.

Once the neighbour map is set up, the non-equivalent interactions 
may be set up using :f:subr:`MMCAL`. The 'equivalence' of two interactions 
is defined in a routine :f:subr:`EQUIV` which again may be replaced by a 
user with specific requirements. :f:subr:`MMCAL` generates a list of 
non-equivalent vectors which may then be used to 
set up the 'interaction blocks' of the
Hamiltonian matrix, for instance as in :f:subr:`SLKODE` to give the Slater-
Koster interaction matrices for the 5 d-orbitals. Subroutines
generating the s, p or d Slater-Koster matrices are supplied in
the library and have names starting SK--. The essential information 
defining the Hamiltonian matrix for a regular lattice may be 
set up by the routine :f:subr:`SETUP` and amendments
to the matrix made using :f:subr:`ADDAT`, :f:subr:`ORPEEL`, 
or :f:subr:`PEEL` (adding or excising atoms from the cluster ). 

For a random system an alternative form of storage of the Hamiltonian 
matrix is probably more efficient and routines :f:subr:`SCAN` 
and :f:subr:`SCAN1` are supplied to facilitate the 
setting up and operation of the matrix of such a problem. 

Generation of Recursion Coefficients 
-------------------------------------

The appropriate results of the above subroutines should then be
loaded by the user into a common block and used by a routine,
such as :f:subr:`HOP` of the example run, to supply the Hamiltonian
operator required by the actual recursion subroutine, :f:subr:`RECAL`. 

Processing the Recursion Coefficients 
----------------------------------------

Two methods are provided for the  calculation of the local
density of states : terminor and quadrature, as described in Nex `[3]`_. 
If integrals are the final object of the calculation then quadrature
(:f:subr:`DENQD`) is appropriate , while to estimate the density itself with
secure knowledge of the band-gaps the analytic terminator 
(:f:subr:`DENSQ` or :f:subr:`DENCRS`) can be used.

If the sum of several densities of states is required 
(e.g. over several orbitals), this again may be done efficiently by 'summing'
their tridiagonalisations using :f:subr:`RECSUM` to produce a resultant set
of coefficients. If the functions summed are very different in
character, this may introduce significant error, in which case the
functions should be tabulated separately , but again for graphical
purposes the results of :f:subr:`RECSUM` are usually adequate.
For quantitative work , only the integrated density of states,
:math:`N(E)`, function values should be relied on, and any other functions
computed from this one. It should be noted that the approximations
obtained from the quadrature formula are not analytically related
in the way one might expect:

.. math::
  :label: intdos

  \int E N(E) dE \neq  E N(E) - \int N(E)dE

where the L.H.S of Eq. :eq:`intdos` and N(E) represent values obtained directly from
:f:subr:`DENQD`. Such identities are usually satisfied approximately, but if
precise analyticity is demanded the approximation to :math:`N(E)` should
be taken and all other results computed directly from it. The
exception is that :math:`\frac{d}{dE}\int N(E)dE = N(E)` by definition. 

