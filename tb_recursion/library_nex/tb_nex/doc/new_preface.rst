Preface
=======

The recursion subroutine library contains routines to enable the
computation of the  integrated  density of states and integrals of
an arbitrary  function times the density of states to be effected. 
The  appropriate  Hamiltonians  are those for a finite (60 - 6000)
cluster of atoms with n orbitals ( 1 - 10 or more) and 'neighbour'
interactions , giving  an effective  large , sparse  matrix eigen-
problem for which is  required  the  distribution  of  eigenvalues
corresponding to the  local density of states. The  mathematics of
the  recursion  method  is  described  in  Heine, Haydock, kelly &
Bullett (1) and of the computation of integrals in Nex (2). For a
description  of  the  basic  mathematical  routines  see  Nex (3),
although there are a few additional specialised  routines included
in this collection. Two independant parameters at the disposal of the user, giving a
control on the  accuracy , are the  size of the  cluster  and  the
'number of levels' in  the  recursion. The  appropriate  values of
these  parameters will vary from  problem to problem  depending on
the Hamiltonian, the  nature of the  lattice, and of the  physical
quantity  of  interest  ultimately  being  computed, and should be
determined by some trial runs of the program. For five d-electrons
per site in fcc and bcc clusters and the calculation of structural
energies it has been  found that the number of levels  appropriate
is  2-5 times the 'radius' of the cluster (in numbers of atoms) in
cuboid  clusters. A  mathematically  exact  error  bound  on  the
integrated   density  of  states   is  available   if  only  exact
coefficients are used , but , as is the nature of such quantities,
this  bound is  pessimistic , usually  by  at  least  an  order of
magnitude. 

References
==========

1. HEINE,HAYDOCK,KELLY & BULLETT  SOLID STATE PHYSICS VOL.35 (1979)
2. NEX C.M.M. J.PHYS. A  VOL.11 653 ET SEQ. (1978)
3. NEX C.M.M. COMP.PHYS.COMM  (1985)
