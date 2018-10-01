Preface
=======
This library of programs for the recursion method was written by Chris
Nex and is made available at no cost. Chris and I only ask that users
acknowledge the recursion library when it has made a significant contribution
to their work, and that they notify me of any errors they discover. I
strongly encourage users to consult with me before embarking on large projects
to avoid the pitfalls others have encountered and to avoid duplication of
effort.

Please enjoy the opportunities these programs offer.

Roger Haydock
Eugene, Oregon
2000 April 5

Project Description
--------------------
The recursion subroutine library contains routines to enable the
computation of the integrated density of states and integrals of
an arbitrary function times the density of states to be effected. 
The appropriate Hamiltonians are those for a finite (60 - 6000)
cluster of atoms with n orbitals ( 1 - 10 or more) and 'neighbour'
interactions, giving an effective large, sparse matrix eigen-
problem for which is required the distribution of eigenvalues
corresponding to the local density of states. 

The mathematics of the recursion method is described in Heine, Haydock, Kelly &
Bullett [1] and of the computation of integrals in Nex [2]. For a
description of the basic mathematical routines see Nex [3],
although there are a few additional specialised routines included
in this collection. Two independant parameters at the 
disposal of the user, giving a control on the accuracy, 
are the size of the cluster and the 'number of levels' in the recursion. 
The appropriate values of these parameters will vary 
from problem to problem depending on
the Hamiltonian, the nature of the lattice, and of the physical
quantity of interest ultimately being computed, and should be
determined by some trial runs of the program. For five d-electrons
per site in fcc and bcc clusters and the calculation of structural
energies it has been found that the number of levels appropriate
is 2-5 times the 'radius' of the cluster (in numbers of atoms) in
cuboid clusters. A mathematically exact error bound on the
integrated density of states is available if only exact
coefficients are used, but, as is the nature of such quantities,
this bound is pessimistic, usually by at least an order of
magnitude. 

References
------------
1. Heine, Haydock, Kelly & Bullett Solid State Physics Vol. 35 (1979)
2. Nex C.M.M. J. Phys. A Vol. 11, 653 et seq. (1978)
3. Nex C.M.M. Comp. Phys. Comm. (1980) 
4. `Nex C.M.M. Comp. Phys. Comm. (1984) <https://doi.org/10.1016/0010-4655(84)90163-2>`_
