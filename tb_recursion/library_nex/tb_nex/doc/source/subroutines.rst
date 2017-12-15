Specification of Subroutines
=============================

.. f:subroutine:: RECAL(HOP,PSI,PMN,M,A,B2,LL)

  DIMENSION PSI(M), PMN(M), A(LL), B2(LL)
  Computes the tri-diagonalization of a large sparse symmetric matrix
  using the recursion (or Lanczos or Paige's) method:

  .. math::

    A_N = <\psi_{N},H \psi_{N}> \\
    B_{N+1}\psi_{N+1} = (H-A_{N})\psi_{N} - B_{N}\psi_{N-1}\\
    <\psi_{N+1},\psi_{N+1}> = 1

  This routine must be accompanied by a routine segment hop to carry
  out the steps of this process.


  
