Specification of Subroutines
=============================

.. f:subroutine:: RECAL(HOP,PSI,PMN,M,A,B2,LL)

  Computes the tri-diagonalization of a large sparse symmetric matrix
  using the recursion (or Lanczos or Paige's) method:

  DIMENSION PSI(M), PMN(M), A(LL), B2(LL)

  .. math::

    a_n = <\psi_{n},H \psi_{n}> \\
    b_{n+1}\psi_{n+1} = (H-a_{n})\psi_{n} - b_{n}\psi_{n-1}\\
    <\psi_{n+1},\psi_{n+1}> = 1


  This routine must be accompanied by a routine segment hop to carry
  out the steps of this process.

.. f:subroutine:: RECSUM(AC,BC,NA,LL,NP,A,B2,EPS,WK,NW)

	DIMENSION AC(NA,NP),BC(NA,NP),WK(NW),A(NA),B(NA)
  AC    AC(N,M) = A(N-1,M),  N=1,M,  M=1,NP
  BC    BC(N,M) = B(N-1,M)**2,  N=1,LL,  M=1,NP
  NA    FIRST DIMENSION OF ARRAYS AC AND BC IN CALLING PROGRAM
  LL*   ON INPUT : THE ABSOLUTE VALUE GIVES LENGTH +1 OF EACH
               TRIDIAGONALISATION. IF >0 M=LL-1 ;  IF <0 M=LL
        ON OUTPUT: LENGTH OF OUTPUT TRIDIAGONALISATION ,
          IF NEGATIVE THEN RECQD FAILED WITH TOO FEW ROOTS
  NP    NUMBER OF CONTINUED FRACTIONS
  A*    A(I) = A(I-1), I=1,M, IN TRIDIAGONALISATION CORRESPONDING TO W(X)
  B2*   B2(I) = B(I-1)**2, I=1,LL, IN TRIDIAGONALISATION CORRESPONDING TO W(X)
  EPS   ACCURACY REQUIRED IN COMPUTATION
  WK*   REAL WORK ARRAY OF LENGTH AT LEAST 5*LL*NP
	NW    LENGTH OF ARRAY WK

	Computes the tridiagonalisation (continued fraction, Jacobi matrix)
	corresponding to the sum of NP tridiagonalisations, :math:`w_{m}(x)`.
	.. math::

     \sqrt{b_{n+1,m}} P_{n+1,m}(x) = (x-a_{n,m})P_{n,m}(x)-\sqrt(b_{n,m})P_{n-1,m}(x)
     w(x) = \sum_{m=1}^{NC} B(0,m) w_{m}(x)

	Note that this routine uses RECQD, CFGEN, RECRTS, NUMC, NUMD.

.. f:function:: FUNCTION RECWT(E,A,B2,LL,EPS,N,P,NS)

	ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  DIMENSION A(LL),B2(LL),P(2,3)
  E    REQUIRED FIXED POINT IN QUADRATURE. IT MAY BE A NODE OF
       THE LL-1 OR LL QUADRATURE IF A(LL) IS APPROPRIATELY DEFINED

  A*   DIAGONAL ELEMENTS OF THE RECURRENCE. IF N IS CHANGED
       FROM -1 INPUT TO 0 ON OUTPUT THEN A(LL) CONTAINS THE ADJUSTED
       VALUE TO ACHIEVE A GAUSSIAN NODE AT E, OTHERWISE A IS
       UNCHANGED.

  B2   SUB-DIAGONAL ELEMENTS OF THE RECURRENCE

  LL   INDEX OF LAST B2 VALUE TO BE USED

  EPS  RELATIVE THRESHOLD VALUE OF THE POLYNOMIAL BELOW WHICH E WILL BE ACCEPTED AS A ZERO

  N*   CODE :

       -1   A(LL) TO BE OVERWRITTEN. N CHANGED TO 0 IF SUCCESSFUL, UNCHANGED OTHERWISE
        0   A(LL) GIVEN (NOT OVERWRITTEN)
        1   A(LL) NOT COMPUTED EXPLICITLY (NOT OVERWRITTEN)

  P* FINAL POLYNOMIAL VALUES USED IN CALCULATION OF WEIGHT TO BE USED UNCHANGED IF ROUTINE IS RE-ENTERED WITH NS=LL

       IF N=LL-IABS(N)
            P(2,1)=P(E,N)       P(1,1)=P(E,N-1)
            P(2,2)=P'(E,N)      P(1,2)=P'(E,N-1)
            P(2,3)=Q(E,N-1)     P(1,3)=Q(E,N-2)
       Q(E,M) IS THE POLYNOMIAL OF THE SECOND KIND OF DEGREE M

  NS   POINT AT WHICH RECURRENCE IS INITIATED . THIS SHOULD BE
       1 INITIALLY , BUT FOR A SUBSEQUENT CALL, WITH E UNCHANGED AND LARGER LL, SHOULD BE SET TO THE CURRENT VALUE OF LL

  Computes the value of the weight at the fixed point in a 1-fixed
  point Gaussian quadrature, given the corresponding 3-term recurrence
  relation:

  .. math::

    P(E,J)= (E-A(J)) * P(E,J-1) - B2(J) * P(E,J-2)

  This routine may be called repeatedly with increasing number
  of levels such that it does not recompute earlier polynomial
  values. If required the value of the last coefficient, A(LL),
  may be computer, or it may be assumed that this has already been
  done and that value used in the calculation of the weight.
  The expression for the weight used is (with N=LL).

  .. math::

    \frac{P(E,N-1)*Q(E,N-1)-P(E,N)*Q(E,N-2)}{P(E,N-1)*P'(E,N)-P'(E,N-1)*P(E,N)+P(E,N)**2/B2(N)}

  As this form is independent of the normalisation of the polynomials. P and Q are the monic
  polynomials of the first and second kinds.

