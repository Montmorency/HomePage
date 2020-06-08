Specification of Subroutines
=============================

.. f:subroutine:: RECAL(HOP,PSI,PMN,M,A,B2,LL)

Computes the tri-diagonalization of a large sparse symmetric matrix
using the recursion (or Lanczos or Paige's) method:

.. math::
  a_n = <\psi_{n},H \psi_{n}> \\
  b_{n+1}\psi_{n+1} = (H-a_{n})\psi_{n} - b_{n}\psi_{n-1}\\
  <\psi_{n+1},\psi_{n+1}> = 1

This routine must be accompanied by a routine segment hop to carry
out the steps of this process:

::

      SUBROUTINE HOP(X,Y,A)
      DIMENSION X( ),Y( )

which, with input :math:`X=X`, :math:`Y=Y`, will produce :math:`Y = H.X+Y` AND 
:math:`A = X.H.X = \langle X,HX\rangle`.

::

  ARGUMENTS : (* DENOTES OVERWRITTEN ARGUMENT)
  HOP   NAME OF SUBROUTINE SEGMENT
  PSI*  INPUT B(0)*PSI(0)  OUTPUT B(LL)*PSI(LL)
  PMN*  INPUT -B(0)*PSI(-1)  OUTPUT -B(LL)*PSI(LL-1) WHERE
        THE NORMS OF ALL VECTORS, PSI(I), ARE EQUAL TO UNITY
  M     DIMENSION OF MATRIX
  A*    OUTPUT A(I) = A(I-1), I=1,LL-1
  B2*   INPUT  B2(1) = B(0)**2
        OUTPUT B2(I) = B(I-1)**2, I=1,LL
  LL    1+LENGTH OF TRIDIAGONALISATION REQUIRED > OR = 2

.. f:subroutine:: HOP(X,Y,A)

This subroutine segment hop, will carry out steps of this process,
which with input X=X, Y=Y will  produce Y = HX + Y and 
A = X.H.X. = :math:`\langle X,HX \rangle`.


.. f:subroutine:: PLYVAL(E,A,B2,LM1,P,Q)

::

  DIMENSION A(LM1),B2(LM1),P(2),Q(2)

Computes relative values of last two polynomials of each
kind of an orthogonal sequence defined by a three-term recurrence
relation.
    
with :math:`P(x,-1)=Q(X,-1) = 0.0`, :math:`P(X,0)=1.0`, and 
:math:`Q(X,1)=B2(1)`

.. math::
    P(X,I) = (X-A(I))P(X,I-1)-B2(I)P(X,I-2)\\
    Q(X,I-1)  = (X-A(I))Q(X,I-2) - B2(I)Q(X,I-3)

The same normalization is used for P and Q, but its value varies
with X in order to maintain accuracy in the relative values of 
the functions.

::

  E    ARGUMENT OF POLYNOMIALS
  A    DIAGONAL RECURRENCE COEFFICIENTS I=1,LM1
  B2   OFF-DIAGONAL RECURRENCE COEFFICIENTS I=1,LM1
  LM1  MAXIMUM DEGREE OF POLYNOMIALS
  P*   P(I) CONTAINS VALUES OF THE POLYNOMIALS OF THE FIRST
       KIND P(E,LM1+I-2) ,I=1,2 , WITH ARBITRARY NORMALISATION
  Q*   Q(I) CONTAINS VALUES OF THE POLYNOMIALS OF THE SECOND
        KIND Q(E,LM1+I-3) ,I=1,2 , WITH ARBITRARY NORMALISATION


.. f:subroutine:: RECQD(A,B2,LM1,X,WT,M,EPS,WK,MU)

::

    DIMENSION A(LM1),B2(LM1),X(LM1),WT(LM1),WK(LM1),MU(LM1)

Calculates the coefficients of the Gaussian quadrature corresponding to a
given J-matrix or three-term recurrence relation. The Sturm sequence property
is used to locate the nodes of the quadrature and the expression given
in :f:subr:`RECWT` to calculate the weights. (Cheney)

This routine calls, :f:subr:`RECWT`, :f:subr:`RECRTS`, :f:subr:`NUMC` &  
:f:subr:`NUMD`.

::

  A    DIAGONAL ELEMENTS OF THE J-MATRIX I=1,LM1
  B2   SQUARES OF THE OFF-DIAGONAL ELEMENTS OF THE J-MATRIX I=2,LM1
         B2(1) GIVES THE NORMALISATION OF THE QUADRATURE
  LM1  DIMENSION OF THE J-MATRIX
  X*   NODES OF GAUSSIAN QUADRATURE
  WT*  WEIGHTS OF GAUSSIAN QUADRATURE
  M*   NUMBER OF QUADRATURE NODES . IF DIFFERENT FROM LM1 THEN
       THE ROUTINE HAS FAULTED.
  EPS  ACCURACY REQUIRED IN NODE CALCULATION
  WK*  WORK ARRAY OF LENGTH AT LEAST LM1
  MU*  WORK ARRAY OF LENGTH AT LEAST LM1 (O/P FROM RECRTS)


.. f:subroutine:: RECSUM(AC,BC,NA,LL,NP,A,B2,EPS,WK,NW)

Computes the tridiagonalisation (continued fraction, Jacobi matrix)
corresponding to the sum of NP tridiagonalisations :math:`w_{m}(x)`.

.. math::
  \sqrt{b_{n+1,m}} P_{n+1,m}(x) = (x-a_{n,m})P_{n,m}(x)-\sqrt{b_{n,m}}P_{n-1,m}(x)\\
  w(x) = \sum_{m=1}^{NP} b(0,m) w_{m}(x)

::

  DIMENSION AC(NA,NP),BC(NA,NP),WK(NW),A(NA),B(NA)
  AC    AC(N,M) = A(N-1,M),  N=1,M,  M=1,NP
  BC    BC(N,M) = B(N-1,M)**2,  N=1,LL,  M=1,NP
  NA    FIRST DIMENSION OF ARRAYS AC AND BC IN CALLING PROGRAM
  LL*   ON INPUT : THE ABSOLUTE VALUE GIVES LENGTH +1 OF EACH
        TRIDIAGONALISATION. IF >0 M=LL-1 ;  IF <0 M=LL
        ON OUTPUT: LENGTH OF OUTPUT TRIDIAGONALISATION ,
        IF NEGATIVE THEN RECQD FAILED WITH TOO FEW ROOTS
  NP    NUMBER OF CONTINUED FRACTIONS
  A*    A(I) = A(I-1), I=1,M, IN TRIDIAGONALISATION 
        CORRESPONDING TO W(X)
  B2*   B2(I) = B(I-1)**2, I=1,LL, IN TRIDIAGONALISATION 
        CORRESPONDING TO W(X)
  EPS   ACCURACY REQUIRED IN COMPUTATION
  WK*   REAL WORK ARRAY OF LENGTH AT LEAST 5*LL*NP
  NW    LENGTH OF ARRAY WK

Note that this routine uses :f:subr:`RECQD`, 
:f:subr:`CFGEN`, :f:subr:`RECRTS`, 
:f:subr:`NUMC`, :f:subr:`NUMD`.

.. f:subroutine:: RECNO(HOP,SOP,U,M,NIT,LS,LL,A,B2,PSI,PMN,EMACH)

::

  DIMENSION U(M),PSI(M),PMN(M),A(LL),B2(LL)
  DOUBLE PRECISION SUM
  COMMON /BLKNNN/SUM

Implements the 'non-orthogonal basis' recursion method
where the eigen-problem takes the form

.. math::
  MX = ESX

S is assumed to have unit diagonal elements and both M and S
to be real symmetric. The inverse of S times a vector is
estimated by taking 'NIT' applications of Gauss=Seidel iteration. If NIT
is set to be zero the routine performs the usual recursion
assuming S is the identity. The number of vectors needed is kept to a minimum
(3) by asking the user to write the S multiplication routine in 
a particular way, so PLEASE NOTE THE SPECIFICATIONS CAREFULLY. 
If S is the identity then NIT should be set to zero and the vectors
U and PSI set to refer to the same real array. The Greenian calculated may
be thought of as 

.. math::
  U^{T}S(SE-M)S^{-1}U,

where U is the starting vector. The vectors
PSI and PMN are in fact S times the usual recursion basis
vectors. The vector SU may be specified rather than U and
this indicated by replacing NIT by -NIT and storing
that vector in PSI on the first call to :f:subr:`RECNO`.

The routine may be restarted to extend the number of levels
as U,PSI and PMN and the NORM sum are all retained (the latter
in the common block /BLKNNN/.

The actual algorithm is as follows:

::

  PMN(0) = 0
  PSI(0) = S U(0)
  B2(1) = U(0)(TRANSPOSED) PSI(0)
  UP TO NUMBER OF LEVELS (LL) DO :
  A(N)   =  U(N)TRANSPOSED M U(N) / U(N)TRANSPOSED PSI(N)
  PSI(N+1) = M U(N) -A(N) PSI(N) - B2(N) PSI(N-1)
  U(N+1) =  S(INVERSE) PSI(N+1)

Here :math:`S^{-1}PSI(N+1)` is calculated by NIT applications
of the formula (I is the iteration number):

.. math::
  U(N+1)(I) = PSI(N+1) - LU(N+1)(I) - RU(N+1)(I-1)

where L and R are the strict left and right triangular parts
of S

with a renormalisation of the vectors for numerical stability.

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  HOP    NAME OF A SUBROUTINE TO PERFORM THE MATRIX MULTIPLICATION
         BY M . ITS SPECIFICATIONS MUST BE AS FOLLOWS:
              SUBROUTINE HOP(X,Y,A)
              DIMENSION X( ),Y( )
        WITH INPUT X=X, Y=Y, WILL PRODUCES Y = M.X+Y AND
        A = X.M.X = <X,MX>
  SOP   NAME OF A SUBROUTINE TO EVALUATE THE PRODUCT OF THE OFF-
        DIAGONAL ELEMENTS OF S WITH A VECTOR. THE SPECIFICATION
        IS AS FOLLOWS

            SUBROUTINE SOP(U,V,W)
            DIMENSION U( ),V( ),W( )
      CALCULATES W = V - OFF-DIAGONAL(S) U

      THIS IS CALLED WITH U AND W REFERRING TO THE SAME ARRAY TO
      ACHIEVE A GAUSS-SEIDEL STEP AND WITH V AND W REFERRING TO THE SAME
      ARRAY TO PERFORM A MORE USUAL MATRIX MULTIPLICATION.
      N.B. NOTE THE MINUS SIGN

      THE IMPLICATION FOR THE USER IS THAT THE ELEMENTS OF THE
      PRODUCT MUST BE EVALUATED AND OVERWRITTEN IN INCREASING
      ORDER, NOT BY A GLOBAL ACCUMULATION TECHNIQUE.
  U*  INPUT AS STARTING STATE IF NIT >OR= 0
        OUTPUT AS LAST RECURSION BASIS VECTOR
  M     DIMENSION OF MATRIX
  NIT   THE ABSOLUTE VALUE DEFINES THE NUMBER OF ITERATIONS IN
        THE GAUSS-SEIDEL PROCESS.  THIS IS  RATHER BETTER THAN
        NIT TERMS IN THE TAYLOR SERIES FOR  THE INVERSE OF S.
        IF NIT IS ZERO S IS ASSUMED TO BE THE IDENTITY MATRIX
        AND U AND PSI MAY REFER TO THE SAME PHYSICAL ARRAY.
        IF NIT IS NEGATIVE IT IS ASSUMED THAT PSI DEFINES THE
        STARTING STATE, AND U WILL BE CALCULATED BY RECNO USING
        GAUSS-SEIDEL ITERATION TO ESTIMATE S INVERSE PSI.
  LS    STARTING LEVEL IN THE RECURSION PROCESS. IF LS=1 THE
        NORMALISATION B2(1) IS CALCULATED AS U(TRANS)S U.
        IF GREATER THAN 1, IT IS ASSUMED THAT ALL ARGUMENTS ARE
        UNCHANGED FROM A PREVIOUS CALL TO RECNO
  LL*   ON INPUT LAST REQUIRED LEVEL IN TRIDIAGONALISATION
        ON OUTPUT, IF NEGATIVE, -LL IS THE NUMBER OF LEVELS
        REACHED BEFORE B2(LL)<EMACH
  A*    DIAGONAL ELEMENTS OF TRIDIAGONALISATION
  B2*   SQUARES OF OFF-DIAGONAL ELEMENTS OF TRIDIAGONALISATION
  PSI*  S TIMES LLTH RECURSION BASIS VECTOR
        IF NIT<0 ON INPUT THIS MUST BE INITIALISED IN THE CALLING
        ROUTINE
        IF NIT=0 THIS MUST REFER TO THE SAME ADDRESS AS U ABOVE
  PMN*  S TIMES LL-1 ST RECURSION BASIS VECTOR
        THE NORMS OF ALL VECTORS, PSI(I), ARE CONSISTENT BUT ARBITRARY
  EMACH THRESHOLD BELOW WHICH A VALUE OF B2(I) WILL TERMINATE THE
        TRIDIAGONALISATION

Note that this routine uses :f:subr:`RECQD`, :f:subr:`CFGEN`, 
:f:subr:`RECRTS`, :f:subr:`NUMC`, :f:subr:`NUMD`.

.. f:subroutine:: TERMGN (A,B2,LL,EPS,ERR,ITMX,AA,RNG,WB,NBP1,AM,BM2,IC,WK,NW,BWK,NBD,IWK)

Generates an analytic terminator to a given continued fraction. The form of the
terminator is a sum of square roots of quadratics, F(E), as in :f:subr:`DENCRS`, 
with parameters to be adjusted to match the apparent bands gaps in the given
continued fraction. The local weight (as calculated in :f:subr:`RECWT` of F(E)
is matched to that of the given continued fraction (A(I), B2(I)) at E values
in the Neighbourhood of band edges and local minima. This routine may
serve as an example for the matching of other forms of terminating function
or matching algorithms.

::

  A    DIAGONAL RECURSION COEFFICIENTS I=1,LL-1
  B2   OFF-DIAGONAL RECURSION COEFFICIENTS I=1,LL
  LL*  LENGTH OF GIVEN RECURSION . ON OUTPUT CONTAINS THE LENGTH
       OF THE COMPUTED CONTINUED FRACTION WHICH IF DIFFERENT
       FROM INPUT INDICATES FAILURE OF CFGPGN
  
  EPS  MACHINE ACCURACY
  ERR* ACCURACY REQUIRED IN LOCATION OF BAND EDGES ,
      ON OUTPUT THE ESTIMATED ACCURACY, SUBJECT TO
  
  ITMX MAXIMUM NUMBER OF ITERATIONS IN LOCATION
  AA*  LIST OF BAND LEFT EXTREMA
  RNG* LIST OF BAND WIDTHS
  WB*  LIST OF BAND WEIGHTS
  NBP1* 1+NUMBER OF BANDS ,MAXIMUM ON INPUT AND
        ON OUTPUT CONTAINS THE 1+NUMBER COMPUTED UNLESS THIS
        EXCEEDS THE INPUT NUMBER WHEN A NEGATIVE VALUE
        INDICATES THE NUMBER OF BANDS IDENTIFIED BUT NOT
        COMPUTED. A ZERO VALUE INDICATES A FAILURE IN THE
        SEARCH PROCEDURE.(INCREASING NW MAY HELP)
  AM*  DIAGONAL C.F. COEFFICIENTS OF MATCHING FUNCTION
  BM2* OFF-DIAGONAL C.F. COEFFICIENTS OF MATCHING FUNCTION
  IC*  WORK ARRAY OF LENGTH AT LEAST NW
  WK*  WORK ARRAY OF LENGTH AT LEAST LL*2*MAX(3,NBP1)
  NW   FIRST DIMENSION OF WK. NO.OF POINTS USED IN INITIAL
       SCAN FOR BAND EXTREMA
  BWK* WORK ARRAY OF MATCHING POINTS OF DIMENSION AT LEAST 8*NBP1
  NBD  FIRST DIMENSION OF BWK : AT LEAST 2*NBP1
  IWK* INTEGER WORK ARRAY OF LENGTH AT LEAST LL

This routine uses: :f:subr:`TABAN`, :f:subr:`BNDCRD`, :f:subr:`BNDREF`, 
:f:subr:`CFGPGN`, :f:subr:`RECWT`, :f:subr:`RECQD`, :f:subr:`RECRTS`,
:f:subr:`NUMC`, :f:subr:`NUMD`, :f:subr:`WTMIN`.


.. f:subroutine:: NUMD(A,B2,ALAM,LM1,DE)

Evaluates the inverse of the logarithmic derivative of the
determinant (symmetric tridiagonal matrix - ALAM*IDENTITY)
when no sub-diagonal matrix element is zero.

::
  ARGUMENTS : ( * INDICATES AN OVERWRITTEN ARGUMENT)
  A    DIAGONAL MATRIX ELEMENTS I=1,LM1
  B2   SQUARES OF SUB-DIAGONAL MATRIX ELEMENTS I=2,LM1
  ALAM ARGUMENT OF DETERMINANT ABOVE
  LM1  DIMENSION OF MATRIX
  DE*  DET(MATRIX - ALAM)/DET'(MATRIX - ALAM)

.. f:function:: NUMC(A,B2,ALAM,LM1)

Evaluates the number of eigenvalues of a symmetric 
tridiagonal matrix strictly greater than ALAM. The sturm
sequency property is used. No sub-diagonal 
element may be zero.

::

     NUMC TAKES THE INTEGER VALUE REQUIRED
     A    DIAGONAL MATRIX ELEMENTS I=1,LM1
     B2   SQUARES OF THE SUBDIAGONAL MATRIX ELEMENTS I=2,LM1
     ALAM POINT OF EVALUATION
     LM1  DIMENSION OF MATRIX

.. f:subroutine:: BNDCRD(ET,IC,WTT,NET,NB,AA,RNG)

Estimates band edge positions from tabular information
coded as in :f:subr:`TABAN`. This consists of a list of
points where the function value either crosses a lower
or upper threshold or has a local extremum between them 
and the estimation is a straightforward analysis of 
the tabular data as generated by :f:subr:`TABAN`.

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  ET   LIST OF ORDINATES OF EXTREMA (AS O/P FROM  TABAN)
  IC*  CODE CHARACTERISING THE EXTREMA (AS O/P FROM TABAN)
       ON INPUT:
       A + SIGN INDICATES INCREASING FUNCTION VALUES AND
       A - SIGN DECREASING ONES. AT LOCAL
       EXTREMA THIS EXTENDS TO + INDICATING A MINIMUM AND - A
       MAXIMUM. THE ABSOLUTE VALUE IS CODED BELOW.
       1    FUNCTION VALUE CROSSES LOWER THRESHOLD
       2    FUNCTION VALUE CROSSES UPPER THRESHOLD
       3    FUNCTION VALUE CROSSES BOTH THRESHOLDS
       4    FUNCTION VALUE IS A LOCAL EXTREMUM

       ON OUTPUT: (ONLY IF -4 ON INPUT)
       0    CORRESPONDING LOCAL MINIMUM HAS BEEN USED IN THE
            CONSTRUCTION OF THE BAND POSITIONS

  WTT  IF IABS(IC(I))=4 THIS GIVES THE LOCAL EXTREMUM VALUE
  NET  NUMBER OF VALUES TABULATED IN ET,IC,WTT
  NB*  MAXIMUM NUMBER OF BANDS EXPECTED
       ON OUTPUT CONTAINS THE NUMBER OF BANDS IDENTIFIED
            OR IF 0 THEN A MISMATCH OF LEFT AND RIGHT EXTREMA
            OR NEGATIVE THEN -NB BANDS HAVE BEEN IDENTIFIED,
            AND THAT NUMBER IS GREATER THAN THE INPUT NB

  AA*  LIST OF LEFT BAND EXTREMA
  RNG* LIST OF BAND WIDTHS

.. f:subroutine:: BNDEST(A,B2,LL,EPS,AA,RNG,NB,EV,FEV,IC,NET,DE,WK,NW)

Estimates the band edge positions of a continued fraction representation of
a density of states, given the maximum number of bands present. This is done
by tabulation of the local quadrature weight followed by an analysis of local
minima and threshold in that table. The lower threshold is the value of the local
weight at the edge of a simple elliptical band. The accuract of the band edge location
is determined by the number, NW, of points in the tabulation.

::

  A    DENOMINATOR CONTINUED FRACTION COEFFICIENTS ; I=1,LL-1
  B2   NUMERATOR CONTINUED FRACTION COEFFICIENTS ; I=1,LL
  LL   LENGTH OF CONTINUED FRACTION
  EPS  THRESHOLD ACCURACY IN COMPUTATIONS
  AA*  LIST OF ESTIMATED BAND LEFT EXTREMA
  RNG* LIST OF ESTIMATED BAND WIDTHS
  NB*  ON INPUT: MAXIMUM NUMBER OF BANDS
       ON OUTPUT:IF >0 NUMBER OF BANDS INDENTIFIED
                 IF <OR= 0 THEN BNDCRD FAILED AND NB CONTAINS
                 ITS FAILURE CODE
  EV*  LIST OF 'MATCHING POSITIONS' FOR USE IN BNDREF
  FEV* LIST OF 'MATCHING VALUES' FOR USE IN BNDREF
  IC*  LIST OF CODES FOR THE ABOVE FOR USE IN BNDREF :
      -1 CORRESPONDS TO A  LEFT BAND EDGE
       0 CORRESPONDS TO A MINIMUM IN A BAND GAP
       1 CORRESPONDS TO A RIGHT BAND EDGE
  NET* ON INPUT : DIMENSION OF EV,FEV,IC , >OR= 2*NB
       ON OUTPUT: IF >0 NUMBER OF MATCHING POINTS
       IF <OR=0 THE TABAN FAILED AND NET HAS ITS FAILURE CODE

  DE*  TABULATION INTERVAL
  WK*  WORK ARRAY USED IN TABULATION
  NW   FIRST DIMENSION OF ARRAY WK, ALSO NUMBER OF POINTS TABULATED
       IN SEARCH FOR BAND EXTREMA

.. f:subroutine:: BNDREF(DEL,AM,BM2,LL,EPS,AA,RNG,NB,BWK,NBD,IC,NET)

Refines the approximation to band extrema in a density of states
function by a single matching to a superposition of bands
with assumed 'square-root' behaviour at their extrema. All
densities of states are represented by their continued fraction
expansion. All band-edge positions are change by + or - DEL.

The local weight of A the model local density of states, of the 
the form F(E) in :f:subr:`DENCRS`, is evaluated at a set of
matching points, and the estimated band extrema are modified
by + or - DEL in order to obtain better agreement between the local
weights evaluated from the two continues fractions.


This routine uses :f:subr:`RECWT`, :f:subr:`WTMIN`.

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  DEL  MAGNITUDE OF CHANGES REQUIRED IN BAND EDGE POSITIONS
  AM   DENOMINATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL-1
  BM2  NUMERATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL
  LL   LENGTH OF RECURSION (A(LL) UNDEFINED)
  EPS  MACHINE ACCURACY
  AA*  ESTIMATED LEFT EXTREMA OF BANDS
  RNG* ESTIMATED WIDTH OF BANDS
  NB   NUMBER OF BANDS
  BWK* MATCHING POINTS FOR MODEL DOS , AS O/P FROM BNDCRD
       FIRST COL CONTAINS ORDINATES
       SECOND COL CONTAINS ABSCISSAE
       LAST COL IS COMPUTED FUNCTION TO BE MADE ZERO, SUCH
       THAT ITS SIGN IS ALSO THE SIGN OF THE IMPOSED CHANGE
       IN ORDINATE
  NBD  FIRST DIMENSION OF BWK
  IC   CODES FOR MATCHING POINTS :0 FOR MIN,-1 FOR LEFT & +1
       FOR RIGHT BAND EDGES
  NET  NUMBER OF MATCHING POINTS

.. f:subroutine:: BNDWT(AA,RNG,WB,NB,A,B2,LL,EPS,WK,NW,IWK)

Estimates the weights of the separate bands in a DOS, given
its continued fraction coefficients, and the band edge locations.
The integrated density of states in each band gap is
calculated as the mean of the lower bound at the top
of the gap and the upper bound at the bottom of the gap,
estimated using quadrature approach (Nex).

This routines uses :f:subr:`DENQD`, :f:subr:`RECWT`, :f:subr:`NUMC`,
:f:subr:`NUMD`.

::

  ARGUMENTS ( * INDICATES AN OVERWRITTEN ARGUMENT)
  AA   LIST OF LEFT BAND EDGES ; I=1,NB
  RNG  LIST OF BAND WIDTHS ; I=1,NB
  WB*  COMPUTED BAND WEIGHTS
  NB   NUMBER OF BANDS
  A*   DENOMINATOR CONTINUED FRACTION COEFFICIENTS I=1,LL-1
        A(LL) OVERWRITTEN
  B2   NUMERATOR CONTINUED FRACTION COEFFICIENTS I=1,LL
  LL   LENGTH OF CONTINUED FRACTION
  EPS  MACHINE ACCURACY
  WK*  WORK ARRAY OF SIZE AT LEAST LL*4
  NW   FIRST DIMENSION OF ARRAY WK
  IWK* INTEGER WORK ARRAY OF LENGTH AT LEAST LL

.. f:subroutine:: CFGPGN(AA,RNG,WB,NBP1,IC,EPS,A,B2,LM1,WK,NW)

Generates coefficients of the J-matrix or 3-term recurrence 
corresponding to a linear combination of the same type
of weight function. i.e. to a set of bands in a 
density of states. The algorithm is the same as that used in 
:f:subr:`CFGEN`, but the summation in the inner product now
runs in addition over the separate bands:

.. math::
  \langle F, G \rangle = \sum_{J=1}^{NB} WB(J) \sum_{I}^{N} W(I)F(X(I,J))G(X(I,J))

where  :math:`X(I,J) = X(I) * RNG(J) + AA(J)`. The X(I) and W(I) may be
set explicitly by the user or there are provided in the routine means
of generating them for weight functions of the form :math:`\sqrt{1-X^{2}}`
(Chebyshev polynomialsof the second kind) or a constant (Legendre polynomials).
Haydock & Nex.

This routine uses :f:subr:`RECQD`, :f:subr:`RECRTS`, :f:subr:`NUMC`, :f:subr:`NUMD` if IC=2.


::

  AA   LIST OF LEFT EXTREMA OF THE BANDS
  RNG  LIST OF THE WIDTHS OF EACH BAND
  WB   LIST OF WEIGHTS OF EACH BAND
  NBP1 1+ NUMBER OF BANDS
  IC*  CODE INDICATING TYPE OF BANDS :
         1    SQRT(1.0-E*E)  -1<E<1 TYPE
         2    CONSTANT       -1<E<1 TYPE
         3    NODES AND WEIGHTS INPUT IN WK(I,J,1),J=1,2 ,I=1,LM1
         0 ON OUTPUT IF REQD FAILED (ONLY IF IC=2 ON INPUT)
  EPS  ZERO THRESHOLD (USED IN RECWTS IF IC=2 & FOR TERMINATION)
  A*   DIAGONAL COEFFICIENTS OF J-MATRIX I=1,LM1
  B2*  OFF-DIAGONAL COEFFICIENTS OF J-MATRIX I=1,LM1
  LM1*  LENGTH OF RECURRENCE OR DIMENSION OF J-MATRIX REQUIRED
        ON OUTPUT CONTAINS LENGTH GENERATED
  WK*  WORK ARRAY OF LENGTH AT LEAST LM1*2*NBP1 ;INPUT WHEN IC=3,
       & ON OUTPUT WK(I,J,1),J=1,2 ,I=1,LM1 CONTAINS THE
       QUADRATURE RULE FOR THE BASIC SINGLE BAND, NORMALISED
       TO UNITY AND DEFINED FOR ARGUMENT IN (0,1)
  NW   FIRST DIMENSION OF WK IN CALLING ROUTINE, AT LEAST LM1


.. f:subroutine:: CFGEN(X,W,N,EPS,A,B2,LM1,WK)

Calculates the three-term recurrence relation corresponding
to a weight function defined in terms of a summation inner 
product:

.. math::
  \langle F,G \rangle = \sum_{I}^{N} F(X(I))G(X(I)) W(I).

Then, with :math:`P(X,0)=1` and :math:`P(X,-1)=0`, for
:math:`N=1,2,...LM1`

.. math::
  A(N) = \frac{\langle X P(X,N-1), P(X,N-1) \rangle}{\langle P(X,N-1), P(X,N-1)\rangle}\\
  B2(N) = \frac{\langle P(X,N-1) , P(X,N-1) \rangle}{\langle P(X,N-2), P(X,N-2) \rangle}\\
  P(X,N) = (X-A(N))P(X,N-1) - B2(N)P(X,N-2)

The polynomials :math:`P(X,N)` are represented by their values
at the points X(I), with an arbitrary normalization to preserve
numerical stability. (Haydock & Nex)

::
  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  X    NODES OF THE SUMMATION INNER PRODUCT
  W    WEIGHTS OF THE SUMMATION INNER PRODUCT
  N    NUMBER OF TERMS IN THE SUMMATION INNER PRODUCT
  EPS  A VALUE OF B2(I) LESS THAN EPS TERMINATES ROUTINE
  A*   DIAGONAL COEFFICIENTS IN THE RECURRENCE A(I),I=1,LM1
  B2*  OFF-DIAGONAL COEFFICIENTS IN THE RECURRENCE B2(I),I=1,LM1
  LM1* ON INPUT  REQUIRED LENGTH OF RECURRENCE
       ON OUTPUT LENGTH OF RECURRENCE ACTUALLY GENERATED
  WK   WORK ARRAY OF LENGTH AT LEAST 2*N

.. f:subroutine:: WTMIN(AA,BB,A,B2,LL,EPS,ACC,ITMX,EM,FEM)

Finds a local minimum in the local weight, as defined in :f:subr:`RECWT`,
of a density of states by simple subdivision of a given interval. The
density of states is defined by its continued fraction or J-matrix.

This routine call :f:subr:`RECWT`.

::

  AA   LEFT END OF INTERVAL IN WHICH A MINIMUM IS LOCATED
  BB   RIGHT END OF INTERVAL IN WHICH A MINIMUM IS LOCATED
  A    DENOMINATOR COEFFICIENTS OF CONTINUED FRACTION; I=1,LL-1
  B2   NUMERATOR COEFFICIENTS OF CONTINUED FRACTION; I=1,LL
  LL   LENGTH OF CONTINUED FRACTION
  EPS  MACHINE ACCURACY
  ACC  ABSOLUTE ACCURACY REQUIRED IN LOCATION OF MINIMUM SUBJECT TO
  ITMX* MAXIMUM NUMBER OF FUNCTION EVALUATIONS
     ON OUTPUT:
     -2 IF MINIMUM OF THREE EQUIDISTANT POINTS IS AT BB
     -1 IF MINIMUM OF THREE EQUIDISTANT POINTS IS AT AA
      0 IF INSUFFICIENT ACCURACY AFTER ITMX FUNCTION EVALUATIONS
     OTHERWISE CONTAINS THE NUMBER OF FUNCTION EVALUATIONS USED
  EM*  COMPUTED LOCATION OF MINIMUM
  FEM* COMPUTED MINIMUM WEIGHT

.. f:subroutine:: SETUP(CRD,ND,NAT,EV,NTYPE,IZP,MM,NN,NND,NM,HCAL,NGBR,IOVPAR,EE,NP,NED,NE,VEC,IW)

Assembles the Hamiltonian matrix from the user supplied routines EV, HCAL, NGBR, IOVPAR and
the library routines :f:subr:`NNCAL` and :f:subr:`MMCAL`.

::

  ARGUMENTS OF SETUP : (* INDICATES OVERWRITTEN BY THE ROUTINE)

  CRD   LATTICE COORDINATES
  ND    FIRST DIMENSION OF CRD
  NAT   NO.OF ATOMS IN THE CLUSTER
  EV    LOGICAL FUNCTION OF 2 ARGUMENTS, BOTH REAL ARRAYS OF LENGTH 3
           RETURNING THE VALUE .TRUE. IF THE ARRAYS 
           ARE EQUIVALENT AND .FALSE. IF NOT.
  
  NTYPE NO. OF DIFFERENT 'TYPES' OF ATOMS
  IZP   'TYPE' OF EACH ATOM
  MM*   IS THE INTERACTION MAP GENERATED BY MMCAL
  NN*   IS THE NEIGHBOUR MAP GENERATED BY NNCAL
  NND   FIRST DIMENSION OF ARRAYS MM & NN
  NM*   MAX NO. OF ATOMS CONNECTED BY INTERACTIONS.  ON OUTPUT
        CONTAINS ACTUAL MAX NO. GENERATED
  
  HCAL  NAME OF A SUBROUTINE TO CALCULATE THE  INTERACTION BETWEEN
        TWO ATOMS. ARGUMENTS ARE
            V     VECTOR POSITION(I) - POSTITION(J)
            II    TYPE AT I
            JJ    TYPE AT J
            E*    OUTPUT INTERACTION MATRIX
                      H OPERATING ON PSI(J) EFFECT AT I
         IOVPAR    NAME OF FUNCTION SUPPLING INFORMATION TO HCAL
  
  NGBR  NAME OF A FUNCTION TO SUPPLY INTERACTION INFORMATION TO NNCAL
        ARGUMENTS ARE :
            II    'TYPE' OF ATOM I
            JJ    'TYPE' OF ATOM J
            R2    SQUARE OF DISTANCE FROM I TO J
            DD    DUMMY ARGUMENT
            NGBR  TAKES THE VALUE 1 IF I & J ARE NEIGHBOURS
                  AND 0 OTHERWISE
  
  EE    LIST OF INTERACTION MATRICES
  NP    FIRST 2 DIMENSIONS OF ARRAY EE
  NED   LAST DIMENSION OF ARRAYS EE,IW,VEC
  NE*   NO. OF DISTINCT DISPLACEMENT VECTORS (MATRICES) FOUND
  VEC*  LIST OF DISTINCT DISPLACEMENT VECTORS FOUND 
        (POSN. J - POSN.I)
  IW*   LIST OF ATOM TYPES AT THE ENDS OF THE VECTORS IN VEC
        IW(1,.) IS TYPE OF I IW(2,.) IS TYPE OF J

.. f:subroutine:: NNCAL(CRD,NDIM,NAT,IZP,NN,ND,NM,NGBR)
  
Calculates the 'NEAREST NEIGHBOUR' map of a lattice, given
a subroutine defining 'neighbour'. It also extends a map
generated by a previous call, in which case added atoms
are indicated by a negative value of IZP.


::

  ARGUMENTS: (* INDICATES OVERWRITING BY THE SUBROUTINE)
  
  CRD(I,J)  LATTICE COORDINATES (I=1,3),J=1,NAT
  NDIM      FIRST DIMENSION OF ARRAY CRD >OR= 3
  NAT       NUMBER OF LATTICE POINTS
  IZP       INTEGER*2 ARRAY LISTING THE 'TYPE' OF EACH SITE (FOR NGBR)
            IF IZP(I) IS NEGATIVE THE ABSOLUTE VALUE IS TAKEN
            AND ONLY THOSE ATOMS WITH NEGATIVE IZP ARE CONSIDERED
            FOR MODIFICATIONS TO NN
  
  NN*       'NEAREST NEIGHBOUR MAP' :
             NN(I,1) = 1+NUMBER OF NEIGHBOURS OF SITE I
             NN(I,J),J=2,NN(I,1) LIST OF SITES CONNECTED TO SITE I
  
  ND        FIRST DIMENSION OF ARRAY NN
  NM*       SECOND DIMENSION OF ARRAY NN (MAX. NO. OF NEIGHBOURS +1)
            ON OUTPUT CONTAINS ACTUAL MAX.NO. OF NEIGHBOURS +1


.. f:subroutine:: ADDAT(CRD,ND,NAT,EV,IZP,MM,NN,NND,NM,NGBR,NE,EE,NP,VEC,IW,NED,OVPAR,HCAL)

Extends the Hamiltonian matrix from the user supplied routines EV, HCAL, NGBR and IOVPAR,
and the library routines :f:subr:`NNCAL` and :f:subr:`MMCAL`. This assumes it has already
been set up by subroutine :f:subr:`SETUP` in the arrays, MM, NN, EE, VEC, and IW.

::

  ARGUMENTS OF ADDAT : (* INDICATES OVERWRITTEN BY THE ROUTINE)
  
  CRD  LATTICE COORDINATES
  ND   FIRST DIMENSION OF CRD
  NAT  NO.OF ATOMS IN THE CLUSTER
  EV   LOGICAL FUNCTION OF 2 ARGUMENTS, BOTH REAL ARRAYS OF LENGTH 3
       RETURNING THE VALUE .TRUE. IF THE ARRAYS ARE EQUIVALENT
       AND .FALSE. IF NOT.
  
  IZP  THE ABSOLUTE VALUE GIVES 'TYPE' OF EACH ATOM
          IF THE SIGN IS + THEN THE ATOM IS ASSUMED 
          PART OF THE ORIGINAL CLUSTER
          IF THE SIGN IS - THEN THE ATOM  HAS ITS 
          CONNECTIVITY AND INTERACTIONS COMPUTED
  
  MM*  IS THE INTERACTION MAP GENERATED BY MMCAL
  NN*  IS THE NEIGHBOUR MAP GENERATED BY NNCAL
  NND  FIRST DIMENSION OF ARRAYS MM & NN
  NM*  MAX NO. OF ATOMS CONNECTED BY INTERACTIONS.  ON OUTPUT
       CONTAINS ACTUAL MAX NO. GENERATED
  NGBR  NAME OF A FUNCTION TO SUPPLY INTERACTION INFORMATION TO NNCAL
        ARGUMENTS ARE :
            II    'TYPE' OF ATOM I
            JJ    'TYPE' OF ATOM J
            R2    SQUARE OF DISTANCE FROM I TO J
            DD    DUMMY ARGUMENT
  
         NGBR  TAKES THE VALUE 1 IF I & J ARE NEIGHBOURS
               AND 0 OTHERWISE
  
  NE*   NO. OF DISTINCT DISPLACEMENT VECTORS (MATRICES) ALREADY FOUND
        ON OUTPUT CONTAINS THE NEW TOTAL NUMBER FOUND
  EE*   LIST OF INTERACTION MATRICES
  NP    FIRST 2 DIMENSIONS OF ARRAY EE
  VEC*  LIST OF DISTINCT DISPLACEMENT VECTORS FOUND (POSN. I - POSN.J)
  IW*   LIST OF ATOM TYPES AT THE ENDS OF THE VECTORS IN VEC
        IW(1,.) IS TYPE OF I IW(2,.) IS TYPE OF J
  NED    LAST DIMENSION OF ARRAYS EE,IW,VEC
  OVPAR  NAME OF A FUNCTION TO SUPPLY OVERLAP PARAMETERS TO HCAL
         ARGUMENTS ARE
            II   'TYPE' OF ATOM I
            JJ   'TYPE' OF ATOM J
            R2    SQUARE OF THE DISTANCE FROM I TO J
            DD*   OVERLAP PARAMETERS AS REQUIRED BY HCAL
                  THE NOTATION USED IS AS FOLLOWS:
                DD(1)   DD SIGMA
                DD(2)   DD PI
                DD(3)   DD DELTA
                DD(4)   PD SIGMA
                DD(5)   PD PI
                DD(6)   PP SIGMA
                DD(7)   PP PI
                DD(8)   SD SIGMA
                DD(9)   SP SIGMA
                DD(10)  SS SIGMA
                DD(11)  D SELF ENERGY
                DD(12)  P SELF ENERGY
  
  HCAL  NAME OF A SUBROUTINE TO CALCULATE THE  INTERACTION BETWEEN
    TWO ATOMS. ARGUMENTS ARE
      V    VECTOR POSITION(I) - POSTITION(J)
      II   TYPE AT I
      JJ   TYPE AT J
      E*   OUTPUT INTERACTION MATRIX
           H OPERATING ON PSI(J) EFFECT AT I
      IOVPAR    NAME OF FUNCTION SUPPLING INFORMATION TO HCAL


.. f:subroutine:: RECRTS(A,B2,LM1,EPS,XLIM,N,X,MULT,BI,NI)

::

  DIMENSION A(LM1),B2(LM1),X(LM1),MULT(LM1),BI(LM1),NI(LM1)

Computes some or all of the eigenvalues of a symmetric tridiagonal
matrix with no zero sub-diagonal elements (i.e. B2(I)>0). The method
used is bisection based on the sturm sequence property followed by Newton's
method for isolated roots. [Wilkinson]. 

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  A    DIAGONAL MATRIX ELEMENTS I=1,LM1
  B2   SQUARES OF SUB-DIAGONAL MATRIX ELEMENTS I=2,LM1
  LM1  DIMENSION OF MATRIX
  EPS  ABSOLUTE ACCURACY REQUIRED IN EIGENVALUES
  XLIM UPPER BOUND ON EIGENVALUES TO BE FOUND (IF RELEVANT)
  N*   IF 0 ON INPUT : ONLY EIGENVALUES LESS THAN XLIM ARE FOUND
       ON OUTPUT : NUMBER OF DISTINCT EIGENVALUES FOUND
  X*   EIGENVALUES IN ASCENDING ORDER
  MULT* MULTIPLICITY OF EACH EIGENVALUE
       IF NEGATIVE THEN THE CORRESPONDING EIGENVALUE WAS FOUND
       WITH LESS ACCURACY THAN EPS
  BI*  REAL WORK ARRAY OF LENGTH AT LEAST LM1
  NI*  INTEGER WORK ARRAY OF LENGTH AT LEAST LM1

This routine uses :f:subr:`NUMC` and :f:subr:`NUMD`.

.. f:subroutine:: MMCAL(CRD,NDIM,NAT,NN,ND,NM,EV,IZP,NMAT,MM,VEC,IW)

Computes an index of distinct vectors linking neighbouring sites
in a given lattice. The vectors are computed and indexed according to
the 'type' (as defined by IZP) of the terminal atoms as well as by the
vector components. Thus if there are 3 types of atoms linked in all
pair combinations by equivalent vectors, all combinations will occur
in the index. (i.e. 12 entries including both senses of the vector)
if any of the 'types' in IZP are negative, it is assumed that
:f:subr:`MMCAL` has already been called for a subcluster of the current cluster
and that those atoms with negative IZP are new additions whose
interactions are to be computed (see :f:subr:`ADDAT` for an example of this
usage).

::

  INTEGER*2 NN(ND,NM),MM(ND,NM),IZP(NAT),IW(2,NMAT)
  DIMENSION CRD(NDIM,NAT),VEC(NDIM,NMAT)
  LOGICAL EV
  COMMON /BLKNNM/NNMAT

  CRD(I,J)  COORDINATES OF THE LATTICE (I=1,NDIM) ,J=1,NAT
  NDIM    FIRST DIMENSION OF ARRAYS CRD AND VEC
  NAT     NUMBER OF SITES IN THE LATTICE
  NN      NEAREST NEIGHBOUR MAP AS CALCULATED BY NNCAL :
          NN(I,1)=1+NO.OF NEIGHBOURS OF SITE I
          NN(I,J),J=2,NN(I,1) LISTS THE NEIGHBOURS OF SITE I

  ND      FIRST DIMENSION OF ARRAY NN
  NM      SECOND DIMENSION OF ARRAY NN
  EV      LOGICAL FUNCTION (DECLARED EXTERNAL IN THE CALLING ROUTINE)
          WITH 2 ARGUMENTS, EACH A REAL ARRAY OF LENGTH NDIM, 
          RETURNING THE VALUE .TRUE. IF ITS ARGUMENTS ARE THE 'SAME'
          AND  .FALSE. IF NOT. THE ARGUMENTS MUST BE UNCHANGED.

  IZP     IZP(I) ABSOLUTE VALUE GIVES 'TYPE' OF I TH LATTICE SITE
          IF ATOMS ARE BEING ADDED TO AN EXISTING CLUSTER THEN A
          NEGATIVE SIGN INDICATES AN ADDED ATOM.

  NMAT*   ON  A FIRST CALL THE MAXIMUM NUMBER OF DISTINCT VECTORS
          ALLOWED. SUBSEQUENTLY THE NUMBER PREVIOUSLY 
          CALCULATED (AS O/P)
          ON OUTPUT THE ACTUAL NUMBER OF VECTORS CALCULATED
          IF 0 THEN NOT ENOUGH STORE HAS BEEN ALLOWED
          AND NMAT MUST BE INCREASED.

  MM*     INDEX OF VECTORS LINKING NEIGHBOURING SITES:
          MM(I,J)= K, THE INDEX OF THE VECTOR STORED IN VEC SUCH
          THAT VEC(K)=SITE VECTOR(NN(I,J)) - SITE VECTOR(I)  ,J=2,NN(I,1)

  VEC(R,K)* LIST OF DISTINCT VECTORS  ,(R=1,NDIM) , K=1,NMAT
  IW(1,K)*  'TYPE' OF ATOM I AT ONE END OF THE K TH VECTOR
  IW(2,K)*  'TYPE' OF ATOM J AT THE OTHER END OF THE K TH VECTOR


.. f:subroutine:: ONION(NN,ND,NM,IZERO,NAT,IST,NNS,IW)

::

      INTEGER*2 NN(ND,NM),IZERO(NAT),IST(NNS),IW(NAT)

Assigns each site in a lattice (defined by a 'connectivity map')
to a shell defined by a 'topological' (number of 'hops') distance from
a given group of sites. The given group is labelled 'SHELL 1'.

::

  NN     NEIGHBOUR MAP AS DEFINED BY NNCAL
  ND     FIRST DIMENSION OF ARRAY NN
  NM     SECOND DIMENSION OF ARRAY NN
  IZERO* INTEGER*2 ARRAY RETURNING THE SHELL NUMBER OF EACH SITE
  NAT    NUMBER OF LATTICE SITES
  IST    INTEGER*2 ARRAY INDEXING THE 'CENTRAL' SITE(S)
  NNS    NUMBER OF CENTRAL SITES
  IW     INTEGER*2 WORK ARRAY OF LENGTH AT LEAST NAT


.. f:subroutine:: ORPEEL(NSTRT,NORB,NO,MM,NN,ND,ID,EE,NP,NE,NED,MEM)

Implements orbital peeling as specified in the PHD thesis
of N.R. Burke. An equivalent (functional) definition is that
the subroutine deletes a row and column of a sparse matrix
as set up using :f:subr:`NNCAL` and :f:subr:`MMCAL`. The matrix
is assumed to be partitioned into NP by NP blocks, of which there 
are only relatively few distinct ones in the overall matrix. To
delete a row and colum, therefore, a copy is made of the blocks 
involved and the list of submatrices modified accordingly.
It is assumed that the overall purpose is to delete rows and 
columns defined by a given diagonal submatrix.

::

  NSTRT THE STARTING ATOM .(DIAGONAL SUBMATRIX TO BE DELETED)
  NORB  ORBITAL TO BE PEELED (ROW & COL. OF SUBMATRX TO BE DELETED)
  NO    CODE :
        IF = 1  THE INTERACTION MATRICES ARE COPIED AND EE EXTENDED
               (I.E.FIRST CALL FOR A GIVEN PEELING SEQUENCE)
        IF BETWEEN 1 & NP THE COPIED INTERACTION MATRICES ARE MODIFIED
        BY DELETION OF THE APPROPRIATE ROW OR COLUMN (THE NORBTH)
        IF = NP THE INTERACTION MATRICES ARE RESTORED TO THOSE
        ORIGINALLY OPERATIVE.(I.E. THE LAST CALL OF A
        SEQUENCE)
  MM*   THE INDEX OF SUBMATRICES CORRESPONDING TO NN
        MM(I,J)  INDEX OF INTERACTION MATRIX BETWEEN ATOM NN(I,J)
        AND ATOM I ; J.NE.1 . IF J=1 THEN = INDEX OF THE SELF
        INTERACTION MATRIX OF ATOM I.
  NN      THE INDEX OF NEIGHBOURS
            NN(I,1) = 1+ NO. OF NEIGHBOURS OF ATOM I
            NN(I,J), J=2,NN(I,1) LIST OF NEIGHBOURS OF ATOM I
  ND      FIRST DIMENSION OF ARRAYS NN & MM
  ID      SECOND DIMENSION OF ARRAYS NN & MM
  EE*     LIST OF INTERACTION MATRICES
  NP      FIRST 2 DIMENSIONS OF ARRAY EE
  NE*     NO. OF INTERACTION MATRICES SO FAR COMPUTED
  NED     MAX NO. OF INTERACTION MATRICES ALLOWED( LAST DIMENSION OF EE)
  MEM*    STORAGE SPACE TO ENABLE RESTORATION OF THE ORIGINAL MATRIX


.. f:function:: DENQD(E,EMX,A,B2,LL,ALP,EPS,TB,NT,NQ,NE,IWK)

::

  DIMENSION A(LL), B2(LL), TB(NT,4), IWK(LL)
  
Evaluates the density of states, :math:`N(E)`, corresponding to a given
continued fraction (J-Matrix) at a given point :math:`E` and returns that value
and also quadrature nodes and weights at a set of points bounded above by
EMX. The table of values TB is output so that the integrated density of states,
densiyt of states, and similar function may be evaluated
at each E(I) not greater than EMX.
e.g. The integral to TB(I,1) of F(x)N(x)dx 
is approximated by the sum J=1,I (last term times alpha)

.. math::

  F(TB(J,1))TB(J,2)

The expressions defining the approximation are as follows
(with N=LL):

.. math::
  A(N) = E- B2(N) \frac{P(E,N-1)}{P(E,N)}
  W(I) = \frac{Q(E(I),N)}{P'(E(I),N+1)}

The term :math:`\frac{DW(I)}{DA(LL)}` is equal to the following expression
evaluated at E(I):

.. math::
  \frac{Q'(N)P(N)-P'(N+1)Q(N-1)+W(I)P'(N)P'(N+1)-P''(N+1)P(N)}{P'(N+1)^{2}}

.. math::
  N(E) = P(E,N+1)/P(E,N){\sum_{I}^{NE} \frac{DW(I)}{DA(LL)} + ALP*\frac{DW(NE)}{DA(LL)}}

P and Q are the monic orthogonal polynomials of the first 
and second kind associated with the weight function
N(E) (see :f:subr:`PLYVAL` for explicit definition of P),
and the E(I) are the eigenvalues of the given Jacobi
matrix with A(LL) appended so that E(NE)=E.
In the actual calculation the values of the polynomials
are renormalised to maintain numerical stability
(only ratios of polynomials appear in the above expressions).

This routines uses :f:subr:`RECWT`, :f:subr:`RECRTS`,
:f:subr:`NUMC`, :f:subr:`NUMD`.

::

  DENQD TAKES THE COMPUTED VALUE OF THE DENSITY OF STATES AT E
  E    VALUE AT WHICH DENSITY OF STATES REQUIRED
  EMX  UPPER LIMIT OF RANGE OF QUADRATURE NODAL VALUES REQUIRED > E
  A*   DIAGONAL J-MATRIX ELEMENTS (A(LL) OVERWRITTEN) I=1,LL-1
  B2   SQUARES OF SUB-DIAGONAL J-MATRIX ELEMENTS I=2,LL
       B2(1) IS THE TOTAL WEIGHT OF THE DENSITY OF STATES
  LL   LENGTH OF TRIDIAGONALISATION
  ALP  PROPORTION OF WEIGHT AT LAST NODE, 0<ALP<1 ,USUALLY =0.5
  EPS  ACCURACY REQUIRED IN ROOT-FINDING
  TB*  TABLE OF QUADRATURE NODES AND DIFFERENTIALS :
       TB(I,1)   NODAL POINTS : E(I)
       TB(I,2)   NODAL WEIGHTS : W(I)
       TB(I,3)   DW(I)/DA(LL)
       TB(I,4)   P'(E(I),LL+1) / P(E(I),LL)
  NT   FIRST DIMENSION OF ARRAY TB IN CALLING SEGMENT
  NQ*  NUMBER OF NODAL VALUES CALCULATED
  NE*  IABS(NE) GIVES INDEX OF NODE CORRESPONDING TO E :TB(NE,1)=E
       IF NEGATIVE THE ACCURACY IS INADEQUATE
       IF = 0 A MULTIPLE ROOT WAS IDENTIFIED IN RECRTS
  IWK* INTEGER WORK SPACE OF LENGTH AT LEAST LL (O/P FROM RECRTS)


.. f:function:: DENSQ(E,A,B2,LL,EI)

::

  DIMENSION A(LL),B2(LL),EI(2),P(2),Q(2)

Computes the value of the local density of states
corresponding to a continue fraction, using the 
square root terminator.(Haydock)

.. math::
  N(E) = \frac{-1}{\pi}{\rm Im}[\frac{Q(E,N-1)-B2(N)T(E)Q(E,N-2)}{P(E,N)-B2(N)T(E)P(E,N-1)}]\\
  T(E) = 0.5E - {EI(1)+EI(2)}*0.5-\sqrt{E-EI(1)}\frac{\sqrt{EI(2)-E}}{B2(LL)}

P and Q are the corresponding orthogonal polynomials of the first
and second kinds.

::
  DENSQ TAKES THE REQUIRED VALUE
  E    ARGUMENT OF CONTINUED FRACTION
  A    DENOMINATOR COEFFICIENTS OF CONTINED FRACTION I=1,LL-1
  B2   NUMERATOR COEFFICIENTS OF CONTINED FRACTION I=1,LL
  LL   LENGTH OF CONTINED FRACTION
  EI   BAND EDGES


.. f:function:: DENCRS(E,A,B2,LL,AA,RNG,WB,NB,AM,BM2)

Computes the value of a continued fraction using a terminator
based on the number, weights and positions of separate bands using
a general prescription (Haydock and Nex- To Appear). 
The matching continued fraction with square 
root band edges may be generated using :f:subr:`CFGPGN` 
or :f:subr:`TERGEN` and should be of the same length 
as the original.

The function:

.. math::
  F(E) = \sum_{K}8.0\frac{WB(K)}{RNG(K)^{2}}(E-(AA(K)+ 0.5\cdot RNG(K))-\sqrt{E-AA(K)}\sqrt{AA(K)+RNG(K)-E}

is assumed to correspond to the supplied coefficients AM(I) and BM2(I).

.. math::
  T(E) = \frac{S(E,N-1)-F(E)R(E,N)}{S(E,N-2)-F(E)R(E,N-1)}\frac{1}{BM2(N)} \\

  N(E) = \frac{-1}{\pi}{\rm Im}[Q(E,N-1)-\frac{B2(N)T(E)Q(E,N-2)}{P(E,N)-B2(N)T(E)P(E,N-1)}]

where N=LL and P,Q and R,S are the orthogonal polynomials of the first and second kinds 
corresponding to A,B2 and AM,BM2 respectively. 

This routine uses :f:subr:`PLYVAL`.

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  
  DENCRS TAKES THE REQUIRED VALUE
  E    ARGUMENT OF CONTINUED FRACTION
  A    DENOMINATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL-1
  B2   NUMERATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL
  LL   LENGTH OF CONTINUED FRACTION
  AA   LIST OF BAND LEFT EXTREMA
  RNG  LIST OF BAND WIDTHS
  WB   LIST OF WEIGHTS OF BANDS
  NB   NUMBER OF BANDS (GREATER THAN 0)
  AM   LL-1 DENOMINATOR COEFFICIENTS OF MATCHING CONTINUED FRACTION
  BM2  LL NUMERATOR COEFFICIENTS OF MATCHING CONTINUED FRACTION


.. f:function:: DENINT(E,A,B2,NA,NP,LL,ALP,EPS,WK,IWK,ICODE)

::

  DIMENSION A(NA,NP),B2(NA,NP),WK(LL,4),IWK(LL)

Evaluates the integrated density of states, N(E), 
corresponding to a given sum of continued fractions
(J-matrices) at a given point E and returns that value,
using the 'quadrature' approach. This routine
uses :f:subr:`DENQD`, :f:subr:`RECWT`,
:f:subr:`RECRTS`, :f:subr:`NUMC`, :f:subr:`NUMD`.

::

  DENINT TAKES THE COMPUTED VALUE OF THE INTEGRATED DENSITY OF STATES AT E
  E      VALUE AT WHICH INTEGRATED DENSITY OF STATES REQUIRED
  A*     DIAGONAL J-MATRIX ELEMENTS (A(LL,K) OVERWRITTEN) I=1,LL-1
  B2     SQUARES OF SUB-DIAGONAL J-MATRIX ELEMENTS I=2,LL
         B2(1,K) IS THE WEIGHT IN THE K TH BAND
  NA     FIRST DIMENSION OF ARRAYS A AND B2 >= LL
  NP     NUMBER OF DENSITY OF STATES TO BE SUMMED
  LL     LENGTH OF TRIDIAGONALISATIONS
  ALP    PROPORTION OF WEIGHT AT LAST NODE, 0<ALP<1 ,USUALLY =0.5
  EPS    ACCURACY REQUIRED IN ROOT-FINDING
  WK*    WORK ARRAY OF LENGTH AT LEAST 4*LL
  IWK*   INTEGER WORK SPACE OF LENGTH AT LEAST LL (O/P FROM RECRTS)
  ICODE* 0 ON A SUCCESSFUL OUTPUT
  NEGATIVE  IF A FAILURE IN DENQD



.. f:subroutine:: DENCRQ(E,A,B2,LL,AA,RNG,WB,NB,AM,BM2)

Computes the value of a Greenian represented by a continued fraction
using a terminator based on the number, weights, and positions of
separate bands using a general prescription 
(Haydock and Nex to appear). The matching continued fraction with
square root band edges may be generated using 
:f:subr:`CFGPGN` or :f:subr:`TERGEN`
and should be of the same length as the original.

The function:

.. math::
  F(E) = \sum 8.0 \frac{WB(K)}{RNG(K)^{2}} E-AA(K) + 0.5RNG(K)
        - \sqrt{E-AA(K)}\sqrt{AA(K)+RNG(K)-E},

is assumed to correspond to the supplied coefficients AM(I),
BM2(I).

.. math::
  T(E) = \frac{S(E,N-1)-F(E)R(E,N)}{S(E,N-2)-F(E)R(E,N-1)}\frac{1}{BM2(N)} \\
  N(E) = \frac{-1}{\pi}Im[G(E)] \\
  G(E) = Q(E,N-1)-B2(N)T(E)Q(E,N-2)/P(E,N)-B2(N)T(E)P(E,N-1),

where N=LL and P,Q and R,S are the orthogonal polynomials of the first and 
second kinds corresponding to A,B2 and AM,BM2 respectively. 

This routine uses :f:subr:`PLYVAL`.

::
  DENCRQ TAKES THE REQUIRED VALUE
  E    ARGUMENT OF CONTINUED FRACTION
  A    DENOMINATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL-1
  B2   NUMERATOR COEFFICIENTS OF CONTINUED FRACTION I=1,LL
  LL   LENGTH OF CONTINUED FRACTION
  AA   LIST OF BAND LEFT EXTREMA
  RNG  LIST OF BAND WIDTHS
  WB   LIST OF WEIGHTS OF BANDS
  NB   NUMBER OF BANDS (GREATER THAN 0)
  AM   LL-1 DENOMINATOR COEFFICIENTS OF MATCHING CONTINUED FRACTION
  BM2  LL NUMERATOR COEFFICIENTS OF MATCHING CONTINUED FRACTION


.. f:function:: RECWT(E,A,B2,LL,EPS,N,P,NS)

Computes the value of the weight at the fixed point in a 1-fixed
point Gaussian quadrature, given the corresponding 3-term recurrence
relation:

.. math::
  P(E,J)= (E-A(J))P(E,J-1) - B2(J)P(E,J-2)

::

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
  
  P* FINAL POLYNOMIAL VALUES USED IN CALCULATION OF WEIGHT TO BE USED 
     UNCHANGED IF ROUTINE IS RE-ENTERED WITH NS=LL
       IF N=LL-IABS(N)
            P(2,1)=P(E,N)       P(1,1)=P(E,N-1)
            P(2,2)=P'(E,N)      P(1,2)=P'(E,N-1)
            P(2,3)=Q(E,N-1)     P(1,3)=Q(E,N-2)
       Q(E,M) IS THE POLYNOMIAL OF THE SECOND KIND OF DEGREE M
  
  NS POINT AT WHICH RECURRENCE IS INITIATED . THIS SHOULD BE
     1 INITIALLY , BUT FOR A SUBSEQUENT CALL, WITH E UNCHANGED AND LARGER LL, 
     SHOULD BE SET TO THE CURRENT VALUE OF LL

This routine may be called repeatedly with increasing number
of levels such that it does not recompute earlier polynomial
values. If required the value of the last coefficient, A(LL),
may be computed, or it may be assumed that this has already been
done and that value used in the calculation of the weight.
The expression for the weight used is (with N=LL):

.. math::
  \frac{P(E,N-1)Q(E,N-1)-P(E,N)Q(E,N-2)}{P(E,N-1)P'(E,N)-P'(E,N-1)P(E,N)+P(E,N)^{2}/B2(N)}.

As this form is independent of the normalisation of the polynomials. P and Q are the monic
polynomials of the first and second kinds.

.. f:subroutine:: SCAN(NN,ND,NNMX,N0,NAT,NON,SUB)

Generates all neighbours (0th, 1st, and 2nd if required) of a subcluster of atoms
(consecutively numbered) of a given cluster. Input is the 'nearest neighbour' map
of the whole cluster and output is via a user supplied subroutine which is called for
each possible neighbour.

::

  NN   NEAREST NEIGHBOUR MAP. (N.B. INTEGER*2 ARRAY)
       NN(I,1) CONTAINS 1+ NO. OF NEIGHOURS OF ATOM I
       NN(I,J),J=2,..,NN(I,1) IS THE LIST OF ATOM NUMBERS
       OF THE NEIGHBOURS OF ATOM I
  
  ND    FIRST DIMENSION OF ARRAY NN
  NNMX  SECOND DIMENSION OF ARRAY NN
  N0    FIRST ATOM OF THE SUBCLUSTER WHOSE NEIGHBOURS ARE TO
        BE GENERATED
  NAT   LAST ATOM OF THAT SUBCLUSTER
  
  NON   'ORDER' OF NEIGHBOURS REQUIRED  I.E.
        1 IF FIRST NEIGHBOURS ONLY
        2 IF FIRST & SECOND NEIGHBOURS
  
  SUB   NAME OF A USER SUPPLIED SUBROUTINE (DECLARED EXTERNAL IN
        THE CALLING ROUTINE)TO PROCESS THE INFORMATION GENERATED.
        ITS ARGUMENTS , WHICH MUST NOT BE MODIFIED, ARE :
        ......... (IA,NA,NOP)
        DIMENSION IA(NOP),NA(NOP)
  
   NOP  CONTAINS THE CODE AS FOLLOWS:
         1   FOR THE SELF INTERACTION
         2   FOR A 1ST. NEIGHBOUR INTERACTION
         3   FOR A 2ND. (NEIGHBOUR OF NEIGHBOUR) INTERACTION
  
   IA(NOP) IS THE INDEX OF THE NEIGHBOUR GENERATED I.E.
       IA(1)=I
       IA(2)=INDEX OF FIRST NEIGHBOUR OF I (IF NOP>OR= 2)
       IA(3)=INDEX OF 2ND. NEIGHBOUR OF I (VIA ATOM IA(2)) IF NOP=3
  
   NA(I) IS THE SUBSCRIPT IN THE NEIGHBOUR MAP NN OF THE
         GENERATED NEIGHBOUR. I.E.
         NA(1)=1
         NA(2)=J  WHERE IA(2)=NN(I,J) (IF NOP>OR= 2)
         NA(3)=K  WHERE IA(3)=NN(J,K) (IF NOP=3)

.. f:subroutine:: RECPER(HOP,VOP,W1,W0,A,B,NW,LLIM,NA,NL,AMAT)

For a discussion of perturbation theory and the recursion method see
`J. Phys. A Vol. 10, No. 4 (1977) <http://iopscience.iop.org/article/10.1088/0305-4470/10/4/009>`_ 
and `R. Haydock, Philos. Mag. [Part] B 37, 97 (1978) <https://doi.org/10.1080/13642817808245310>`_.
See pg. 283 in SSPV 35 for a formal discussion of perturbations to the chain model and the 
change in the coefficients of the continued fraction. 

::

  ARGUMENTS (* INDICATES OVERWRITTEN BY THE ROUTINE)
  
  HOP      NAME OF A SUBROUTINE SUPPLIED BY THE USER (AND DECLARED
           EXTERNAL IN THE CALLING ROUTINE) TO CALCULATE HX+Y
           AND Y(TRANSPOSED)HX, FOR ARBITRARY MATRICES X AND Y.THE
           ARGUMENTS OF HOP MUST BE AS FOLLOWS:
  
             SUBROUTINE HOP(X,Y,A,NW,NA,LL)
             DIMENSION X(NW,LL),Y(NW,LL),A(NA,LL)
  
             X   AN NW BY LL ARRAY TO BE PROCESSED
             Y*  AN NW BY LL ARRAY TO BE PROCESSED CONTAINING Y
                 ON INPUT AND HX+Y ON OUTPUT.
             A*  THE COMPUTED MATRIX Y(TRANSPOSED)HX
             NW  FIRST DIMENSION OF MATRICES X AND Y
             NA  FIRST DIMENSION OF ARRAY A
             LL  NO. OF COLUMNS IN MATRICES X AND Y
  
         NOTE THAT ONLY THE STARRED (*) ITEMS ARE TO BE SET BY THE USER.
  
  VOP    NAME OF A SUBROUTINE  SATISFYING THE SAME CONDITIONS AS HOP
         BUT WITH V REPLACING H.
  W1*    SQRT(B(0,0))*W0 : THE STARTING VECTORS OF THE
         RECURRENCE (UNNORMALISED).THE FIRST SUBSCRIPT RUNS
         OVER THE VECTOR COMPONENTS AND THE SECOND OVER THE
         PERTURBATION SERIES.
  W0*    W(-1,K) THE NORMALISED (-1) STARTING VECTORS STORED AS W1
  A*     OUTPUT AS THE ARRAY OF A COEFFICIENTS, THE FIRST SUBSCRIPT
         RUNNING OVER THE  RECURRENCE RELATION AND THE SECOND OVER THE
         PERTURBATION SERIES.
  B*     THE SQRT(B(N,0)*B(N,K)) COEFFICIENTS STORED AS THE AS.
         B(1,K) MUST BE SET AND CONSISTENT WITH W1.
  NW     DIMENSION OF MATRICES H AND V
  LLIM   LENGTH OF PERTURBATION SERIES REQUIRED.
  NA     FIRST DIMENSION OF ARRAYS A AND B
  NL     NO. OF 'LEVELS' IN THE RECURRENCE
  AMAT*  WORK ARRAY OF AT LEAST LLIM*LLIM ELEMENTS

.. f:subroutine:: BCCLAT(CRD,NDIM,IZP,NAT,NX,NY,NZ,NTYPE)

Generates a BCC lattice on a positive integer grid, 
enclosed by a cuboid of a given size.

::

  ARGUMENTS:( * INDICATES AN OVERWRITTEN ARGUMENT)
  CRD*    LATTICE COORDINATES ((I,J),I=1,3),J=1,NAT
  NDIM    FIRST FIRST DIMENSION OF ARRAY COORD >OR= 3
  IZP*    INTEGER*2 ARRAY RETURNING THE VALUE NTYPE IN EACH ELEMENT
  NAT*    ON INPUT THE MAXIMUM NUMBER OF LATTICE POINTS ALLOWED
          ON OUTPUT THE ACTUAL NUMBER OF POINTS GENERATED
          NX,NY,NZ  INTEGER DIMENSIONS OF THE CUBOID TO CONTAIN THE LATTICE
          NTYPE   'TYPE' CODE FOR EACH LATTICE SITE

.. f:function:: BCCBFE(I,J,R2,DD)

Determines whether a distance is a 'nearest neighbour' or 'next nearest neighbour'
distance in the BCC lattice generated by :f:subr:`BCCLAT`, and if so outputs the DD 
parameters for iron according to D.G. Pettifor. 

::

  ARGUMENTS:
  I   'TYPE' OF ONE LATTICE SITE
  J   'TYPE' OF THE OTHER LATTICE
  R2   SQUARE OF THE DISTANCE BETWEEN THE TWO LATTICE SITES
  DD*  OUTPUT AS THE DD PARAMETERS OF D.G.PETTIFOR (SIGMA,PI,DELTA)
       AND DD(11)=0.0 OF R2<1.0E-4 AS THE SELF ENERGY
       BCCBFE TAKES THE VALUE 0 IF THE SITES ARE NOT NEIGHBOURS
       AND 1 IF THEY ARE NEIGHBOURS

.. f:function:: EQUIV(V,W)

::

  DIMENSION V(3),W(3)

Determines whether two vectors are 'equivalent' (for the purposes of the
subroutine :f:subr:`MMCAL`).

::
  
  INPUT:
  V    REAL ARRAY OF LENGTH 3
  W    REAL ARRAY OF LENGTH 3

  OUTPUT :
  EQUIV TAKES THE VALUE .TRUE. IF THE VECTORS ARE 'EQUIVALENT' 
  AND .FALSE. OTHERWISE



.. f:function:: FCCBND(I,J,R2,DD)

Determines whether a distance is a 'Nearest Neighbour' Distance 
in the FCC lattice generated by :f:subr:`FCCLAT`, and if so outputs
the DD parameters according to D.G. Pettifor. (In Rydbergs).

::

  I    'TYPE' OF ONE LATTICE SITE
  J    'TYPE' OF THE OTHER LATTICE
  R2   SQUARE OF THE DISTANCE BETWEEN THE TWO LATTICE SITES
  DD*  OUTPUT AS THE DD PARAMETERS OF D.G.PETTIFOR (SIGMA,PI,DELTA)
       AND DD(11) IS OUTPUT AS 0.0  (THE SELF-ENERGY) IF R2<1.0E-4
       FCCBND TAKES THE VALUE 0 IF THE SITES ARE NOT NEIGHBOURS
       AND 1 IF THEY ARE NEIGHBOURS

 
.. f:subroutine:: TABAN(E,WT,NPTS,THU,THL,ET,IC,WTT,NET)

::

  DIMENSION E(NPTS),WT(NPTS),ET(NET),IC(NET),WTT(NET)

Identifies extremal values of a tabulated function, within a
user-defined range of function values. This is effected by
simple comparison of the tabular values.

::

  ARGUMENTS : (* INDICATES AN OVERWRITTEN ARGUMENT)
  E    LIST OF ORDINATE VALUES
  WT   LIST OF ABSCISSAE
  NPTS NUMBER OF COORDINATE PAIRS
  THU  UPPER VALUE OF 'WINDOW' ON FUNCTION
  THL  LOWER VALUE OF 'WINDOW' ON FUNCTION
  ET*  LIST OF ORDINATES OF EXTREMA WITHIN WINDOW
  IC*  CODE OF TYPE OF EXTREMA : A + SIGN INDICATES INCREASING
          FUNCTION VALUES AND -, DECREASING ONES. AT LOCAL
          EXTREMA THIS EXTENDS TO + INDICATING A MINIMUM AND - A
          MAXIMUM. THE ABSOLUTE VALUE IS CODED BELOW.
          1    FUNCTION VALUE CROSSES LOWER THRESHOLD
          2    FUNCTION VALUE CROSSES UPPER THRESHOLD
          3    FUNCTION VALUE CROSSES BOTH THRESHOLDS
          4    FUNCTION VALUE IS A LOCAL EXTREMUM
  WTT* IF IABS(IC(I))=4 THIS GIVES THE LOCAL EXTREMUM VALUE
  NET* NUMBER OF VALUES TABULATED IN ET,IC,WTT
       IF NEGATIVE ON OUTPUT, THEN THERE IS NOT ENOUGH SPACE
       AND THE ABSOLUTE VALUE GIVES THE INDEX OF THE LAST
       COORDINATE PAIR EXAMINED. THE RESULTS TO THAT POINT
       ARE STORED AS ABOVE.

.. f:subroutine:: FCCLAT(CRD,NDIM,IZP,NAT,NX,NY,NZ,NTYPE)

Generates a FCC lattice on positive integer grid, enclosed by a cuboid
of a given size.

::

  ARGUMENTS:( * INDICATES AN OVERWRITTEN ARGUMENT)
  CRD*    LATTICE COORDINATES ((I,J),I=1,3),J=1,NAT
  NDIM    FIRST FIRST DIMENSION OF ARRAY COORD >OR= 3
  IZP*    INTEGER*2 ARRAY RETURNING THE VALUE NTYPE IN EACH ELEMENT
  NAT*    ON INPUT THE MAXIMUM NUMBER OF LATTICE POINTS ALLOWED
          ON OUTPUT THE ACTUAL NUMBER OF POINTS GENERATED
  NX,NY,NZ  INTEGER DIMENSIONS OF THE CUBOID TO CONTAIN THE LATTICE
  NTYPE   'TYPE' CODE FOR EACH LATTICE SITE

.. f:subroutine:: PEEL(CRD,NDIM,NAT,NN,ND,NM,IST,NS,IZP,IZERO,NSH,IW)

Given a cluster of 'sites' and its 'shell' structure, retains only those 'sites'
within a given number of 'shells'. The accepted sites are then renumbered as are
other relevant reference arrays.

::

  ARGUMENTS :  (* INDICATES AN OVERWRITTEN ARGUMENT)
  CRD*      LIST OF LATTICE COORDINATES
  NDIM      FIRST DIMENSION OF ARRAY CRD
  NAT*      ON INPUT THE NUMBER OF LATTICE SITES
            ON OUTPUT THE NUMBER OF RETAINED LATTICE SITES
  NN*       INTEGER*2 NEAREST NEIGHBOUR MAP (AS O/P FROM NNCAL)
  ND        FIRST DIMENSION OF ARRAY NN
  NM        SECOND DIMENSION AF ARRAY NN
  IST*      INTEGER*2  ARRAY LISTING THE 'CENTRAL' SITES
  NS        NUMBER OF SITES LISTED IN IST
  IZP*      INTEGER*2 'TYPE' OF EACH SITE
  IZERO*    INTEGER*2 'SHELL NUMBER' OF EACH SITE
  NSH       NUMBER OF SHELLS TO BE RETAINED
  IW*       INTEGER*2 IW(I) =NEW INDEX OF OLD SITE I

.. f:subroutine:: SLKODE(DUM,I,J,EM,IBONDS)

::
  
  DIMENSION DUM(3),EM(5,5)

Calculates the 5x5 2 centre interaction matrices between the 
D-electrons on two sites a given vector apart. The prescription used
is that of Slater & Koster in Phys. Rev. 94 Vol. 6 P. 1498 et. seq..
The orbitals are ordered: :math:`xy, yz, zx, x^{2}-y^{2,} 3x^{2}-r^{2}`.
If DUM has squared modulus :math:`< 10^{-4}` then a diagonal
self interaction matrix is returned. 

This routines call :f:subr:`SKDD` and :f:subr:`SELFD`.

::

  ARGUMENTS (* INDICATES OVERWRITTEN BY THE ROUTINE)
  DUM    THE VECTOR B-A
  I      THE 'TYPE' OF ATOM A
  J      THE 'TYPE' OF ATOM B
  EM*    THE 5 X 5 MATRIX SUCH THAT
         (HY(B))(A) = EM Y(B)

  IBONDS THE NAME OF AN INTEGER FUNCTION SUPPLIED BY THE CALLING
  ROUTINE AND DECLARED EXTERNAL WITH THE FOLLOWING
  ARGUMENTS :
    I    'TYPE' OF ATOM A
    J    'TYPE' OF ATOM B
   R2     SQUARE OF DISTANCE FROM A TO B
   DD*    DD(SIGMA,PI,DELTA) FOR THE ABOVE ARGUMENTS
   DD(11) CONTAINS THE SELF ENERGY IF R2 < 1.0E-4


