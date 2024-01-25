!******************************************************************************
!                                                                             *
!                                                                             *
!                                                                             *
!               DIFFERENTIAL ALGEBRA PACKAGE OF M. BERZ                       *
!                     ****************************                            *
!                               holy3                                         *
!                                                                             *
!                                                                             *
!                                                                             *
!         VERSION FOR MACHINE IN LINE THAT IS NOT COMMENTED OFF               *
!        TO CREATE DIFFERENT VERSIONS, USE THE PROGRAM 'VERSION'              *
!                                                                             *
!                                                                             *
!                                                                             *
!                                                                             *
!        THIS PACKAGE WAS INITIALLY WRITTEN BY PROF. M. BERZ WHILE AT         *
!        THE LAWRENCE BERKELEY LABORATORY.                                    *
!        IT HAS BEEN EXTENSIVELY MODIFIED BY THE MEMBERS OF THE ESG GROUP.    *
!        THEREFORE PROF. BERZ SHOULD NOT BE HELD RESPONSIBLE FOR ANY BUGS.    *
!                                                                             *
!                  NEW RULES OF THE GAME (EXHAUSTIVE)                         *
!                 **********************************                          *
!                         THERE ARE NONE                                      *
!                                                                             *
!******************************************************************************
!
!
!     THIS FILE CONTAINS ROUTINES TO PERFORM DIFFERENTIAL ALGEBRA (DA)
!     AS AN OPTION, ALSO COMPONENTWISE ALGEBRA (CA) CAN BE PERFORMED.
!     A DESCRIPTION OF THE INTERNAL ARRAYS USED BY THE ROUTINES CAN
!     BE FOUND IN BLOCKDATA DABLD.
!
!
!     SHORT REFERENCE CHART
!     *********************
!
!     THE PARAMETERS USED BELOW HAVE THE FOLLOWING MEANING:
!
!     A,B:                NAME OF INPUT DA VECTORS   (INTEGER)
!     C:                  NAME OF OUTPUT DA VECTOR   (INTEGER)
!     X,Y:                NAME OF INPUT DA MATRIX    (INTEGER(...))
!     Z:                  NAME OF OUTPUT DA MATRIX   (INTEGER(...))
!
!     F:                  NAME OF A DA FUNCTION      (CHARACTER*4)
!     G:                  NAME OF EXTERNAL FUNCTION  (DOUBLE PRECISION)
!     JJ:                 ARRAY OF EXPONENTS         (INTEGER(20))
!     O:                  ORDER                      (INTEGER)
!     N:                  NUMBER OF VARIABLES        (INTEGER)
!     I,J,K:              INTEGER NUMBER             (INTEGER
!     R,RA,RB:            REAL NUMBERS               (DOUBLE PRECISION)
!     H:                  ARRAY OF LENGTH LH         (DOUBLE PRECISION)
!     U:                  OUTPUT UNIT NUMBER         (INTEGER)
!     T:                  COMMENT TEXT               (CHARACTER*10)
!
!
!               SUBROUTINES AND THEIR CALLING PARAMETERS
!               ****************************************
!
!     DAINI(O,N,U):       INITIALIZES CONTROL ARRAYS AND SETS MAX. ORDER O AND
!                         MAX. NUMBER OF VARIABLES N. MUST BE CALLED BEFORE ANY
!                         OTHER DA ROUTINE CAN BE USED.
!
!     DAALL(A,I,T,O,N):   ALLOCATES SPACE FOR I VECTORS A. T: CHARACTER NAME
!     DADAL(A,I):         DEALLOCATES THE I VECTORS A.
!!     DAVAR(A,R,I):       MAKES A INDEPENDENT VARIABLE # I WITH INITIAL VALUE R
!!     DACON(A,R):         SETS A TO CONSTANT R
!     DANOT(O):           SETS NEW TRUNCATION ORDER O FOR DA OPERATIONS
!     DAEPS(R):           SETS NEW ZERO TOLERANCE EPSILON
!
!!     DAPEK(A,JJ,R):      RETURNS COEF R OF MONOMIAL WITH EXPONENTS JJ OF A
!!     DAPOK(A,JJ,R):      SETS COEF OF MONOMIAL WITH EXPONENTS JJ OF A TO R
!
!!     DACOP(A,C):         PERFORMS C = A
!!     DAADD(A,B,C):       PERFORMS C = A + B
!!    DASUB(A,B,C):       PERFORMS C = A - B
!!     DAMUL(A,B,C):       PERFORMS C = A * B
!!     DADIV(A,B,C):       PERFORMS C = A / B
!!     DASQR(A,C):         PERFORMS C = A^2           (SQUARE OF A)
!
!!     DACAD(A,RA,C):      PERFORMS C = A + RA
!!     DACSU(A,RA,C):      PERFORMS C = A - RA
!!     DASUC(A,RA,C):      PERFORMS C = RA - A
!!     DACMU(A,RA,C):      PERFORMS C = A * RA
!!    DACDI(A,RA,C):      PERFORMS C = A / RA
!!     DADIC(A,RA,C):      PERFORMS C = RA / A
!!     DACMA(A,B,RB,C):    PERFORMS C = A + RB*B
!!DAMULIN(A,B,RA,C,D,RB,C):    PERFORMS C = A*B*RA + C*D*RB
!!     DALIN(A,RA,B,RB,C): PERFORMS C = A*RA + B*RB
!!     DAFUN(F,A,C):       PERFORMS C = F(A)          (DA FUNCTION)
!
!!     DAABS(A,R):         PERFORMS R = |A|           (NORM OF A)
!!     DACOM(A,B,R):       PERFORMS R = |A-B|         (NORM OF A-B)
!!     DAPOS(A,C):         PERFORMS C(I) = |A(I)|     (MAKE SIGNS POSITIVE)
!
!!     DACCT(X,I,Y,J,Z,K)  CONCATENATES Z = X o Y;   I,J,K: # OF VECTORS IN X,Y,
!!     DAINV(X,I,Z,K)      INVERTS Z = X^-1;           I,J: # OF VECTORS IN X,Y
!!     DAPIN(X,I,Z,K,JJ)   PARTIALLY INVERTS Z = X^-1; I,J: # OF VECTORS IN X,Y,
!                         JJ: ARRAY; NONZERO ENTRIES DENOTE TO BE INVERTED LINES
!
!!     DADER(I,A,C):       PERFORMS C = DA/DI (DERIV. WITH RESPECT TO VARIABLE I
!!     DAPOI(A,B,C,I):     PERFORMS C = [A,B] (POISSON BRACKET, 2*I: # PHASEVARS
!!     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
!     DAIMP(H,LH,A):      "IMPORTS" THE ARRAY H WITH LENGTH LH TO DA VAR A
!     DAEXP(H,LH,A):      "EXPORTS" THE DA VAR A TO ARRAY H WITH LENGTH LH
!     DAPRI(A,U):         PRINTS DA VECTOR A TO UNIT U
!     DAREA(A,U):         READS DA VECTOR A FROM UNIT U
!     DADEB(U,T,I):       DEBUGGER, DUMPS TO U. T: MEMO, I=0: RETURN, I=1:STOP
!!     DARAN(A,R,seed):         FILLS A WITH RANDOM NUMBERS. R: FILLFACTOR
!     DANUM(O,N,I):       COMPUTES NUMBER OF MONOMIALS IN N VAR THROUGH ORDER O
!
!
!     ADDITIONAL ROUTINES THE USER DOES NOT NEED TO CALL:
!
!     DAINF: RETURNS INFOS ABOUT A DA VECTOR PREVIOUSLY DECLARED
!     DAPAC: PACKS DA VECTORS
!     DACHK: CHECKS IF DA VECTORS HAVE COMPATIBLE ATTRIBUTES
!     DCODE: TRANSFORMS DIGITS IN A CERTAIN BASE TO A DECIMAL INTEGER
!     NCODE: EXTRACTS DIGITS IN A CERTAIN BASE FROM A DECIMAL INTEGER
!
!
!     FURTHER WISHES
!     **************
!
!     - CHECK DAREA AND DAPRI FOR CA VECTORS
!     - MAKE DAFUN USE DASQR
!
!
!      BLOCKDATA DABLD
!     ***************
!
!
!     PARAMETERS:
!
!     LDA: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
!     LST: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
!     LEA: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO,NV
!     LIA: DIMENSION OF IA1,IA2;            CAN BE INCREASED FOR LARGE NO,NV
!     LNO: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
!     LNV: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
!
!     ALL THE CHANGES IN THE VALUES OF PARAMETERS HAVE TO BE MADE BY GLOBAL
!     SUBSTITUTIONS IN ALL SUBROUTINES.
!
!     DANAME:   NAME OF DA VECTOR
!
!     CC:       STACK OF DOUBLE PRECISON COEFFICIENTS
!     I1:       FIRST CHARACTERISTIC INTEGER (CF DAINI)
!     I2:       SECOND CHARACTERISTIC INTEGER (CF DAINI)
!
!     IE1:      CHARACTERISTIC INTEGER 1 OF UNPACKED REPRESENTATION (CF DAINI)
!     IE2:      CHARACTERISTIC INTEGER 2 OF UNPACKED REPRESENTATION (CF DAINI)
!     IEO:      ORDER OF ENTRY IN UNPACKED REPRESENTATION
!     IA1:      REVERSE TO IE1 (CF DAINI)
!     IA2:      REVERSE TO IE2 (CF DAINI)
!
!     IDANO:    ORDER OF DA VECTOR; IN CA, NUMBER OF COMPONENTS
!     IDANV:    NUMBER OF VARIABLES; IF 0, INDICATES CA VECTOR
!     IDAPO:    FIRST ADDRESS IN STACK
!     IDALM:    NUMBER OF RESERVED STACK POSITIONS
!     IDALL:    NUMBER OF MOMENTARILY REQUIRED STACK POSITIONS
!
!     NDA:      NUMBER OF DA VECTORS MOMENTARILY DEFINED
!     NST:      NUMBER OF STACK POSITIONS MOMENTARILY ALLOCATED
!     NOMAX:    MAXIMUM REQUESTED ORDER  (CF DAINI)
!     NVMAX:    MAXIMUM REQUESTED NUMBER OF VARIABLES (CF DAINI)
!     NMMAX:    MAXIMUM NUMBER OF MONOMIALS FOR NOMAX, NVMAX (CF DAINI)
!     NOCUT:    MOMENTARY TRUNCATION ORDER
!     EPS:      TRUNCATION ACCURACY (CAN BE SET BY USER)
!     EPSMAC:   MANTISSA LENGTH OF MACHINE (PESSIMISTIC ESTIMATE)
!
!-----------------------------------------------------------------------------1
!

      module tpsa
      contains

      subroutine daini(no,nv,iunit) bind(C, name="daini_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) no, nv, iunit

      integer i,iall,ibase,ic1,ic2,icmax,io1,io2,iout,j,jd,jj,jjj,      &
     &jjjj,jl,js,k,n
      integer(8) nn
!     *****************************
!
!     THIS SUBROUTINE SETS UP THE MAJOR ORDERING AND ADDRESSING ARRAYS IN
!     COMMON BLOCK DAINI. IF IUNIT > 0, THE ARRAYS WILL BE PRINTED TO UNIT
!     NUMBER IUNIT. AN EXAMPLE FOR THE ARRAYS GENERATED BY DAINI CAN BE
!     FOUND AFTER THE ROUTINE.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
!      COMMON / DASCR /  IS(20), RS(20)                                        1
      integer idao,is,iscrri
      double precision rs
      common/dascr/is(100),rs(100),iscrri(100),idao
!-----------------------------------------------------------------------------2

      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

!
      character aa*10
      dimension n(lnv+1),k(0:lnv),j(lnv),jj(lnv)
!
      write(*, *) "daini:", no, nv
      if(eps.le.0.d0) eps=1.d-38
!      if(EPS.le.0.d0) eps=1.d-90
      epsmac=1.d-7
      if(nv.eq.0) return
      ndamaxi=0
!
      do i=1, lda
        allvec(i) = .false.
      enddo
      nhole=0
!*****************************************

!     INITIALIZING VARIABLES IN COMMON / DA /
!     ***************************************
!
      nda   = 0
      nst   = 0
      nomax = no
      nvmax = nv
      call danum(no,nv,nmmax)
      nocut = no
      lfi   = 0
!
      do i=0,lia
        ia1(i) = 0
        ia2(i) = 0
      enddo
!
      do i=1,100
        is(i) = 0
      enddo
!
      if(nv.gt.lnv.or.no.gt.lno) then
        write(6,*)'ERROR IN SUBROUTINE DAINI, NO, NV = ',no,nv
        call dadeb(31,'ERR DAINI ',1)
      endif
!
      ibase = no+1
      js    = nv/2
      if(float(ibase)**((nv+1)/2).gt.float(lia)) then
        write(6,*)'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR',      &
     &  ' LIA = ',lia
        call dadeb(31,'ERR DAINI ',1)
      endif
!
      icmax = 0
      nn    = 0
      k(0)  = 0
!
      do io2=0,no
!     ***************
!
        n(1)  = io2
        jl    = 0
        jd    = 1
!
  50    jl    = jl + jd
!

!old
!      IF(JL.EQ.0) THEN
!old
!     modified according to Wu Ying
!
        if(jl.le.0) then
!
          goto 100
        elseif(jd.eq.1) then
          j(jl) = 0
        else
          j(jl) = j(jl) + 1
        endif
!
        k(jl)    = k(jl-1)*ibase + j(jl)
        n(jl+1)  = n(jl) - j(jl)
!
        if(j(jl).gt.n(jl)) then
          jd    = -1
          goto 50
        elseif(jl.lt.js) then
          jd    = 1
          goto 50
        else
          j(jl) = n(jl)
          k(jl) = k(jl-1)*ibase + j(jl)

          ic2   = k(jl)
          icmax = max(icmax,ic2)
          k(jl) = 0
!
          ia2(ic2) = nn
!
          do io1=0,no-io2
!        ******************
!
            n(js+1) = io1
            jd      = 1
!
  70        jl      = jl + jd
!
            if(jl.eq.js) then
              goto 80
            elseif(jd.eq.1) then
              j(jl) = 0
            else
              j(jl) = j(jl) + 1
            endif
!
            k(jl)    = k(jl-1)*ibase + j(jl)
            n(jl+1)  = n(jl) - j(jl)
!
            if(j(jl).gt.n(jl)) then
              jd    = -1
              goto 70
            elseif(jl.lt.nv) then
              jd    = 1
              goto 70
            else
              jd    = -1
              j(jl) = n(jl)
              k(jl) = k(jl-1)*ibase + j(jl)
              ic1   = k(jl)
              icmax = max(icmax,ic1)
              nn = nn + 1
!
              ie1(nn) = ic1
              ie2(nn) = ic2
              i1 (nn) = ic1
              i2 (nn) = ic2
              if(ic2.eq.0) ia1(ic1) = nn
              ieo(nn) = io1 + io2
!
              goto 70
            endif
!
   80       continue
          enddo
!
          jd = -1
          goto 50
        endif
!
  100   continue
      enddo
!
      if(nn.gt.lea) then
        write(6,*)'ERROR IN DAINI, NN = ',nn,' EXCEEDS LEA'
        call dadeb(31,'ERR DAINI ',1)
      endif
!
!     ALLOCATING SCRATCH VARIABLES
!     ****************************
!
      iall = 0
      call daall1(iall,'$$UNPACK$$',nomax,nvmax)
!
      do i=0,nomax
        aa = '$$MUL   $$'
        write(aa(6:10),'(I5)') i
        iall = 0
!      CALL DAALL(IALL,1,AA,I,NVMAX)
        call daall1(iall,aa,nomax,nvmax)
      enddo
!
      idall(1) = nmmax
!
!     DOUBLE CHECKING ARRAYS IE1,IE2,IA1,IA2
!     **************************************
!
      do i=1,nmmax
!
        jjj = ia1(ie1(i)) + ia2(ie2(i))
        if(jjj.ne.i) then
          write(6,*)'ERROR IN DAINI IN ARRAYS IE1,IE2,IA1,IA2 AT I = ',i
          call dadeb(31,'ERR DAINI ',1)
        endif
!
      enddo
!
      if(iunit.eq.0) return
!
      write(6,*)'ARRAY SETUP DONE, BEGIN PRINTING'
!
      iout = 32
!      open(iout,file='DAINI.DAT',status='NEW')
      open(iout,file='DAINI.DAT',status='UNKNOWN')
!CRAY OPEN(IOUT,FILE='DAINI',STATUS='UNKNOWN',FORM='FORMATTED')          *CRAY
!CRAY REWIND IOUT                                                        *CRAY
!
      write(iout,'(/A/A/)') ' ARRAYS I1 THROUGH I20, IE1,IE2,IEO',      &
     &' **********************************'
      do i=1,nmmax
        call dancd(ie1(i),ie2(i),jj)
        write(iout,'(1X,I5,2X,4(5I2,1X),3I6)') i,(jj(jjjj),jjjj=1,lnv), &
     &  ie1(i),ie2(i),ieo(i)
      enddo
!
      write(iout,'(/A/A/)') ' ARRAYS IA1,IA2',' **************'
      do i=0,icmax
        write(iout,'(3I10)') i,ia1(i),ia2(i)
      enddo
!
      return
      end subroutine

      subroutine daexter() bind(C, name="daexter_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none

      integer i
!     *****************************
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

!
      do i=1, lda
        allvec(i)=.false.
      enddo

      return
      end subroutine

      subroutine dallsta(ldanow)
      implicit none
      integer i,ldanow
!     *****************************
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

!
      ldanow=0
      do i=1, lda
        if(allvec(i)) ldanow=ldanow+1
      enddo

      write(6,*) ' ALLOCATED ',ldanow

      return
      end subroutine

!
! EXAMPLE: ARRAYS I1 THROUGH I20, IE1,IE2,IEO (NOMAX=3,NVMAX=4)
! *************************************************************
!     I   I1               THROUGH               I20     IE1   IE2   IEO
!     1   0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      0     0     0
!     2   1 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      1     0     1
!     3   0 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      4     0     1
!     4   2 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      2     0     2
!     5   1 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      5     0     2
!     6   0 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      8     0     2
!     7   3 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      3     0     3
!     8   2 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      6     0     3
!     9   1 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      9     0     3
!    10   0 3 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0     12     0     3
!    11   0 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      0     1     1
!    12   1 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      1     1     2
!    13   0 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      4     1     2
!    14   2 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      2     1     3
!    15   1 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      5     1     3
!    16   0 2 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      8     1     3
!    17   0 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      0     4     1
!    18   1 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      1     4     2
!    19   0 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      4     4     2
!    20   2 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      2     4     3
!    21   1 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      5     4     3
!    22   0 2 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      8     4     3
!    23   0 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      0     2     2
!    24   1 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      1     2     3
!    25   0 1 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      4     2     3
!    26   0 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      0     5     2
!    27   1 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      1     5     3
!    28   0 1 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      4     5     3
!    29   0 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      0     8     2
!    30   1 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      1     8     3
!    31   0 1 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      4     8     3
!    32   0 0 0 0 0  0 0 0 0 0  3 0 0 0 0  0 0 0 0 0      0     3     3
!    33   0 0 0 0 0  0 0 0 0 0  2 1 0 0 0  0 0 0 0 0      0     6     3
!    34   0 0 0 0 0  0 0 0 0 0  1 2 0 0 0  0 0 0 0 0      0     9     3
!    35   0 0 0 0 0  0 0 0 0 0  0 3 0 0 0  0 0 0 0 0      0    12     3
!
!    ARRAYS IA1,IA2
!    **************
!    I        IA1       IA2
!    0         1         0   IE1,IE2 AND IA1,IA2 ALLOW THE EASY COMPUTATION
!    1         2        10   OF THE ADDRESS OF THE PRODUCT OF TWO MONOMIALS.
!    2         4        22   LET IX AND IY BE THE POSITIONS OF THE TWO
!    3         7        31   FACTORS. THEN THE POSITION IZ OF THE PRODUCT OF
!    4         3        16   THE TWO FACTORS IS GIVEN BY
!    5         5        25
!    6         8        32   IZ = IA1(IE1(IX)+IE1(IY)) + IA2(IE2(IX)+IE2(IY))
!    7         0         0
!    8         6        28
!    9         9        33   THE OTHER VARIABLES SET BY DAINI WOULD HAVE THE
!   10         0         0   VALUES
!   11         0         0
!   12        10        34   NOMAX = 3,  NVMAX = 4, NMMAX = 35
!

      subroutine daallno1(ic,ccc) bind(C, name="daallno1_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ic
      character(c_char) ccc(*)

      integer ind,l,ndanum,no,nv
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NOmax AND NUMBER OF VARIABLES NVmax
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      logical allvec(lda)
      integer nhole,j
      common /hole/nhole
      common /alloc/ allvec

!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      logical incnda
      character c*10
!
      no=nomax
      nv=nvmax
      ind = 1
      if(ic.gt.0.and.ic.le.nda) then
!     DANAME(IC) = C
!     IF(IDANO(IC).EQ.NO.AND.IDANV(IC).EQ.NV) THEN
      else
         if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
            write(6,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',     &
     &           no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
         endif
!     
         if(nhole.gt.0) then
            ind=nda
 20         if (allvec(ind)) then
               ind = ind - 1
               goto 20
            endif
            incnda = .false.
            nhole=nhole-1
         else
            incnda = .true.
            nda = nda + 1
            ind=nda
            if(nda.gt.lda) then
               write(6,*)'ERROR IN DAALL, MAX NUMBER OF DA VECTORS ',   &
     &              'EXHAUSTED'
               call dadeb(31,'ERR DAALL ',1)
            endif
         endif

         allvec(ind) = .true.
         ic = ind
!     
         if(nv.ne.0) then
            call danum(no,nv,ndanum)
         else
            ndanum = no
         endif

 !        c = ccc
         do j=1, 4
            c(j:j) = ccc(j)
         enddo
         write(c(5:10),'(I6)') ic
         c(5:10) = adjustl(c(5:10))
         daname(ind) = c

         if (incnda) then
            if(ind.gt.nomax+2) then
               idano(ind) = nomax
               idanv(ind) = nvmax
               idapo(ind) = nst + 1
               idalm(ind) = nmmax
               idall(ind) = 0
               nst = nst + nmmax
            else
               idano(ind) = no
               idanv(ind) = nv
               idapo(ind) = nst + 1
               idalm(ind) = ndanum
               idall(ind) = 0
               nst = nst + ndanum
            endif
         endif
!     
         if(nst.gt.lst) then
            x=-1.d0
            write(6,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(6,*) ' NST,LST '
            write(6,*)  nst,lst
            write(6,*) ' NDA,NDANUM,NDA*NDANUM '
            write(6,*)  nda,ndanum,nda*ndanum
!     X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
         endif
!     
         if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic)
            idall(ic) = idalm(ic)
         endif
      endif
!     
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end subroutine


      subroutine daallno(ic,l,ccc) bind(C, name="daallno_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ic(*), l
      character(c_char) ccc(*)

      integer i,j,ind,ndanum,no,nv
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NOmax AND NUMBER OF VARIABLES NVmax
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      logical incnda
      character c*10
!
      no=nomax
      nv=nvmax
      ind = 1
      do i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
            write(6,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',     &
     &      no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif
!
          if(nhole.gt.0) then
            ind=nda
20          if (allvec(ind)) then
              ind = ind - 1
              goto 20
            endif
            incnda = .false.
            nhole=nhole-1
          else
            incnda = .true.
            nda = nda + 1
            ind=nda
            if(nda.gt.lda) then
              write(6,*)'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHA',&
     &        'USTED'
              call dadeb(31,'ERR DAALL ',1)
            endif
          endif

          allvec(ind) = .true.
          ic(i) = ind
!
          if(nv.ne.0) then
            call danum(no,nv,ndanum)
          else
            ndanum = no
          endif

!          c = ccc
          do j=1, 4
             c(j:j) = ccc(j)
          enddo
          write(c(5:10),'(I6)') ic(i)
          c(5:10) = adjustl(c(5:10))
          daname(ind) = c

          if (incnda) then
            if(ind.gt.nomax+2) then
              idano(ind) = nomax
              idanv(ind) = nvmax
              idapo(ind) = nst + 1
              idalm(ind) = nmmax
              idall(ind) = 0
              nst = nst + nmmax
            else
              idano(ind) = no
              idanv(ind) = nv
              idapo(ind) = nst + 1
              idalm(ind) = ndanum
              idall(ind) = 0
              nst = nst + ndanum
            endif
          endif
!
          if(nst.gt.lst) then
            x=-1.d0
            write(6,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(6,*) ' NST,LST '
            write(6,*)  nst,lst
            write(6,*) ' NDA,NDANUM,NDA*NDANUM '
            write(6,*)  nda,ndanum,nda*ndanum
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif
!
          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
      enddo
!
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end subroutine


      subroutine daall1(ic,ccc,no,nv) bind(C, name="daall1_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ic, no, nv
      character(c_char) ccc(*)

      integer j,ind,l,ndanum
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NO AND NUMBER OF VARIABLES NV
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec
      integer ndat
      common /alloctot/ ndat

!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      logical incnda
!      character c*10,ccc*10
      character c*10
!
      ind = 1

      if(ic.gt.0.and.ic.le.nda) then
!     DANAME(IC) = C
!     IF(IDANO(IC).EQ.NO.AND.IDANV(IC).EQ.NV) THEN
      else
         if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
            write(6,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',     &
     &           no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
         endif
!     
         if(nhole.gt.0) then
            ind=nda
 20         if (allvec(ind)) then
               ind = ind - 1
               goto 20
            endif
            incnda = .false.
            nhole=nhole-1
         else
            incnda = .true.
            nda = nda + 1
            ind=nda
            if(nda.gt.lda) then
               write(6,*)'ERROR IN DAALL, MAX NUMBER OF DA VECTORS ',   &
     &              'EXHAUSTED'
               call dadeb(31,'ERR DAALL ',1)
            endif
         endif

         allvec(ind) = .true.

         ic = ind
!     
         if(nv.ne.0) then
            call danum(no,nv,ndanum)
         else
            ndanum = no
         endif

!     c = ccc
         do j=1, 4
            c(j:j) = ccc(j)
         enddo
         write(c(5:10),'(I6)') ic
         c(5:10) = adjustl(c(5:10))
         daname(ind) = c

         if (incnda) then
            if(ind.gt.nomax+2) then
               idano(ind) = nomax
               idanv(ind) = nvmax
               idapo(ind) = nst + 1
               idalm(ind) = nmmax
               idall(ind) = 0
               nst = nst + nmmax
            else
               idano(ind) = no
               idanv(ind) = nv
               idapo(ind) = nst + 1
               idalm(ind) = ndanum
               idall(ind) = 0
               nst = nst + ndanum
            endif
         endif
!     
         if(nst.gt.lst) then
            x=-1.d0
            write(6,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(6,*) ' NST,LST '
            write(6,*)  nst,lst
            write(6,*) ' NDA,NDANUM,NDA*NDANUM '
            write(6,*)  nda,ndanum,nda*ndanum
!     X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
         endif
!     
!     IF(NV.EQ.0) THEN
         if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic)
            idall(ic) = idalm(ic)
         endif
      endif
!     
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end subroutine
!
      subroutine daall(ic,l,ccc,no,nv) bind(C, name="daall_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ic(*), l, no, nv
      character(c_char) ccc(*)

      integer i,j,ind,ndanum
      double precision x
!     ********************************
!
!     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
!     ORDER NO AND NUMBER OF VARIABLES NV
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec
      integer ndat
      common /alloctot/ ndat

!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      logical incnda
!      character c*10,ccc*10
      character c*10
!
      ind = 1

      do i=1,l
        if(ic(i).gt.0.and.ic(i).le.nda) then
!         DANAME(IC(I)) = C
!         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
        else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
            write(6,*)'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',     &
     &      no,nv,' NOMAX, NVMAX = ',nomax,nvmax
            call dadeb(31,'ERR DAALL ',1)
          endif
!
          if(nhole.gt.0) then
            ind=nda
20          if (allvec(ind)) then
              ind = ind - 1
              goto 20
            endif
            incnda = .false.
            nhole=nhole-1
          else
            incnda = .true.
            nda = nda + 1
            ind=nda
            if(nda.gt.lda) then
              write(6,*)'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHA',&
     &        'USTED'
              call dadeb(31,'ERR DAALL ',1)
            endif
          endif

          allvec(ind) = .true.

          ic(i) = ind
!
          if(nv.ne.0) then
            call danum(no,nv,ndanum)
          else
            ndanum = no
          endif

!          c = ccc
          do j=1, 4
             c(j:j) = ccc(j)
          enddo
          write(c(5:10),'(I6)') ic(i)
          c(5:10) = adjustl(c(5:10))
          daname(ind) = c

          if (incnda) then
            if(ind.gt.nomax+2) then
              idano(ind) = nomax
              idanv(ind) = nvmax
              idapo(ind) = nst + 1
              idalm(ind) = nmmax
              idall(ind) = 0
              nst = nst + nmmax
            else
              idano(ind) = no
              idanv(ind) = nv
              idapo(ind) = nst + 1
              idalm(ind) = ndanum
              idall(ind) = 0
              nst = nst + ndanum
            endif
          endif
!
          if(nst.gt.lst) then
            x=-1.d0
            write(6,*)'ERROR IN DAALL, STACK EXHAUSTED '
            write(6,*) ' NST,LST '
            write(6,*)  nst,lst
            write(6,*) ' NDA,NDANUM,NDA*NDANUM '
            write(6,*)  nda,ndanum,nda*ndanum
!            X=DSQRT(X)
            call dadeb(31,'ERR DAALL ',1)
          endif
!
!          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.nomax.eq.1) then
            call daclr(ic(i))
            idall(ic(i)) = idalm(ic(i))
          endif
        endif
      enddo
!
      if(nda.gt.ndamaxi) ndamaxi=nda

      return
      end subroutine
!
      subroutine dadal1(idal) bind(C, name="dadal1_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) idal

!     ************************
!
!     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

      if(idal.le.nomax+2.or.idal.gt.nda) then
         write(6,*)'ERROR IN ROUTINE DADAL, IDAL,NDA = ',idal,nda
         call dadeb(31,'ERR DADAL ',1)
      endif
      if(idal.eq.nda) then
!     deallocate
         nst = idapo(nda) - 1
         nda = nda - 1
!     else
!     write(6,'(a10)')daname(i)
!     write(6,*)' etienne',idal,nda
!     write(6,*) sqrt(-1.d0)
      else
         nhole=nhole+1
      endif

      allvec(idal) = .false.

!     IDANO(IDAL) = 0
!     IDANV(IDAL) = 0
!     IDAPO(IDAL) = 0
!     IDALM(IDAL) = 0
      idall(idal) = 0

      idal = 0

      return
      end subroutine


      subroutine dadal(idal,l) bind(C, name="dadal_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) idal(*), l

      integer i
!     ************************
!
!     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
      logical allvec(lda)
      integer nhole
      common /hole/nhole
      common /alloc/ allvec

      do i=l,1,-1
        if(idal(i).le.nomax+2.or.idal(i).gt.nda) then
          write(6,*)'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda
          call dadeb(31,'ERR DADAL ',1)
        endif
        if(idal(i).eq.nda) then
!       deallocate
          nst = idapo(nda) - 1
          nda = nda - 1
!        else
!        write(6,'(a10)')daname(i)
!        write(6,*)' etienne',idal(i),nda
!        write(6,*) sqrt(-1.d0)
        else
          nhole=nhole+1
        endif

        allvec(idal(i)) = .false.

!        IDANO(IDAL(I)) = 0
!        IDANV(IDAL(I)) = 0
!        IDAPO(IDAL(I)) = 0
!        IDALM(IDAL(I)) = 0
        idall(idal(i)) = 0

        idal(i) = 0
      enddo

      return
      end subroutine


      subroutine davar(ina,ckon,i) bind(C, name="davar_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, i
      real(c_double) ckon

      integer ibase,ic1,ic2,illa,ilma,inoa,inva
      integer(8) ipoa
!     ****************************
!
!     THIS SUBROUTINE DECLARES THE DA VECTOR
!     AS THE INDEPENDENT VARIABLE NUMBER I.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      if(i.gt.inva) then
        write(6,*)'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
        call dadeb(31,'ERR DAVAR ',1)
      endif
!
      if(nomax.eq.1) then
        if(i.gt.inva) then
          print*,'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
!           CALL DADEB(31,'ERR DAVAR3',1)
        endif
        call daclr(ina)
        cc(ipoa) = ckon
        cc(ipoa+i) = 1d0
        return
      endif

      ibase = nomax+1
!
      if(i.gt.(nvmax+1)/2) then
        ic1 = 0
        ic2 = ibase**(i-(nvmax+1)/2-1)
      else
        ic1 = ibase**(i-1)
        ic2 = 0
      endif
!
      if(dabs(ckon).gt.eps) then
        idall(ina) = 2
        cc(ipoa) = ckon
        i1(ipoa) = 0
        i2(ipoa) = 0
!
        cc(ipoa+1) = 1.d0
        i1(ipoa+1) = ic1
        i2(ipoa+1) = ic2
      else
        idall(ina) = 1
        cc(ipoa) = 1.d0
        i1(ipoa) = ic1
        i2(ipoa) = ic2
      endif
!
      return
      end subroutine
!
      subroutine dacon(ina,ckon) bind(C, name="dacon_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina
      real(c_double) ckon

      integer illa,ilma,inoa,inva
      integer(8) ipoa
!     **************************
!
!     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      if(nomax.eq.1) then
        call daclr(ina)
        cc(ipoa) = ckon
        return
      endif

      idall(ina) = 1
      cc(ipoa) = ckon
      i1(ipoa) = 0
      i2(ipoa) = 0
      if(dabs(ckon).lt.eps) idall(ina) = 0
!
      return
      end subroutine
!
      subroutine danot(not) bind(C, name="danot_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) not

!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      if(not.gt.nomax) then
        write(6,*)'ERROR, NOCUT = ',not,' EXCEEDS NOMAX = ',nomax
        call dadeb(31,'ERR DANOT ',1)
      endif
!
      nocut = not
!
      return
      end subroutine

      function getno() bind(C, name="getno_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) getno

      include "TPSALib_prm.f"
      getno = nocut
      return
      end function

      subroutine getdanot(not)
      implicit none
      integer not
!     *********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      if(not.gt.nomax) then
        write(6,*)'ERROR, NOCUT = ',not,' EXCEEDS NOMAX = ',nomax
        call dadeb(31,'ERR DANOT ',1)
      endif
!
      not=nocut
!
      return
      end subroutine

      subroutine daeps(deps) bind(C, name="daeps_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      real(c_double) deps

!     **********************
!
!     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      if(deps.ge.0.d0) then
        eps = deps
      else
        deps=eps
      endif
!
      return
      end subroutine
!
      subroutine dapek(ina,jj,cjj) bind(C, name="dapek_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      include "TPSALib_prm.f"
      integer(c_long) ina, jj(lnv)
      real(c_double) cjj

      integer ibase,ic,ic1,ic2,icu,icz,ii1,ikk,illa,ilma,               &
     &inoa,inva,jj1,mchk
      integer(8) ipoa, iu, iz, i, ipek
!     ****************************
!
!     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ AND RETURNS IT IN CJJ
!
!-----------------------------------------------------------------------------1
!
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      jj1 = 1
      if(inva.eq.0.or.nomax.eq.1) then
        if(inva.ne.0.and.nomax.eq.1) then
          if(illa.ge.2) then
            do i=1,illa - 1
              if(jj(i).eq.1) jj1 = i + 1
            enddo
          else
            jj1 = jj(1) + 1
          endif
        else
          jj1 = jj(1)
        endif
        if(jj1.lt.1.or.jj1.gt.illa) then
          print*,'ERROR IN DAPEK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
!           CALL DADEB(31,'ERR DAPEK1',1)
        endif
        ipek = ipoa + jj1 - 1
        cjj = cc(ipek)
        return
      endif

      ii1 = (nvmax+1)/2
      ibase = nomax+1
!
!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************
!
      call dadcd(jj,ic1,ic2)
!
!ETIENNE
      if(ic1.gt.lia.or.ic2.gt.lia) then
        write(6,*) 'DISASTER IN DAPEK, INA= ',ina
        write(32,*) 'DISASTER IN DAPEK, INA= ',ina
        write(32,*) ic1,ic2
        write(32,*) (jj(ikk),ikk=1,lnv)
      endif
!ETIENNE
      ic = ia1(ic1) + ia2(ic2)
!
!     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
!     ****************************************************************
!
!      IF(ICO.GT.INOA.OR.ICV.GT.INVA.OR.ICO.GT.NOCUT) THEN
!         CJJ = 0
!         RETURN
!      ENDIF
!
!     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
!     *************************************************************
!
      iu = ipoa
      iz = ipoa + illa - 1
      icu = ia1(i1(iu))+ia2(i2(iu))
      icz = ia1(i1(iz))+ia2(i2(iz))
!
      if(illa.eq.0) then
        cjj = 0
        return
      elseif(ic.eq.icu) then
        cjj = cc(iu)
        return
      elseif(ic.eq.icz) then
        cjj = cc(iz)
        return
      elseif(ic.lt.icu.or.ic.gt.icz) then
        cjj = 0
        return
      endif
!
!     SEARCHING PROPER MONOMIAL
!     *************************
!
 10   continue
      if(iz-iu.le.1) then
        cjj = 0
        return
      endif
      i = (iu+iz)/2
!
!     if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
      mchk=ia1(i1(i))+ia2(i2(i)) - ic
      if(mchk.lt.0) goto 20
      if(mchk.eq.0) goto 30
      if(mchk.gt.0) goto 40
 20   iu = i
      goto 10
 30   cjj = cc(i)
      return
 40   iz = i
      goto 10
!
      end subroutine
!
      subroutine dapok(ina,jj,cjj) bind(C, name="dapok_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      include "TPSALib_prm.f"
      integer(c_long) ina, jj(lnv)
      real(c_double) cjj

      integer ic,ic1,ic2,icu,icz,illa,ilma,inoa,inva,                   &
     &jj1,mchk
      integer(8) i, ii, ipoa, iu, iz, ipok
!     ****************************
!
!     THIS SUBROUTINE SETS THE COEFFICIENT OF THE ARRAY
!     OF EXPONENTS JJ TO THE VALUE CJJ
!
!-----------------------------------------------------------------------------1
!
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
!
      jj1 = 1
      if(inva.eq.0.or.nomax.eq.1) then
        if(inva.ne.0.and.nomax.eq.1) then
          if(illa.ge.2) then
            do i=1,illa - 1
              if(jj(i).eq.1) jj1 = i + 1
            enddo
          else
            jj1 = jj(1) + 1
          endif
        else
          jj1 = jj(1)
        endif
        if(jj1.lt.1.or.jj1.gt.illa) then
          print*,'ERROR IN DAPOK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
!           CALL DADEB(31,'ERR DAPOK1',1)
        endif
        ipok = ipoa + jj1 - 1
        cc(ipok) = cjj
        return
      endif

!     DETERMINE INDEX TO BE SEARCHED FOR
!     **********************************
!
      call dadcd(jj,ic1,ic2)
!
      ic = ia1(ic1) + ia2(ic2)
!
!     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
!     ****************************************************************
!
!      IF(ICO.GT.INOA.OR.ICV.GT.INVA) THEN
!         write(6,*)'ERROR IN DAPOK, MONOMIAL NOT ALLOWED FOR ',A
!         CALL DADEB(31,'ERR DAPOK ',1)
!      ENDIF
!      IF(ICO.GT.NOCUT) RETURN
!
      iu = ipoa
      iz = ipoa + illa - 1
!
!     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
!     *************************************************************
!
      icu = ia1(i1(iu))+ia2(i2(iu))
      icz = ia1(i1(iz))+ia2(i2(iz))
      if(illa.eq.0) then
        i = ipoa
        goto 100
      elseif(ic.eq.icu) then
        cc(iu) = cjj
        i = iu
        goto 200
      elseif(ic.eq.icz) then
        cc(iz) = cjj
        i = iz
        goto 200
      elseif(ic.lt.icu) then
        i = iu
        goto 100
      elseif(ic.gt.icz) then
        i = iz + 1
        goto 100
      endif
!
!
!     SEARCHING PLACE TO POKE INTO OR BEFORE WHICH TO POKE
!     ****************************************************
!
      iu = ipoa
      iz = ipoa + illa
!
 10   continue
      if(iz-iu.le.1) then
        i = iz
        goto 100
      endif
      i = (iu+iz)/2
!
!      if(ia1(i1(i))+ia2(i2(i)) - ic) 20,30,40
      mchk=ia1(i1(i))+ia2(i2(i)) - ic
      if(mchk.lt.0) goto 20
      if(mchk.eq.0) goto 30
      if(mchk.gt.0) goto 40
 20   iu = i
      goto 10
 30   cc(i) = cjj
      goto 200
 40   iz = i
      goto 10
!
!     INSERTING THE MONOMIAL, MOVING THE REST
!     ***************************************
!
 100  continue
!
      if(dabs(cjj).lt.eps) return
!
      do ii=ipoa+illa,i+1,-1
        cc(ii) = cc(ii-1)
        i1(ii) = i1(ii-1)
        i2(ii) = i2(ii-1)
      enddo
!
      cc(i) = cjj
      i1(i) = ic1
      i2(i) = ic2
!
      idall(ina) = illa + 1
      if(idall(ina).gt.idalm(ina)) then
        write(6,*)'ERROR IN DAPOK '
        call dadeb(31,'ERR DAPOK ',1)
      endif
!
      return
!
!     CASE OF CJJ = 0 WHICH MEANS MOVING THE REST
!     *********************************************
!
 200  continue
      if(dabs(cjj).lt.eps) then
        do ii=i,ipoa+illa-2
          cc(ii) = cc(ii+1)
          i1(ii) = i1(ii+1)
          i2(ii) = i2(ii+1)
        enddo
        idall(ina) = illa - 1
      endif
      return
!
      end subroutine
!
      subroutine daclr(inc) bind(C, name="daclr_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) inc

      integer illc,ilmc,inoc,invc
      integer(8) i, ipoc
!     *********************
!
!     THIS SUBROUTINE SETS ALL THE STACK SPACE RESERVED FOR VARIABLE
!     C TO ZERO
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
      do i=ipoc,ipoc+ilmc-1
!
        cc(i) = 0.d0
!
      enddo
!
      return
      end subroutine
!
      subroutine dacop(ina,inb) bind(C, name="dacop_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina,inb

      integer iif,illa,illb,ilma,ilmb,inoa,inob,inva,invb
      integer(8) ia, ib, ipoa, ipob
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      ipob = idapo(inb)
      ipoa = idapo(ina)
      illa = idall(ina)
      ib = ipob - 1
!
!      iif = 0
!      if(nomax.eq.1.or.inva.eq.0) iif = 1

      do ia = ipoa,ipoa+illa-1
!
        if(nomax.gt.1) then
          if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
        endif
        ib = ib + 1
        cc(ib) = cc(ia)
        i1(ib) = i1(ia)
        i2(ib) = i2(ia)
!
 100    continue
      enddo
!
      idall(inb) = ib - ipob + 1
!      if(idall(inb).gt.idalm(inb)) then
!         write(6,*)'ERROR IN DACOP'
!         call dadeb(31,'ERR DACOP ',1)
!      endif
!
      return
      end subroutine

!
      subroutine datrashn(idif,ina,inbb)
      implicit none
      integer i,ia,idif,illa,ilma,ina,inb,inbb,inoa,inva
      integer(8) ipoa
      double precision rr
!     *************************
!
!     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      integer jd(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      inb=0

      if(inbb.eq.ina) then
        call daall1(inb,'$$DAADD $$',inoa,inva)
      else
        inb=inbb
      endif

      call daclr(inb)
!
!
      do ia = ipoa,ipoa+illa-1
!
        if(nomax.ne.1) then
          call dancd(i1(ia),i2(ia),jd)
        else
          do i=1,lnv
            jd(i)=0
          enddo
          if(ia.ne.ipoa) then
            jd(ia-ipoa+1)=1
          endif
        endif

        call dapek(ina,jd,rr)
        jd(idif)=0
        if(dabs(rr).gt.0.d0) call dapok(inb,jd,rr)
!
      enddo
!
!
      if(inbb.eq.ina) then
        call dacop(inb,inbb)
        call dadal1(inb)
      endif

      return
      end subroutine
!
      subroutine daadd(ina,inb,inc) bind(C, name="daadd_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb, inc

      integer i,ic,ic1,ic2,icu,icz,ii,illa,ilma,inoa,inva,              &
     &ipok,iu,iz,jj,jj1
      integer idaadd,illc,ilmc,inoc,invc
      integer(8) ipoc, ipoa, ipob

      include "TPSALib_prm.f"

!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA ADDITION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
        ipob = idapo(inb)
!         minv = min(inva,invb,invc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   + cc(ipob+i)
        enddo
        return
      endif

      if(ina.ne.inc.and.inb.ne.inc) then
        call dalin(ina,+1.d0,inb,+1.d0,inc)
      else
        idaadd = 0
        call daall1(idaadd,'$$DAADD $$',nomax,nvmax)
        call dalin(ina,+1.d0,inb,+1.d0,idaadd)
        call dacop(idaadd,inc)
        call dadal1(idaadd)
      endif
!
      return
      end subroutine
!
      subroutine dasub(ina,inb,inc) bind(C, name="dasub_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb, inc

      integer idasub
      integer i,ic,ic1,ic2,icu,icz,ii,illa,ilma,inoa,inva,              &
     &ipoa,ipok,iu,iz,jj,jj1
      integer idaadd,illc,ilmc,inoc,invc
      integer(8) ipoc

      include "TPSALib_prm.f"

      integer ipob
!     THIS SUBROUTINE PERFORMS A DA SUBTRACTION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
        ipob = idapo(inb)
!         minv = min(inva,invb,invc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   - cc(ipob+i)
        enddo
        return
      endif
      if(ina.ne.inc.and.inb.ne.inc) then
        call dalin(ina,+1.d0,inb,-1.d0,inc)
      else
        idasub = -1
!         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        call daall1(idasub,'$$DASUB $$',nomax,nvmax)
        call dalin(ina,+1.d0,inb,-1.d0,idasub)
        call dacop(idasub,inc)
        call dadal1(idasub)
      endif
!
      return
      end subroutine
!
      subroutine damulin(ina,inb,coe1,inc,ind,coe2,ine)
      implicit none
      integer ina,inb,inc,incc,ind,ine,inoc,invc
      double precision coe1,coe2
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!

      call daall1(incc,'$$DAJUNK$$',inoc,invc)
      call damul(ina,inb,incc)
      call damul(inc,ind,ine)
      call dalin(incc,coe1,ine,coe2,ine )
      call dadal1(incc)

      return
      end subroutine

!
!
! ANFANG UNTERPROGRAMM
      subroutine daexx(ina,inb,inc)
      implicit none
      integer illc,ilmc,ina,inaa,inb,inbb,inc,inoc,invc
      integer(8) ipoc

      include "TPSALib_prm.f"

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      write(6,*) "daexx"
      if(ina.ne.inc.and.inb.ne.inc) then
        call daexxt(ina,inb,inc)
      else
        inaa = 0
        inbb = 0
!         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        call daall1(inaa,'$$DAADD $$',nomax,nvmax)
        call daall1(inbb,'$$DAADD $$',nomax,nvmax)
        call dacop(ina,inaa)
        call dacop(inb,inbb)
        call daexxt(inaa,inbb,inc)
        call dadal1(inaa)
        call dadal1(inbb)
      endif

      return
      end subroutine

! ANFANG UNTERPROGRAMM
      subroutine daexxt(ina,inb,inc)
      implicit none
      integer idaexx,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,    &
     &inob,inoc,inva,invb,invc
      integer(8) ipoa,ipob,ipoc
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INA WITH INB
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
      idaexx = 0
      call daall1(idaexx,'$$DAEXX $$',nomax,nvmax)
      call dafun('LOG   ',ina,inc)
      call damul(inb,inc,idaexx)
      call dafun('EXP   ',idaexx,inc)
      call dadal1(idaexx)
!
      return
      end subroutine

! ANFANG UNTERPROGRAMM
      subroutine dacex(ina,ckon,inb)
      implicit none
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc
      integer(8) ipoc
      double precision ckon

      include "TPSALib_prm.f"

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      write(6,*) "dacex"
      if(ina.eq.inb) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call dacext(ina,ckon,incc)
        call dacop(incc,inb)
        call dadal1(incc)
      else
        call dacext(ina,ckon,inb)
      endif

      return
      end subroutine
      subroutine dacext(ina,ckon,inb)
      implicit none
      integer idacex,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8) ipoa,ipob
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES THE CONSTANT CKON WITH INA
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      if(ckon.le.0) then
        print*,'ERROR IN DACEX, CKON NOT POSITIVE'
!        CALL DADEB(31,'ERR DACEX1',1)
      endif
!
      idacex = 0
      call daall1(idacex,'$$DACEX $$',nomax,nvmax)
      ckon = dlog(ckon)
      call dacmu(ina,ckon,idacex)
      call dafun('EXP   ',idacex,inb)
      call dadal1(idacex)
!
      return
      end subroutine

      subroutine daexc(ina,ckon,inb)
      implicit none
      integer illc,ilmc,ina,inb,inc,incc,inoc,invc
      integer(8) ipoc
      double precision ckon

      include "TPSALib_prm.f"

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
!        write(6,*) "daexc"

      if(ina.eq.inb) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call daexct(ina,ckon,incc)
        call dacop(incc,inb)
        call dadal1(incc)
      else
        call daexct(ina,ckon,inb)
      endif

      return
      end subroutine

      subroutine daexct(ina,ckon,inb)
      implicit none
      integer i,ic,idaexc,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,   &
     &invb
      integer(8) ipoa, ipob
      double precision ckon,xic

      include "TPSALib_prm.f"

!     ******************************
!
!     THIS SUBROUTINE EXPONENTIATES INE WITH THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      idaexc = 0
      call daall1(idaexc,'$$DAEXC $$',nomax,nvmax)

      if(ckon.lt.0.d0) then
        call dafun('LOG   ',ina,inb)
        call dacmu(inb,ckon,idaexc)
        call dafun('EXP   ',idaexc,inb)
      else
        xic=dabs(ckon-dble(idint(ckon)))
        if(xic.gt.eps) then
          call dafun('LOG   ',ina,inb)
          call dacmu(inb,ckon,idaexc)
          call dafun('EXP   ',idaexc,inb)
        else
          ic=idint(ckon)
          call dacon(idaexc,1.d0)
          do i=1,ic
            call damul(idaexc,ina,idaexc)
          enddo
          call dacop(idaexc,inb)
        endif
      endif
      call dadal1(idaexc)
!
      return
      end subroutine

      subroutine damul(ina,inb,inc) bind(C, name="damul_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb, inc

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoc, ipoa, ipob
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1

      include "TPSALib_prm.f"

      double precision ccipoa,ccipob
      integer i
!
      if(nomax.eq.1) then
        ipoa=idapo(ina)
        ipob=idapo(inb)
        ipoc=idapo(inc)
!         minv = min(inva,invb,invc)
        ccipoa = cc(ipoa)
        ccipob = cc(ipob)
        cc(ipoc) = ccipoa*ccipob
        do i=1,nvmax
          cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
        enddo
!         do 30 i=ipoc+minv+1,ipoc+invc
!  30     cc(i) = 0d0
        return
      endif

      if(ina.eq.inc.or.inb.eq.inc) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call damult(ina,inb,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call damult(ina,inb,inc)
      endif

      return
      end subroutine

      subroutine damult(ina,inb,inc)
      implicit none
      integer i,i1ia,i2ia,illa,illb,illc,ilma,ilmb,ilmc,ina,            &
     &inb,inc,inoa,inoc,inva,invb,invc,ioffb,ipno,                      &
     &ipos,minv,noff,noib,nom
      integer(8) inob, ia, ib, ic, ipoa, ipob, ipoc
      double precision ccia,ccipoa,ccipob
!     *****************************
!
!     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
!     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension ipno(0:lno),noff(0:lno)
!
!
!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
!
!     CASE OF FIRST ORDER ONLY
!     ************************


      if(nomax.eq.1) then
        ipoa=idapo(ina)
        ipob=idapo(inb)
        ipoc=idapo(inc)
!         minv = min(inva,invb,invc)
        ccipoa = cc(ipoa)
        ccipob = cc(ipob)
        cc(ipoc) = ccipoa*ccipob
        do i=1,nvmax
          cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
        enddo
!         do 30 i=ipoc+minv+1,ipoc+invc
!  30     cc(i) = 0d0
        return
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)

!     GENERAL CASE
!     ************
!
      do i=0,nomax
        noff(i) = idapo(i+2)
        ipno(i) = 0
      enddo
!
      call daclr(1)
!
!     RE-SORTING THE VECTOR B INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************
!
      do ib=ipob,ipob+illb-1
!
        noib = ieo(ia1(i1(ib))+ia2(i2(ib)))
        ipos = ipno(noib) + 1
        ipno(noib) = ipos
        inob = noff(noib) + ipos
!
        cc(inob) = cc(ib)
        i1(inob) = i1(ib)
        i2(inob) = i2(ib)
!
      enddo
!
      do i=0,nomax
        idall(i+2) = ipno(i)
      enddo
!
!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************
!
      nom = min(nocut,inoc)
!
      do ia=ipoa,ipoa+illa-1
!
        i1ia = i1(ia)
        i2ia = i2(ia)
        ccia = cc(ia)
!
        do noib = 0,nom-ieo(ia1(i1(ia))+ia2(i2(ia)))
!
          ioffb = noff(noib)
!
          do ib = ioffb+1,ioffb+ipno(noib)
!
            ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
            cc(ic) = cc(ic) + ccia*cc(ib)
!
          enddo
        enddo
      enddo
!
      call dapac(inc)
!
      return
      end subroutine
!
      subroutine dadiv(ina,inb,inc) bind(C, name="dadiv_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb, inc

      integer idadiv
      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoa, ipob, ipoc
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      double precision ck,ck1
      integer i
!
!     THIS SUBROUTINE PERFORMS A DA DIVISION OF THE DA VECTORS A AND B.
!     THE RESULT IS STORED IN C.
!

      if(nomax.eq.1) then
!         minv = min(inva,invb)
        ipoa = idapo(ina)
        ipob = idapo(inb)
        ipoc = idapo(inc)
        ck=1.d0/cc(ipob)
        ck1=cc(ipoa)*ck
        do i=1,nvmax
          cc(ipoc+i) = (cc(ipoa+i)-cc(ipob+i)*ck1)*ck
        enddo
        cc(ipoc)=ck1
        return
      endif

      idadiv = 0
!      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall1(idadiv,'$$DADIV $$',nomax,nvmax)
      call dafun('INV ',inb,idadiv)
      call damul(ina,idadiv,inc)
      call dadal1(idadiv)
!
      return
      end subroutine

!
      subroutine dasqr(ina,inc) bind(C, name="dasqr_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inc

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoa, ipoc
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer i
      double precision ccipoa
!
!     CASE OF FIRST ORDER ONLY
!     ************************
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
!         minv = min(inva,invc)
        ccipoa = cc(ipoa)
        cc(ipoc) = ccipoa*ccipoa
        do i=1,nvmax
          cc(ipoc+i) = 2.d0*ccipoa*cc(ipoa+i)
        enddo
        return
      endif

      if(ina.eq.inc) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call dasqrt(ina,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dasqrt(ina,inc)
      endif

      return
      end subroutine

      subroutine dasqrt(ina,inc)
      implicit none
      integer i,i1ia,i2ia,ib1,illa,illc,ilma,ilmc,ina,inc,              &
     &inoc,inva,invc,ioffa,ioffb,ipno,ipos,                             &
     &minv,noff,noia,noib,nom
      integer(8) ia, ib, ic, inoa, ipoa, ipoc
      double precision ccia,ccipoa
!     *************************
!
!     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension ipno(0:lno),noff(0:lno)
!
!
!      CALL DACHK(INA,INOA,INVA,'          ',-1,-1,INC,INOC,INVC)
!
!
!     CASE OF FIRST ORDER ONLY
!     ************************
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
!         minv = min(inva,invc)
        ccipoa = cc(ipoa)
        cc(ipoc) = ccipoa*ccipoa
        do i=1,nvmax
          cc(ipoc+i) = 2d0*ccipoa*cc(ipoa+i)
        enddo
        return
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!     GENERAL CASE
!     ************
!
      do i=0,nomax
        noff(i) = idapo(i+2)
        ipno(i) = 0
      enddo
!
      call daclr(1)
!
!     RESORTING THE VECTOR A INTO PIECES THAT ARE OF ONLY ONE ORDER
!     *************************************************************
!
      do ia=ipoa,ipoa+illa-1
!
        noia = ieo(ia1(i1(ia))+ia2(i2(ia)))
        ipos = ipno(noia) + 1
        ipno(noia) = ipos
        inoa = noff(noia) + ipos
!
        cc(inoa) = cc(ia)
        i1(inoa) = i1(ia)
        i2(inoa) = i2(ia)
!
      enddo
!
      do i=0,nomax
        idall(i+2) = ipno(i)
      enddo
!
!     PERFORMING ACTUAL MULTIPLICATION
!     ********************************
!
      nom = min(nocut,inoc)
!
      do noia = 0,nom/2
!
        ioffa = noff(noia)
!
        do ia=ioffa+1,ioffa+ipno(noia)
!
          i1ia = i1(ia)
          i2ia = i2(ia)
          ccia = cc(ia)
!
          ic = ia2(i2ia+i2ia) + ia1(i1ia+i1ia)
          cc(ic) = cc(ic) + ccia*ccia
          ccia = ccia + ccia
!
          do noib = noia,nom-noia
!
            ioffb = noff(noib)
            if(noib.eq.noia) then
              ib1 = ia + 1
            else
              ib1 = ioffb + 1
            endif
!
            do ib = ib1,ioffb+ipno(noib)
!
              ic = ia2(i2ia+i2(ib)) + ia1(i1ia + i1(ib))
              cc(ic) = cc(ic) + ccia*cc(ib)
!
            enddo
          enddo
        enddo
      enddo
!
      call dapac(inc)
!
      return
      end subroutine
!
      subroutine dacad(ina,ckon,inb) bind(C, name="dacad_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb
      real(c_double) ckon

      integer illa,illb,ilma,ilmb,inoa,inob,inva,invb,ipoa,ipob
      double precision const
!     ******************************
!
!     THIS SUBROUTINE ADDS THE CONSTANT CKON TO THE VECTOR A
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer jj(lnv)
      data jj / lnv*0 /
!
!
!
      call dacop(ina,inb)
      if(nomax.eq.1) then
        cc(idapo(inb)) = cc(idapo(inb)) + ckon
        return
      endif
!
      call dapek(inb,jj,const)
      call dapok(inb,jj,const+ckon)
!
      return
      end subroutine
!
      subroutine dacsu(ina,ckon,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,ipob
      double precision ckon,const
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE CONSTANT CKON FROM THE VECTOR A
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer jj(lnv)
      data jj / lnv*0 /
!
!
      call dacop(ina,inb)
!
      if(nomax.eq.1) then
        cc(idapo(inb)) = cc(idapo(inb)) - ckon
        return
      endif
!
      call dapek(inb,jj,const)
      call dapok(inb,jj,const-ckon)
!
      return
      end subroutine
!
      subroutine dasuc(ina,ckon,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8) ipoa, ipob
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE SUBTRACTS THE VECTOR INA FROM THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer i
!
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      ipob=idapo(inb)
      ipoa=idapo(ina)
      if(nomax.eq.1) then
        cc(ipob) = ckon - cc(ipoa)
        do i=1,nvmax
          cc(ipob+i) =-cc(ipoa+i)
        enddo
        return
      endif

      call dacsu(ina,ckon,inb)
      call dacmu(inb,-1.d0,inb)
!
      return
      end subroutine
!
      subroutine dacmu(ina,ckon,inc) bind(C, name="dacmu_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inc
      real(c_double) ckon

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoa, ipoc
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer i
!
      if(nomax.eq.1) then
!         minv = min(inva,invb)
        ipoa = idapo(ina)
        ipoc = idapo(inc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * ckon
        enddo
        return
      endif

      if(ina.eq.inc) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daallno1(incc,'$$DAJUNK$$')
        call dacmut(ina,ckon,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dacmut(ina,ckon,inc)
      endif

      return
      end subroutine

      subroutine dacmut(ina,ckon,inb)
      implicit none
      integer i,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,        &
     &minv
      integer(8) ia, ib, ipoa, ipob
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      if(nomax.eq.1) then
!         minv = min(inva,invb)
        ipoa = idapo(ina)
        ipob = idapo(inb)
        do i=0,nvmax
          cc(ipob+i) = cc(ipoa+i) * ckon
        enddo
        return
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      if(dabs(ckon).lt.eps) then
        idall(inb) = 0
        return
      endif
!
      ib = ipob - 1
!
      do ia=ipoa,ipoa+illa-1
!
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
        ib = ib + 1
        cc(ib) = cc(ia)*ckon
        i1(ib) = i1(ia)
        i2(ib) = i2(ia)
!
 100    continue
      enddo
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
        write(6,*)'ERROR IN DACMU '
        call dadeb(31,'ERR DACMU ',1)
      endif
!
      return
      end subroutine
!
      subroutine dacdi(ina,ckon,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8) ipoa,ipob
      double precision ckon
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE VECTOR INA BY THE CONSTANT CKON
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer i
!
!      if(ckon.eq.0.d0) then
!         write(6,*)'ERROR IN DACDI, CKON IS ZERO'
!         call dadeb(31,'ERR DACDI ',1)
!      endif
!
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
      if(nomax.eq.1) then
!         minv = min(inva,invb)
        ipoa = idapo(ina)
        ipob = idapo(inb)
        do i=0,nvmax
          cc(ipob+i) = cc(ipoa+i)/ ckon
        enddo
        return
      endif

      call dacmu(ina,1.d0/ckon,inb)
!
      return
      end subroutine
!
!
      subroutine dadic(ina,ckon,inc)
      implicit none
      integer idadic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc
      integer(8) ipoa, ipoc
      double precision ckon,zero
!     ******************************
!
!     THIS SUBROUTINE DIVIDES THE CONSTANT CKON BY THE VECTOR INA
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      double precision ck
      integer i
!
!
!
      zero=0.d0
      if(nomax.eq.1) then
!         minv = min(inva,invb)
        ipoa = idapo(ina)
        ipoc = idapo(inc)
        cc(ipoc)=ckon/cc(ipoa)
        ck=cc(ipoc)/cc(ipoa)
        do i=1,nvmax
          cc(ipoc+i) = -cc(ipoa+i)* ck
        enddo
        return
      endif


      if(abs(ckon).lt.eps) then
        call dacon(inc,zero)
        return
      endif

      idadic = 0
      call daall1(idadic,'$$DADIC $$',nomax,nvmax)

      if(ckon.eq.0.d0) then
        write(6,*)'ERROR IN DACDI and DADIC, CKON IS ZERO'
        call dadeb(31,'ERR DACDI ',1)
      endif
      call dacdi(ina,ckon,idadic)
      call dafun('INV ',idadic,inc)
      call dadal1(idadic)
!
      return
      end subroutine
!
      subroutine dacma(ina,inb,bfac,inc)
      implicit none
      integer idacma,illc,ilmc,ina,inb,inc,inoc,invc
      integer(8) ipoc, ipob, ipoa
      double precision bfac
!     **********************************
!
!     THIS SUBROUTINE PERFORMS THE OPERATIONS C = A + B*BFAC, WHERE A,B,C ARE
!     DA VECTORS AND BFAC IS A DOUBLE PRECISION. A AND C CAN BE IDENTICAL.
!     CAN LATER BE REPLACED BY SOMETHING LIKE DAADD WITH MINOR CHANGES.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer  i
!
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
        ipob = idapo(inb)
!         minv = min(inva,invb,invc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   + cc(ipob+i) * bfac
        enddo
!         do 8 i=ipoc+minv+1,ipoc+invc
! 8       cc(i) = 0d0
        return
      endif
!      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      idacma = 0
      call daall1(idacma,'$$DACMA $$',nomax,nvmax)
      call dalin(ina,+1.d0,inb,bfac,idacma)
      call dacop(idacma,inc)
      call dadal1(idacma)
!
      return
      end subroutine
!
      subroutine dalin(ina,afac,inb,bfac,inc) bind(C, name="dalin_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, inb, inc
      real(c_double) afac, bfac

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoc, ipob, ipoa
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer  i
!
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
        ipob = idapo(inb)
!         minv = min(inva,invb,invc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
        enddo
!         do 8 i=ipoc+minv+1,ipoc+invc
! 8       cc(i) = 0d0
        return
      endif

      if(ina.eq.inc.or.inb.eq.inc) then
!        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call dalint(ina,afac,inb,bfac,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dalint(ina,afac,inb,bfac,inc)
      endif

      return
      end subroutine


      subroutine dalint(ina,afac,inb,bfac,inc)
      implicit none
      integer i,iamax,ibmax,icmax,illa,illb,illc,ilma,ilmb,             &
     &ilmc,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,                   &
     &ismax,ismin,ja,jb,minv,mchk
      integer(8) ia, ib, ic, is, ipoa, ipob, ipoc
      double precision afac,bfac,ccc,copf
!     ***************************************
!
!     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
!     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!
!      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
!
      if(nomax.eq.1) then
        ipoc = idapo(inc)
        ipoa = idapo(ina)
        ipob = idapo(inb)
!         minv = min(inva,invb,invc)
        do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
        enddo
!         do 8 i=ipoc+minv+1,ipoc+invc
! 8       cc(i) = 0d0
        return
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      ia = ipoa
      ib = ipob
      ic = ipoc - 1
      iamax = ipoa+illa-1
      ibmax = ipob+illb-1
      icmax = ipoc+ilmc-1
      ja = ia1(i1(ia)) + ia2(i2(ia))
      jb = ia1(i1(ib)) + ia2(i2(ib))
!
      if(ia.gt.iamax) then
        ismin = ib
        ismax = ibmax
        copf  = bfac
        goto 50
      endif
      if(ib.gt.ibmax) then
        ismin = ia
        ismax = iamax
        copf  = afac
        goto 50
      endif
!
!     COMPARING
!     *********
!
  10  continue
!      if(ja-jb) 30,20,40
      mchk=ja-jb
      if(mchk.lt.0) goto 30
      if(mchk.eq.0) goto 20
      if(mchk.gt.0) goto 40
!
!     ADDING TWO TERMS
!     ****************
!
  20  continue
      ccc = cc(ia)*afac + cc(ib)*bfac
      if(dabs(ccc).lt.eps) goto 25
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 25
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
  25  continue
      ia = ia + 1
      ib = ib + 1
      if(ia.gt.iamax) then
        ismin = ib
        ismax = ibmax
        copf  = bfac
        goto 50
      endif
      if(ib.gt.ibmax) then
        ismin = ia
        ismax = iamax
        copf  = afac
        goto 50
      endif
      ja = ia1(i1(ia)) + ia2(i2(ia))
      jb = ia1(i1(ib)) + ia2(i2(ib))
      goto 10
!
!     STORING TERM A
!     **************
!
  30  continue
      if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 35
      ccc = cc(ia)*afac
      if(dabs(ccc).lt.eps) goto 35
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
  35  continue
      ia = ia + 1
      if(ia.gt.iamax) then
        ismin = ib
        ismax = ibmax
        copf  = bfac
        goto 50
      endif
      ja = ia1(i1(ia)) + ia2(i2(ia))
      goto 10
!
!     STORING TERM B
!     **************
!
  40  continue
      if(ieo(ia1(i1(ib))+ia2(i2(ib))).gt.nocut) goto 45
      ccc = cc(ib)*bfac
      if(dabs(ccc).lt.eps) goto 45
      ic = ic + 1
      cc(ic) = ccc
      i1(ic) = i1(ib)
      i2(ic) = i2(ib)
  45  continue
      ib = ib + 1
      if(ib.gt.ibmax) then
        ismin = ia
        ismax = iamax
        copf  = afac
        goto 50
      endif
      jb = ia1(i1(ib)) + ia2(i2(ib))
      goto 10
!
!     COPYING THE REST
!     ****************
!
  50  continue
      do is=ismin,ismax
        if(ieo(ia1(i1(is))+ia2(i2(is))).gt.nocut) goto 60
        ccc = cc(is)*copf
        if(dabs(ccc).lt.eps) goto 60
        ic = ic + 1
        cc(ic) = ccc
        i1(ic) = i1(is)
        i2(ic) = i2(is)
  60    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
!
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DALIN, RESULT HAS TOO MANY TERMS '
        call dadeb(31,'ERR DALIN ',1)
      endif
!
      return
      end subroutine
!
      subroutine dafun(cf1,ina,inc) bind(C, name="dafun_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ina, inc
      character(c_char) cf1(4)

      integer illc,ilmc,incc,inoc,invc,j
      integer(8) ipoc
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      character cf*4

      do j=1, 4
         cf(j:j) = cf1(j)
      enddo
      if(ina.eq.inc) then
!       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',nomax,nvmax)
        call dafunt(cf,ina,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dafunt(cf,ina,inc)
      endif

      return
      end subroutine

      subroutine dafunt(cf,ina,inc)
      implicit none
      integer i,illa,illc,ilma,ilmc,ina,inc,ind,inoa,inoc,inon,inva,    &
     &invc,ipow,iscr,jj,lfun,no
      integer(8) ipoa, ipoc
      double precision a0,a1,a2,a3,a4,a5,ca,e1,                         &
     &e2,ea,era,p,ra,rpi4,sa,scr,                                       &
     &t,xf
!     ****************************
!
!     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
!     AND STORES THE RESULT IN C.
!     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
!     THIS HAS TO BE FIXED IN THE FUTURE.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      character cf*4,cfh*4,abcs*26,abcc*26
      dimension xf(0:lno),jj(lnv)
!
      data jj /lnv*0/
      data abcs /'abcdefghijklmnopqrstuvwxyz'/
      data abcc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
!
      if(cf(1:1).eq.' ') then
        cfh(1:3) = cf(2:4)
        cfh(1:4) = ' '
        cf = cfh
      endif
!
      do i=1,4
        ind = index(abcs,cf(i:i))
        if(ind.ne.0) cf(i:i) = abcc(ind:ind)
      enddo
!      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!     CASE OF NV = 0 WHICH MEANS COORDINATEWISE OPERATION
!     ***************************************************
!
!     CASE OF NV > 0 WHICH MEANS DIFFERENTIAL ALGEBRAIC OPERATION
!     ***********************************************************
!
      if(cf.eq.'SQR ') then
        call dasqr(ina,inc)
        return
      endif
!
!     ALLOCATE VARIABLES, PICK ZEROTH ORDER TERM
!     ******************************************
!
      ipow = 0
      inon = 0
      iscr = 0
!
      call daall1(ipow,'$$DAFUN1$$',nomax,nvmax)
      call daall1(inon,'$$DAFUN2$$',nomax,nvmax)
      call daall1(iscr,'$$DAFUN3$$',nomax,nvmax)
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      call dapek(ina,jj,a0)
!
!      no = min(nocut,inoa,inoc)
      no = min(nocut,nomax)
!
!     BRANCHING TO DIFFERENT FUNCTIONS
!     ********************************
!
      if(cf.eq.'INV ') then
!        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
        if(a0.eq.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        xf(0) = 1.d0/a0
        do i=1,no
          xf(i) = -xf(i-1)/a0
        enddo
!
      elseif(cf.eq.'SQRT') then
!        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)
        if(a0.le.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        ra = dsqrt(a0)
        xf(0) = ra
        do i=1,no
          xf(i) = -xf(i-1)/a0/dble(2*i)*dble(2*i-3)
        enddo
!
      elseif(cf.eq.'ISRT') then
!        1/SQRT(A0+P) = 1/SQRT(A0)*(1-1/2(P/A0)+3/8*(P/A0)**2-...)
        if(a0.le.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        era = 1.d0/dsqrt(a0)
        xf(0) = era
        do i=1,no
          xf(i) = -xf(i-1)/a0/dble(2*i)*dble(2*i-1)
        enddo
!
      elseif(cf.eq.'EXP ') then
!        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
        ea  = dexp(a0)
        xf(0) = ea
        do i=1,no
          xf(i) = xf(i-1)/dble(i)
        enddo
!
      elseif(cf.eq.'LOG ') then
!        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)
        if(a0.le.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        ea  = dlog(a0)
        xf(0) = ea
        xf(1) = 1.d0/a0
        do i=2,no
          xf(i) = -xf(i-1)/a0/dble(i)*dble(i-1)
        enddo
!
      elseif(cf.eq.'SIN ') then
!        SIN(A0+P) = SIN(A0)*(1-P**2/2!+P**4/4!) + COS(A0)*(P-P**3/3!+P**5/5!)
        sa  = dsin(a0)
        ca  = dcos(a0)
        xf(0) = sa
        xf(1) = ca
        do i=2,no
          xf(i) = -xf(i-2)/dble(i*(i-1))
        enddo
!
      elseif(cf.eq.'COS ') then
!        COS(A0+P) = COS(A0)*(1-P**2/2!+P**4/4!) - SIN(A0)*(P-P**3/3!+P**5/5!)
        sa  = dsin(a0)
        ca  = dcos(a0)
        xf(0) = ca
        xf(1) = -sa
        do i=2,no
          xf(i) = -xf(i-2)/dble(i*(i-1))
        enddo
!
      elseif(cf.eq.'SIRX') then
!        SIN(SQRT(P))/SQRT(P) = 1 - P/3! + P**2/5! - P**3/7! + ...
        if(a0.ne.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        xf(0)=1.d0
        do i=1,no
          xf(i) = -xf(i-1)/dble(2*i*(2*i+1))
        enddo
!
      elseif(cf.eq.'CORX') then
!        COS(SQRT(P)) = 1 - P/2! + P**2/4! - P**3/6! + ...
        if(a0.ne.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        xf(0)=1.d0
        do i=1,no
          xf(i) = -xf(i-1)/dble(2*i*(2*i-1))
        enddo
!
      elseif(cf.eq.'SIDX') then
!        SIN(P)/P = 1 - P**2/3! + P**4/5! - P**6/7! + ...
        if(a0.ne.0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        xf(0)=1.d0
        xf(1)=0.d0
        do i=2,no
          xf(i) = -xf(i-2)/dble(i*(i+1))
        enddo
!
      elseif(cf.eq.'TAN ') then
        if(dabs(dcos(a0)).lt.epsmac) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        sa  = dsin(a0)
        ca  = dcos(a0)
        xf(0) = sa/ca
        xf(1) = 1.d0/ca/ca
        xf(2) = 2.d0*sa/ca/ca/ca/2.d0
        xf(3) = (2.d0*ca*ca+6.d0*sa*sa)/ca/ca/ca/ca/6.d0
        xf(4) = (16*sa+8.d0*sa*sa*sa)/ca/ca/ca/ca/ca/24.d0
        xf(5) = (16.d0*ca*ca+24.d0*ca*ca*sa*sa+80.d0*sa*sa+             &
     &  40.d0*sa*sa*sa*sa)/ca/ca/ca/ca/ca/ca/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
          stop
        endif
      elseif(cf.eq.'COT ') then
        if(dabs(dsin(a0)).lt.epsmac) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        sa  = dsin(a0)
        ca  = dcos(a0)
        xf(0) = ca/sa
        xf(1) = -1.d0/sa/sa
        xf(2) = 2.d0*ca/sa/sa/sa/2.d0
        xf(3) = -(2.d0*sa*sa+6.d0*ca*ca)/sa/sa/sa/sa/6.d0
        xf(4) = (16*ca+8.d0*ca*ca*ca)/sa/sa/sa/sa/sa/24.d0
        xf(5) = -(16.d0*sa*sa+24.d0*sa*sa*ca*ca+80.d0*ca*ca+            &
     &  40.d0*ca*ca*ca*ca)/sa/sa/sa/sa/sa/sa/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
          stop
        endif
      elseif(cf.eq.'ASIN') then
        if((1.d0-dabs(a0)).lt.0.d0) then
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          lfun = 0
          return
        endif
        xf(0) = dasin(a0)
        xf(1) = (1.d0-a0*a0)**(-0.5d0)
        xf(2) = a0*xf(1)**3.d0/2.d0
        xf(3) = (1+2.d0*a0*a0)*xf(1)**5.d0/6.d0
        xf(4) = (9.d0*a0+6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
        xf(5) = (9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          stop
        endif
      elseif(cf.eq.'ACOS')then
        if((1.d0-dabs(a0)).lt.0.d0) then
          call dadeb(31,'ERR DAFUN ',1)
          write(*,1000) cf,ina,a0
          lfun = 0
          return
        endif
        xf(0) =  dacos(a0)
        scr =  (1.d0-a0*a0)**(-0.5d0)
        xf(1) =  -scr
        xf(2) = -a0*scr**3.d0/2.d0
        xf(3) = -(1+2.d0*a0*a0)*scr**5.d0/6.d0
        xf(4) = -(9.d0*a0+6.d0*a0*a0*a0)*scr**7.d0/24.d0
        xf(5) = -(9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*scr**9.d0/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ATAN') then
!        ATAN(A0+P) = ATAN(A0)+1/(1+A0**2)*P-A0/(1+A0**2)**2*P**2+....)
        xf(0) = datan(a0)
        xf(1) = 1.d0/(1.d0+a0*a0)
        xf(2) = -a0*(xf(1)*xf(1))
        xf(3) = (a0*a0-1.d0/3.d0)*xf(1)**3
        xf(4) = (a0-a0*a0*a0)*xf(1)**4
        xf(5) = (1.d0/5.d0+a0**4-2.d0*a0*a0)*xf(1)**5
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ACOT') then
        xf(0) = 2.d0*datan(1.d0)-datan(a0)
        scr = 1.d0/(1.d0+a0*a0)
        xf(1) = -scr
        xf(2) = a0*(scr*scr)
        xf(3) = -(a0*a0-1.d0/3.d0)*scr**3
        xf(4) = -(a0-a0*a0*a0)*scr**4
        xf(5) = -(1.d0/5.d0+a0**4-2.d0*a0*a0)*scr**5
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'SINH') then
        sa  = dsinh(a0)
        ca  = dcosh(a0)
        xf(0) = sa
        xf(1) = ca
        xf(2) = sa/2.d0
        xf(3) = ca/6.d0
        xf(4) = sa/24.d0
        xf(5) = ca/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'COSH') then
        sa  = dsinh(a0)
        ca  = dcosh(a0)
        xf(0) = ca
        xf(1) = sa
        xf(2) = ca/2.d0
        xf(3) = sa/6.d0
        xf(4) = ca/24.d0
        xf(5) = sa/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'TANH') then
        sa  = dsinh(a0)
        ca  = dcosh(a0)
        xf(0) = sa/ca
        xf(1) = 1.d0/ca/ca
        xf(2) = -2.d0*sa/ca/ca/ca/2.d0
        xf(3) = (-2.d0*ca*ca+6.d0*sa*sa)/ca/ca/ca/ca/6.d0
        xf(4) = (16*sa-8.d0*sa*sa*sa)/ca/ca/ca/ca/ca/24.d0
        xf(5) = (16.d0*ca*ca-24.d0*ca*ca*sa*sa-80.d0*sa*sa+             &
     &  40.d0*sa*sa*sa*sa)/ca/ca/ca/ca/ca/ca/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'COTH') then
        if(dabs(dsinh(a0)).lt.epsmac) then
          lfun = 0
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          return
        endif
        sa  = dsinh(a0)
        ca  = dcosh(a0)
        xf(0) = ca/sa
        xf(1) = -1.d0/sa/sa
        xf(2) =  2.d0*ca/sa/sa/sa/2.d0
        xf(3) = (2.d0*sa*sa-6.d0*ca*ca)/sa/sa/sa/sa/6.d0
        xf(4) = (16*ca+8.d0*ca*ca*ca)/sa/sa/sa/sa/sa/24.d0
        xf(5) = (16.d0*sa*sa+24.d0*sa*sa*ca*ca-80.d0*ca*ca-             &
     &  40.d0*ca*ca*ca*ca)/sa/sa/sa/sa/sa/sa/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ASNH') then
        xf(0) = dlog(a0+dsqrt(a0*a0+1.d0))
        xf(1) = (1.d0+a0*a0)**(-0.5d0)
        xf(2) = -a0*xf(1)**3.d0/2.d0
        xf(3) = (2.d0*a0*a0-1.d0)*xf(1)**5.d0/6.d0
        xf(4) = (9.d0*a0-6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
        xf(5) = (9.d0-72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ACSH') then
        if((1.d0-a0).ge.0.d0) then
          lfun = 0
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          return
        endif
        xf(0) = dlog(a0+dsqrt(a0*a0-1.d0))
        xf(1) = (a0*a0-1.d0)**(-0.5d0)
        xf(2) = -a0*xf(1)**3.d0/2.d0
        xf(3) = (2.d0*a0*a0+1.d0)*xf(1)**5.d0/6.d0
        xf(4) = (-9.d0*a0-6.d0*a0*a0*a0)*xf(1)**7.d0/24.d0
        xf(5) = (9.d0+72.d0*a0*a0+24.d0*a0*a0*a0*a0)*xf(1)**9.d0/120.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ATNH') then
        if((dabs(a0)-1.d0).ge.0.d0) then
          lfun = 0
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          return
        endif
        xf(0) =  0.5d0*dlog((1+a0)/(1-a0))
        xf(1) =  1.d0/(1.d0-a0*a0)
        xf(2) =  a0*(xf(1)*xf(1))
        xf(3) = (a0*a0+1.d0/3.d0)*xf(1)**3
        xf(4) = (a0+a0*a0*a0)*xf(1)**4
        xf(5) = (1.d0/5.d0+a0**4+2.d0*a0*a0)*xf(1)**5
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      elseif(cf.eq.'ACTH') then
        if(1.d0-dabs(a0).ge.0.d0) then
          lfun = 0
          write(*,1000) cf,ina,a0
          call dadeb(31,'ERR DAFUN ',1)
          return
        endif
        xf(0) =  0.5d0*dlog((a0+1)/(a0-1))
        scr =  1.d0/(-1.d0+a0*a0)
        xf(1) = -scr
        xf(2) =  a0*(scr*scr)
        xf(3) = (-a0*a0-1.d0/3.d0)*scr**3.d0
        xf(4) = (a0+a0*a0*a0)*scr**4.d0
        xf(5) = (-1.d0/5.d0-a0**4-2.d0*a0*a0)*scr**5.d0
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
!      ELSEIF(CF.EQ.'ABF ') THEN
!
!     DIESE FUNKTION BESCHREIBT DEN FELDABFALL BEI IONENOPTISCHEN ELEMENTEN
!     ABF=1/(1+EXP(A0+X))
!        =1/(1+EXP(A0)*(1-EXP(A0)/(1+EXP(A0))*X+....)
!         XF(0) = 1.D0/(1+DEXP(A0))
!         E1  = DEXP(A0)*X1
!         E1  = DEXP(A0)*Xf(0)
!         E2  = E1 * E1
!         E3  = E2 * E1
!         E4  = E3 * E1
!         E5  = E4 * E1
!         XF(1) = X1*(-E1)
!         XF(2) = X1*(-0.5D0* E1 + E2)
!         XF(3) = X1*(-E1/6.D0 + E2 - E3)
!         XF(4) = X1*(-E1/24.D0 + E2*7.D0/12.D0 - E3*3.D0/2.D0 + E4)
!         XF(5) = X1*(-E1/120.D0 + E2/4.D0 - E3*5.D0/4.D0 +
!     *         E4*2.D0 - E5)
!         IF(NO.GT.5) THEN
!            write(6,*)'ERROR IN DAFUN, ',CF, ' ONLY UP TO NO = 5'
!            CALL DADEB(31,'ERR DAFUN ',1)
!         ENDIF
!      ELSEIF(CF.EQ.'GAUS') THEN
!
!     DIESE FUNKTION BESCHREIBT DIE ENTWICKLUNG VON EXP(-X*X)
!
!         XF(0) = DEXP(-A0*A0)
!         XF(1) = -2.D0*A0*X1
!         XF(2) = (-1.D0+2.D0*A0*A0)*X1
!         XF(3) = (12.D0*A0-8.D0*A0*A0*A0)/6.D0*X1
!         XF(4) = (16.D0*A0*A0*A0*A0-48.D0*A0*A0+12.D0)/24.D0*X1
!         XF(5) = (-32.D0*A0*A0*A0*A0*A0+160.D0*A0*A0*A0-120.D0*A0)/
!     *           120.D0*X1
!         IF(NO.GT.5) THEN
!            write(6,*)'ERROR IN DAFUN, ',CF, ' ONLY UP TO NO = 5'
!            CALL DADEB(31,'ERR DAFUN ',1)
!         ENDIF
      elseif(cf.eq.'ERF ') then
!
!    ERF(X) STELLT DAS INTEGRAL VON 0 BIS X VON [ 2/SQRT(PI) * EXP(-X*X) ]
!    DAR
!
        e1 = dexp(-a0*a0)
        a1 = .254829592d0
        a2 = -.284496736d0
        a3 = 1.421413741d0
        a4 = -1.453152027d0
        a5 = 1.061405429d0
        p  = .3275911d0
        rpi4 = sqrt(datan(1.d0))
        t  = 1.d0/(1.d0+p*a0)
        e2 = 1.d0-t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))*e1
        xf(0)= e2
        xf(1) = e1/rpi4
        xf(2) = -a0*e1/rpi4
        xf(3) = (-2.d0+4.d0*a0*a0)/6.d0*e1/rpi4
        xf(4) = (12.d0*a0-8.d0*a0*a0*a0)/24.d0*e1/rpi4
        xf(5) = (16.d0*a0*a0*a0*a0-48.d0*a0*a0+12.d0)/120.d0*e1/rpi4
        if(no.gt.5) then
          write(6,*)'ERROR IN DAFUN, ',cf, ' ONLY UP TO NO = 5'
          call dadeb(31,'ERR DAFUN ',1)
        endif
      else
        write(6,*)'ERROR, UNSOPPORTED FUNCTION ',cf
      endif
!
      call dacon(inc,xf(0))
      call dacop(ina,inon)
      call dapok(inon,jj,0.d0)
      call dacon(ipow,1.d0)
!
      do i=1,min(no,nocut)
!
        call damul(inon,ipow,iscr)
        call dacop(iscr,ipow)
        call dacma(inc,ipow,xf(i),inc)
!
      enddo
!
 1000 format('ERROR IN DAFUN, ',a4,' DOES NOT EXIST FOR VECTOR ',i10,   &
     &'CONST TERM  = ',e12.5)
!
      call dadal1(iscr)
      call dadal1(inon)
      call dadal1(ipow)
!
      return
      end subroutine
!

      subroutine daabs(ina,anorm) bind(C, name="daabs_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina
      real(c_double) anorm

      integer illa,ilma,inoa,inva
      integer(8) i, ipoa
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
      anorm = 0.d0
      do i=ipoa,ipoa+illa-1
        anorm = anorm + dabs(cc(i))
      enddo
!
      return
      end subroutine
!
      subroutine daabs2(ina,anorm) bind(C, name="daabs2_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina
      real(c_double) anorm

      integer illa,ilma,inoa,inva
      integer(8) i, ipoa
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
      anorm = 0.d0
      do i=ipoa,ipoa+illa-1
        anorm = anorm + dabs(cc(i)**2)
      enddo
!
      return
      end subroutine
!

      subroutine dacom(ina,inb,dnorm)
      implicit none
      integer idacom,illc,ilmc,ina,inb,inc,inoc,invc,ipoc
      double precision dnorm
!     *******************************
!
!     THIS SUBROUTINE COMPARES TWO DA VECTORS BY RETURNING THE NORM
!     OF THE DIFFERENCE
!
      idacom = 0
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      call daall1(idacom,'$$DACOM $$',inoc,invc)
      call dasub(ina,inb,idacom)
      call daabs(idacom,dnorm)
      call dadal1(idacom)
!
      return
      end subroutine
!

      subroutine dapos(ina,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8)  ia, ib, ipoa, ipob
!     *************************
!
!     THIS SUBROUTINE MAKES THE SIGNS OF ALL THE COEFFICIENTS OF A POSITIVE
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
      ib = ipob - 1
!
      do ia = ipoa,ipoa+illa-1
!
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
        ib     = ib + 1
        cc(ib) = dabs(cc(ia))
        i1(ib) = i1(ia)
        i2(ib) = i2(ia)
!
 100    continue
      enddo
!
      idall(inb) = ib - ipob + 1
      if(idall(inb).gt.idalm(inb)) then
        write(6,*)'ERROR IN DAPOS '
        call dadeb(31,'ERR DAPOS ',1)
      endif
!
      return
      end subroutine
!
      subroutine dacct(ma,ia,mb,ib,mc,ic) bind(C, name="dacct_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ma(*), mb(*), mc(*), ia, ib, ic

      integer i,ij,illc,ilmc,inoc,invc
      integer(8) ipoc
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer mon(lnv)

      if(ma(1).eq.mc(1).or.mb(1).eq.mc(1)) then
        call dainf(mc(1),inoc,invc,ipoc,ilmc,illc)
        do ij=1,ic
          mon(ij)=0
        enddo
        call daall(mon,ic,'$$DAJUNK$$',inoc,invc)
        call dacctt(ma,ia,mb,ib,mon,ic)
        do i=1,ic
          call dacop(mon(i),mc(i))
        enddo
        call dadal(mon,ic)
      else
        call dacctt(ma,ia,mb,ib,mc,ic)
      endif

      return
      end subroutine

      subroutine dacctt(mb,ib,mc,ic,ma,ia)
      implicit none
      integer i,ia,ib,ic,iia,iib,iic,illa,illb,illc,ilma,ilmb,ilmc,inoa,&
     &inob,inoc,inva,invb,invc,iv,jl,jv
      integer(8) ipoa,ipob,ipoc
      double precision ccf
!     ***********************************
!
!     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
!     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
!     DA VECTORS, RESPECTIVELY.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!      INTEGER MON(LNO+1),ICC(LNV),MB(*),MC(*),MA(*)
!ETIENNE
!      integer mon(lno+1),icc(lno),mb(*),mc(*),ma(*)
      integer mon(lno+1),icc(lnv),mb(*),mc(*),ma(*)
!ETIENNE
!
!     CONSISTENCY CHECKS
!     ******************
!
      iia = ma(1)
      iib = mb(1)
      iic = mc(1)
      call dainf(iia,inoa,inva,ipoa,ilma,illa)
      call dainf(iib,inob,invb,ipob,ilmb,illb)
      call dainf(iic,inoc,invc,ipoc,ilmc,illc)
!
      call damch(ma,ia)
      call damch(mb,ib)
!
      if(ia.ne.ib) then
        write(6,*)'ERROR IN DACCT, IA .NE. IB'
        call dadeb(31,'ERR DACCT1',1)
      elseif(ic.ne.invb) then
        write(6,*)'ERROR IN DACCT, IC.NE.INVB'
        call dadeb(31,'ERR DACCT2',1)
      endif
!
!     ALLOCATING LOCAL VECTORS AND CALLING MTREE
!     ******************************************
!
      do i=1,ib
        icc(i) = 0
      enddo
!
      do i=1,nomax+1
        mon(i) = 0
      enddo
!
      call daall(icc,ib,'$$DACCT $$',nomax,nvmax)
      call daall(mon,nomax+1,'$$DAMON $$',inoc,invc)
!
      call mtree(mb,ib,icc,ib)
!
!     PERFORMING CONCATENATION
!     ************************
!
      do i=1,ia
        call dacon(ma(i),cc(idapo(icc(i))))
      enddo
!
      call dacon(mon(1),1.d0)
!
      do i=1,idall(icc(1))-1
!
        jl = i1(idapo(icc(1))+i)
        jv = i2(idapo(icc(1))+i)
!
        call damul(mon(jl),mc(jv),mon(jl+1))
!
        do iv=1,ia
!
          ccf = cc(idapo(icc(iv))+i)
          if(dabs(ccf).gt.eps) call dacma(ma(iv),mon(jl+1),ccf,ma(iv))
!
        enddo
      enddo
!
      call dadal(mon,nomax+1)
      call dadal(icc,ib)
!
      return
      end subroutine

      subroutine mtree(mb,ib,mc,ic) bind(C, name="mtree_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) mb(*), ib, mc(*), ic

      integer ic1,ic2,ichk,ii,iib,iic,illb,illc,                        &
     &ilmb,ilmc,inob,inoc,invb,invc,j,jl,jnon,nterm,ntermf
      integer(8) i, ibi, ib1, icc, ipob, ipoc
      double precision apek,bbijj,chkjj
!     *****************************
!
!     THIS SUBROUTINE IS USED FOR CONCATENATION AND TRACKING OF VECTORS
!     THROUGH A DA MAP. IT COMPUTES THE TREE THAT HAS TO BE TRANSVERSED
!     MB IS THE DA MATRIX WITH IA TERMS. THE OUTPUT MC IS A CA MATRIX WHICH
!     CONTAINS COEFFICIENTS AND CONTROL INTEGERS USED FOR THE TRAVERSAL.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv),jv(0:lno)
!
!     CONSISTENCY CHECKS
!     ******************
!
      iib = mb(1)
      iic = mc(1)
      call dainf(iib,inob,invb,ipob,ilmb,illb)
      call dainf(iic,inoc,invc,ipoc,ilmc,illc)
!
      call damch(mb,ib)
      call damch(mc,ic)
!
      if(ib.ne.ic) then
        write(6,*)'ERROR IN MTREE, IB .NE. IC'
        call dadeb(31,'ERR MTREE1',1)
      endif
!
!     ALLOCATING LOCAL VECTORS
!     ************************
!
      ichk = 0
      call daall1(ichk,'$$MTREE $$',nomax,nvmax)
!
!     FIND ALL THE ENTRIES TO BE LOOKED FOR
!     *************************************
!
      call daclr(1)
!
      cc(1) = 1.d0
!
      do i=1,ib
        if(nomax.eq.1) then
          do ib1 = 2,7
            cc(ib1) = 1d0
          enddo
        else
          do ibi = idapo(mb(i)),idapo(mb(i))+idall(mb(i))-1
            icc = ia1(i1(ibi)) + ia2(i2(ibi))
            if(ieo(icc).gt.inob) goto 90
            cc(icc) = 1.d0
   90       continue
          enddo
        endif
      enddo
!
      do ii=1,inob
!
!     SEARCHING FOR FATHER FOR EACH TERM
!
        do i=1,nmmax
          if(cc(i).lt.0.5d0) goto 140
!
          jnon = 0
          call dancd(i1(i),i2(i),jj)
          do j=1,invb
            if(jj(j).eq.0) goto 130
            jnon = j
            jj(j) = jj(j) - 1
            call dadcd(jj,ic1,ic2)
            apek = cc(ia1(ic1)+ia2(ic2))
            jj(j) = jj(j) + 1
            if(apek.ge.0.5d0) goto 140
  130       continue
          enddo
!
          if(jnon.eq.0) goto 140
!
!     TERM IS AN ORPHAN, SO CREATE FOSTER FATHER
!
          jj(jnon) = jj(jnon) - 1
          call dadcd(jj,ic1,ic2)
          cc(ia1(ic1)+ia2(ic2)) = 1.d0
!
  140     continue
        enddo
      enddo
!
      call dapac(ichk)
!ETIENNE      CALL DAPRI(ICHK,32)
!
!     SETTING UP TREE STRUCTURE
!     *************************
!
      ntermf = idall(ichk)
!
!     ZEROTH ORDER TERMS
!     ******************
!
      do i=1,lnv
        jj(i) = 0
      enddo
!
      do i=1,ib
        call dapek(mb(i),jj,bbijj)
        i1(idapo(mc(i))) = 0
        i2(idapo(mc(i))) = 0
        cc(idapo(mc(i))) = bbijj
      enddo
!
      call dapek(ichk,jj,chkjj)
      if(chkjj.gt.0.5d0) then
        call dapok(ichk,jj,-1.d0)
      else
        write(6,*)'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS ZERO'
        call dadeb(31,'ERR MTREE2',1)
      endif
!
      nterm = 1
!
!     HIGHER ORDER TERMS
!     ******************
!
      do jl=1,inob
        jv(jl) = 0
      enddo
!
      jl = 0
      chkjj = 1.d0
!
 200  continue
      if(jl.eq.0.and.chkjj.le.0.5d0) goto 250
      if(jl.lt.inob.and.chkjj.gt.0.5d0) then
        jl = jl + 1
        jj(1) = jj(1) + 1
        jv(jl) = 1
      elseif(jv(jl).eq.invb) then
        jj(jv(jl)) = jj(jv(jl)) - 1
        jv(jl) = 0
        jl = jl - 1
        chkjj = 0.d0
        goto 200
      else
        jj(jv(jl)) = jj(jv(jl)) - 1
        jv(jl) = jv(jl) + 1
        jj(jv(jl)) = jj(jv(jl)) + 1
      endif
!
      call dapek(ichk,jj,chkjj)
!
      if(chkjj.le.0.5d0) goto 200
!
      nterm = nterm + 1
      if(nterm.gt.idalm(mc(1))) then
        write(6,*)'ERROR IN MTREE, NTERM TOO LARGE'
        call dadeb(31,'ERR MTREE3',1)
      endif
!
      call dapok(ichk,jj,-1.d0)
!
!     write(6,*)'JL,JV = ',JL,JV(JL)
      do i=1,ib
        call dapek(mb(i),jj,bbijj)
        i1(idapo(mc(i))+nterm-1) = jl
        i2(idapo(mc(i))+nterm-1) = jv(jl)
        cc(idapo(mc(i))+nterm-1) = bbijj
      enddo
!
      goto 200
!
 250  continue
!
      do i=1,ib
        idall(mc(i)) = nterm
      enddo
!
!     PERFORMING CROSS CHECKS
!     ***********************
!
      if(nterm.ne.ntermf.or.nterm.ne.idall(ichk)) then
        write(6,*)'ERROR IN MTREE, NTERM, NTERMF, IDALL(ICHK) =  '      &
     &  ,nterm,ntermf,idall(ichk)
        call dadeb(31,'ERR MTREE4',1)
      endif
!
      do i=idapo(ichk),idapo(ichk)+nterm-1
        if(dabs(cc(i)+1.d0).gt.epsmac) then
          write(6,*)'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
          call dadeb(31,'ERR MTREE5',1)
        endif
      enddo
!
      call dadal1(ichk)
!
      return
      end subroutine

      subroutine ppushprint(mc,ic,mf,jc,line)
      implicit none
      integer i,ic,iv,jc,jl,jv,mc,mf
      include "TPSALib_prm.f"
      dimension mc(*)
      character*20 line
      if(mf.le.0) return
      write(mf,*) 0,0,jc+1,0,line
      write(mf,*) 0,0,jc+1,0
      do i=1,ic
        jc=1+jc
        write(mf,*) jc,jl,jv,cc(idapo(mc(i)))
      enddo
!     xf(i) = cc(idapo(mc(i)))
!      xm(1) = 1.d0
      do i=1,idall(mc(1))-1
        jl = i1(idapo(mc(1))+i)
        jv = i2(idapo(mc(1))+i)
!      xx = xm(jl)*xi(jv)
!      xm(jl+1) = xx
        do iv=1,ic
          jc=1+jc
          write(mf,*) jc,jl,jv,cc(idapo(mc(iv))+i)
!        xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
        enddo
      enddo
      return
      end subroutine
!

      subroutine ppush(mc,ic,xi,xf) bind(C, name="ppush_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) mc, ic
      real(c_double) xi, xf

      integer i,iv,jl,jv
      double precision xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension mc(*),xf(*),xi(*),xm(lno+1) ,xt(lno)
!
      do i=1,ic
        xt(i)=xi(i)
      enddo
      do i=1,ic
        xf(i) = cc(idapo(mc(i)))
      enddo
!
      xm(1) = 1.d0
!
      do i=1,idall(mc(1))-1
!
        jl = i1(idapo(mc(1))+i)
        jv = i2(idapo(mc(1))+i)
        xx = xm(jl)*xt(jv)
        xm(jl+1) = xx
!
        do iv=1,ic
          xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
        enddo
      enddo
!
      return
      end subroutine

      subroutine ppush1(mc,xi,xf)
      implicit none
      integer i,iv,jl,jv,mc
      double precision xf,xi,xm,xt,xx
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
!     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension xi(*),xm(lno+1) ,xt(lno)
!
      do i=1,nvmax
        xt(i)=xi(i)
      enddo

      xf = cc(idapo(mc))
!
      xm(1) = 1.d0
!
      do i=1,idall(mc)-1
!
        jl = i1(idapo(mc)+i)
        jv = i2(idapo(mc)+i)
        xx = xm(jl)*xt(jv)
        xm(jl+1) = xx
!
        xf = xf + cc(idapo(mc)+i) * xx
      enddo
!
      return
      end subroutine

      subroutine dainv(ma,ia,mb,ib) bind(C, name="dainv_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ma(*), ia, mb(*), ib

      integer i,ij,illb,ilmb,inob,invb
      integer(8) ipob
      double precision x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv),ml(lnv)
!
      dimension x(lnv)
!

      do i=1,lnv
        jj(i)=0
      enddo

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)
        do i=1,ia
          call dapok(ma(i),jj,0.d0)
        enddo
        do ij=1,ib
          ml(ij)=0
        enddo
        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dainvt(ma,ia,ml,ib)
        do i=1,ib
          call dacop(ml(i),mb(i))
        enddo
        call dadal(ml,ib)
      else
        do i=1,ia
          call dapek(ma(i),jj,x(i))
          call dapok(ma(i),jj,0.d0)
        enddo
        call dainvt(ma,ia,mb,ib)
        do i=1,ia
          call dapok(ma(i),jj,x(i))
        enddo
      endif

      return
      end subroutine

      subroutine dainvt(ma,ia,mb,ib)
      implicit none
      integer i,ia,ib,ie,ier,illa,illb,ilma,ilmb,inoa,inob,inva,invb,   &
     &j,k,nocut0
      integer(8) ipoa, ipob
      double precision aa,ai,amjj,amsjj,prod
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv),ms(lnv),ml(lnv),ma(*),mb(*)
!
      dimension aa(lnv,lnv),ai(lnv,lnv)
!
      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
      call dainf(mb(1),inob,invb,ipob,ilmb,illb)
!
!     CONSISTENCY CHECKS
!     ******************
!
      call damch(ma,ia)
      call damch(mb,ib)
!etienne
      do ie=1,ib
        call dacon(mb(ie),0.d0)
      enddo
!etienne
!
      if(ia.ne.ib) then
        write(6,*)'ERROR IN DAINV, IA .NE. IB'
        call dadeb(31,'ERR DAINV1',1)
      elseif(ia.ne.inva.or.ib.ne.invb) then
        write(6,*)'ERROR IN DAINV, IA.NE.INVA.OR.IB.NE.INVB'
        call dadeb(31,'ERR DAINV2',1)
      endif
!
!     ALLOCATING LOCAL VECTORS
!     ************************
!
      do i=1,ia
        ms(i) = 0
        ml(i) = 0
      enddo
!
      call daall(ms,ia,'$$INV   $$',inoa,inva)
      call daall(ml,ia,'$$INVL  $$',inoa,inva)
!
!     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
!     ********************************************************
!
      do i=1,ib
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          enddo
          jj(j) = 1
          call dapek(ma(i),jj,amjj)
          if(dabs(amjj).gt.eps) call dapok(ma(i),jj,0.d0)
          aa(i,j) = amjj
        enddo
        call dacmu(ma(i),-1.d0,ma(i))
      enddo
!
!     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
!     **********************************************************
!
      call matinv(aa,ai,ia,lnv,ier)
!
      if(ier.eq.132) then
        write(6,*)'ERROR IN ROUTINE DAINV'
        call dadeb(31,'ERR DAINV3',1)
      endif
!
      ier = 0
      do i=1,ib
        do j=1,ib
          prod = 0.d0
          do k=1,ib
            jj(k) = 0
            prod = prod + aa(i,k)*ai(k,j)
          enddo
          if(i.eq.j) prod = prod - 1.d0
          if(dabs(prod).gt.100*epsmac) then
            write(6,*)'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD',&
     &      ' = ',                                                      &
     &      i,j,prod,epsmac,eps
            ier = 1
!ETIENNE
            return
!ETIENNE
          endif
          jj(j) = 1
          call dapok(mb(i),jj,ai(i,j))
          call dapok(ml(i),jj,ai(i,j))
        enddo
      enddo
!
      if(ier.eq.1) call dadeb(31,'ERR DAINV4',1)
!
!     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
!     ****************************************************
!
!     MB (OF ORDER I) = A1^-1 o [ E - ANL (NONLINEAR) o MB (OF ORDER I) ]
!
      nocut0 = nocut
!
      do i=2,nocut
!
        nocut = i
!
        call dacct(ma,ia,mb,ib,ms,ia)
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          enddo
          jj(j) = 1
          call dapek(ms(j),jj,amsjj)
          call dapok(ms(j),jj,amsjj+1.d0)
        enddo
!
        call dacct(ml,ia,ms,ia,mb,ib)
!
      enddo
!
      nocut = nocut0
!
!     FLIPPING BACK SIGN OF A, FILLING UP FIRST ORDER PART AGAIN
!     **********************************************************
!
      do i=1,ib
        call dacmu(ma(i),-1.d0,ma(i))
        do j=1,ib
          do k=1,ib
            jj(k) = 0
          enddo
          jj(j) = 1
          call dapok(ma(i),jj,aa(i,j))
        enddo
      enddo
!
      call dadal(ml,ia)
      call dadal(ms,ia)
!
      return
      end subroutine
!
      subroutine matinv(a,ai,n,nmx,ier) bind(C, name="matinv_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) nmax
      parameter (nmax=400)
      integer(c_long) n, nmx, ier
      real(c_double) a(nmx,nmx),ai(nmx,nmx)

      integer i,indx,j
      double precision aw,d
!     *********************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX A AND STORES THE RESULT IN AI
!     INPUT  A   - SAVED
!            N   - ORDER OF MATRIX < 100
!     OUTPUT AI  - A INVERSE
!            IER - 0 NO ERROR
!                  132 ZERO DETERMINANT
!
      dimension aw(nmax,nmax),indx(nmax)

      do i=1,n
        do j=1,n
          aw(i,j) = a(i,j)
          ai(i,j) = 0.0
        enddo
        ai(i,i) = 1.d0
      enddo

      call ludcmp(aw,n,nmax,indx,d,ier)
      if (ier .eq. 132) return
      do j=1,n
        call lubksb(aw,n,nmax,indx,ai(1,j),nmx)
      enddo
!
      return
      end subroutine
!
      subroutine ludcmp(a,n,np,indx,d,ier)
      implicit none
      integer i,ier,imax,indx,j,k,n,nmax,np
      double precision a,aamax,d,dum,sum,tiny,vv
!     ************************************
!
!     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
!     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
!           NP: PHYSICAL DIMENSION OF A
!           INDX: ROW PERMUTATION VECTOR
!           D: EVEN OR ODD ROW INTERCHANGES
!
!     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
!
      parameter (nmax = 400, tiny = 1.0e-20)
      dimension a(np,np), indx(np), vv(nmax)
      ier=0
      d=1.d0
      do i=1,n
        aamax=0.d0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        enddo
        if(aamax.eq.0.d0) then
          ier=132
          return
        endif
        vv(i)=1.d0/aamax
      enddo
      do j=1,n
        if(j.gt.1) then
          do i=1,j-1
            sum=a(i,j)
            if(i.gt.1) then
              do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
              enddo
              a(i,j)=sum
            endif
          enddo
        endif
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          if (j.gt.1) then
            do k=1,j-1
              sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n) then
          if(a(j,j).eq.0.d0) a(j,j)=tiny
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      if(a(n,n).eq.0.d0) a(n,n)=tiny
      return
      end subroutine
!
      subroutine lubksb(a,n,np,indx,b,nmx)
      implicit none
      integer i,ii,indx,j,ll,n,nmx,np
      double precision a,b,sum
!     ************************************
!
!     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
!     INPUT A: NXN MATRIX IN lu FORM GIVEN BY LUDCMP
!           NP: PHYSICAL DIMENSION OF A
!           INDX: ROW PERMUTATION VECTOR
!           D: EVEN OR ODD ROW INTERCHANGES
!           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
!
!     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
!
      dimension a(np,np), indx(np), b(nmx)
      ii = 0
      do i=1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if(ii.ne.0) then
          do j=ii,i-1
            sum = sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.d0) then
          ii = i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        if(i.lt.n) then
          do j=i+1,n
            sum = sum-a(i,j)*b(j)
          enddo
        endif

        b(i)=sum/a(i,i)

      enddo
      return
      end subroutine
!

      subroutine dapin(ma,ia,mb,ib,jx) bind(C, name="dapin_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ma(*), ia, mb(*), ib, jx(*)

      integer i,ij,illb,ilmb,inob,invb
      integer(8) ipob
      double precision x
!     *****************************
!
!     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
!     STORES THE RESULT IN MI
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv),ml(lnv)
!
      dimension x(lnv)
!

      do i=1,lnv
        jj(i)=0
      enddo

      if(ma(1).eq.mb(1)) then
        call dainf(mb(1),inob,invb,ipob,ilmb,illb)
        do i=1,ia
          call dapok(ma(i),jj,0.d0)
        enddo
        do ij=1,ib
          ml(ij)=0
        enddo
        call daall(ml,ib,'$$DAJUNK$$',inob,invb)
        call dapint(ma,ia,ml,ib,jx)
        do i=1,ib
          call dacop(ml(i),mb(i))
        enddo
        call dadal(ml,ib)
      else
        do i=1,ia
          call dapek(ma(i),jj,x(i))
          call dapok(ma(i),jj,0.d0)
        enddo
        call dapint(ma,ia,mb,ib,jx)
        do i=1,ia
          call dapok(ma(i),jj,x(i))
        enddo
      endif

      return
      end subroutine

      subroutine dapint(ma,ia,mb,ib,jind)
      implicit none
      integer i,ia,ib,illa,ilma,inoa,inva,k
      integer(8) ipoa
!     **********************************
!
!     THIS SUBROUTINE PERFORMS A PARTIAL INVERSION OF THE ROWS MARKED WITH
!     NONZERO ENTRIES IN JJ OF THE MATRIX A. THE RESULT IS STORED IN B.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv),jind(*),ma(*),mb(*),mn(lnv),mi(lnv),me(lnv)
!
      call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
!

      do i=1,ia
        mn(i) = 0
        mi(i) = 0
        me(i) = 0
      enddo
!
      call daall(mn,ia,'$$PIN1  $$',inoa,inva)
      call daall(mi,ia,'$$PIN2  $$',inoa,inva)
      call daall(me,ia,'$$PIN3  $$',inoa,inva)
!
      do i=1,ia
        do k=1,nvmax
          jj(k) = 0
        enddo
        jj(i) = 1
        call dapok(me(i),jj,1.d0)
      enddo
!
      do i=1,ia
        call dacop(ma(i),mn(i))
        if(jind(i).eq.0) call dacop(me(i),mn(i))
      enddo
!
      call dainv(mn,ia,mi,ia)
!
      do i=1,ia
        if(jind(i).eq.0) call dacop(ma(i),me(i))
      enddo
!
      call dacct(me,ia,mi,ia,mb,ib)
!
      call dadal(me,ia)
      call dadal(mi,ia)
      call dadal(mn,ia)
!
      return
      end subroutine
!
      subroutine dader(idif,ina,inc) bind(C, name="dader_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) idif, ina, inc

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoc
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call dadert(idif,ina,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dadert(idif,ina,inc)
      endif

      return
      end subroutine

      subroutine dadert(idif,ina,inc)
      implicit none
      integer ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,    &
     &illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,jj
      integer(8) i, ipoa, ipoc
      double precision rr,x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      integer jd(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
!         PRINT*,'ERROR, DADER CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DADER1',1)
!         stop 666
        do i=1,lnv
          jd(i)=0
        enddo
        jd(idif)=1
        call dapek(ina,jd,rr)
        call dacon(inc,rr)
        return
      endif
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
        ider1  = 0
        ider1s = 0
        ider2  = idif-(nvmax+1)/2
        ider2s = 1
        do jj=1,ider2-1
          ider2s = ider2s*ibase
        enddo
        xdivi  = ider2s*ibase
      else
        ider1  = idif
        ider1s = 1
        do jj=1,ider1-1
          ider1s = ider1s*ibase
        enddo
        ider2  = 0
        ider2s = 0
        xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do i=ipoa,ipoa+illa-1
!
        if(ider1.eq.0) then
          iee = i2(i)
        else
          iee = i1(i)
        endif
!
        x = iee/xdivi
        ifac = int(ibase*(x-int(x+epsmac)+epsmac))
!
        if(ifac.eq.0) goto 100
!
        ic = ic + 1
        cc(ic) = cc(i)*ifac
        i1(ic) = i1(i) - ider1s
        i2(ic) = i2(i) - ider2s
!
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DADER '
        call dadeb(31,'ERR DADER2',1)
      endif
!
      return
      end subroutine
!
      subroutine dapoi(ina,inb,inc,n)
      implicit none
      integer i,ina,inb,inc,n
!     *******************************
!
!     THIS SUBROUTINE COMPUTES THE POISSON BRACKET OF THE VECTORS A AND
!     B AND STORES THE RESULT IN C. N IS THE DEGREE OF FREEDOM OF THE SYSTEM.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer is(4)
!
      is(1) = 0
      is(2) = 0
      is(3) = 0
      is(4) = 0
      call daall(is,4,'$$DAPOI $$',nomax,nvmax)
!
!
      do i=1,n
!
        call dader(2*i-1,ina,is(1))
        call dader(2*i,  inb,is(2))
        call damul(is(1),is(2),is(3))
        call daadd(is(4),is(3),is(1))
        call dacop(is(1),is(4))
!
        call dader(2*i,  ina,is(1))
        call dader(2*i-1,inb,is(2))
        call damul(is(1),is(2),is(3))
        call dasub(is(4),is(3),is(1))
        call dacop(is(1),is(4))
!
      enddo

      call dacop(is(4),inc)
!
      call dadal(is,4)
!
      return
      end subroutine
!
      subroutine dacfuR(ina,fun,inc)
      implicit none
      integer illc,ilmc,ina,inc,incc,inoc,invc
      integer(8) ipoc
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call dacfuRt(ina,fun,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dacfuRt(ina,fun,inc)
      endif

      return
      end subroutine

      subroutine dacfuRt(ina,fun,inc)
      implicit none
      integer i,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,        &
     &j
      integer(8) ia, ic, ipoa, ipoc
      double precision cfac,rr
      double COMPLEX fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
          j(i)=0
        enddo
        call dapek(ina,j,rr)
        cfac = DREAL(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        do i=1,lnv
          j(i)=1
          call dapek(ina,j,rr)
          cfac = DREAL(fun(j))
          rr=cfac*rr
          call dapok(inc,j,rr)
          j(i)=0
        enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
!         stop 667
        return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
        call dancd(i1(ia),i2(ia),j)
        cfac = DREAL(fun(j))
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
        if(dabs(cfac*cc(ia)).lt.eps.or.dabs(cc(ia)).lt.eps) goto 100
!
        ic = ic + 1
        cc(ic) = cc(ia)*cfac
        i1(ic) = i1(ia)
        i2(ic) = i2(ia)
!
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DACFU '
        call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end subroutine
!
      subroutine dacfu(ina,fun,inc) bind(C, name="dacfu_")
      use iso_c_binding, only: c_char, c_long, c_double, c_funptr
      implicit none
      integer(c_long) ina, inc

      abstract interface
        function fun(j) bind(C)
          import :: c_long, c_double
          real(c_long), intent(in) :: j(*)
          real(c_double) :: fun
        end function
      end interface

      integer illc,ilmc,incc,inoc,invc
      integer(8) ipoc
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call dacfut(ina,fun,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dacfut(ina,fun,inc)
      endif

      return
      end subroutine




      subroutine dacfuI(ina,fun,inc)
      implicit none
      integer illc,ilmc,ina,inc,incc,inoc,invc
      integer(8) ipoc
      double complex fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call dacfuIt(ina,fun,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call dacfuIt(ina,fun,inc)
      endif

      return
      end subroutine

      subroutine dacfuIt(ina,fun,inc)
      implicit none
      integer i,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,        &
     &j
      integer(8) ia, ic, ipoa, ipoc
      double precision cfac,rr
      double COMPLEX fun
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
          j(i)=0
        enddo
        call dapek(ina,j,rr)
        cfac = DIMAG(fun(j))
        rr=cfac*rr
        call dapok(inc,j,rr)
        do i=1,lnv
          j(i)=1
          call dapek(ina,j,rr)
          cfac = DIMAG(fun(j))
          rr=cfac*rr
          call dapok(inc,j,rr)
          j(i)=0
        enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
!         stop 667
        return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
        call dancd(i1(ia),i2(ia),j)
        cfac = DIMAG(fun(j))
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
        if(dabs(cfac*cc(ia)).lt.eps.or.dabs(cc(ia)).lt.eps) goto 100
!
        ic = ic + 1
        cc(ic) = cc(ia)*cfac
        i1(ic) = i1(ia)
        i2(ic) = i2(ia)
!
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DACFU '
        call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end subroutine
!

      subroutine dacfut(ina,fun,inc)
      implicit none
      integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,  &
     &j
      integer(8) ipoa, ipoc
      double precision cfac,fun,rr
      external fun
!     *****************************
!
!     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE PRECISION FUNCTION
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension j(lnv)
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
      if(nomax.eq.1) then
        do i=1,lnv
          j(i)=0
        enddo
        call dapek(ina,j,rr)
        cfac = fun(j)
        rr=cfac*rr
        call dapok(inc,j,rr)
        do i=1,lnv
          j(i)=1
          call dapek(ina,j,rr)
          cfac = fun(j)
          rr=cfac*rr
          call dapok(inc,j,rr)
          j(i)=0
        enddo
!         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
!        CALL DADEB(31,'ERR DACFU ',1)
!         stop 667
        return
      endif
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      do ia=ipoa,ipoa+illa-1
!
        call dancd(i1(ia),i2(ia),j)
        cfac = fun(j)
!      IF(dABS(CFAC).LT.EPS) GOTO 100
!      IF(dABS(CFAC*CC(IA)).LT.EPS) GOTO 100
        if(dabs(cfac*cc(ia)).lt.eps.or.dabs(cc(ia)).lt.eps) goto 100
!
        ic = ic + 1
        cc(ic) = cc(ia)*cfac
        i1(ic) = i1(ia)
        i2(ic) = i2(ia)
!
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DACFU '
        call dadeb(31,'ERR DACFU ',1)
      endif
!
      return
      end subroutine
!

      subroutine daimp(r, ic1, ic2, ina) bind(C, name="daimp_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ic1, ic2, ina
      real(c_double) r

      integer    lh
      integer(8) i, ic
*     **************************
*
*     THIS SUBROUTINE "IMPORTS" THE ARRAY H WITH LENGTH LH AND PUTS ITS
*     ENTRIES INTO THE DA VECTOR A
*
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9

      DIMENSION r(nmmax+1), ic1(nmmax), ic2(nmmax)

      lh = nint(r(1))
      call daclr(1)
      ic = idapo(ina)

      do i = 1, lh
        cc(ic) = r(i+1)
        i1(ic) = ic1(i)
        i2(ic) = ic2(i)
        ic = ic + 1
      enddo

      idall(ina) = lh

      return
      end subroutine

      subroutine daexp(ina, r, ic1, ic2, name) bind(C, name="daexp_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long)   ina, ic1, ic2
      real(c_double)    r
      character(c_char) name(*)

      integer    i, lh, k, jj
      integer(8) ic
*     **************************
*
*     THIS SUBROUTINE "EXPORTS" THE DA VACTOR A AND STORES ITS COEFFICIENTS
*     IN THE ARRAY H WITH LENGTH LH
*
      integer j
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
*
      DIMENSION r(nmmax+1), ic1(nmmax), ic2(nmmax), jj(lnv)

!      call dapri77(ina, 6)
      ic = idapo(ina)
      lh = idall(ina)
      r(1) = lh
!      name = daname(ina)
      do k=1, 10
         name(k) = daname(ina)(k:k)
      enddo

      do k = 1, lnv
         jj(k) = 0
      enddo
      do i = 1, lh
         r(i+1) = cc(ic)
         if (nomax .eq. 1) then
            if (i > 1) jj(i-1) = 1;
            call hash(nomax, lnv, jj, ic1(i), ic2(i))
            if (i > 1) jj(i-1) = 0;
         else
            ic1(i) = i1(ic)
            ic2(i) = i2(ic)
         endif
         ic = ic + 1
      enddo

      return
      end subroutine

      subroutine dapri(ina,iunit) bind(C, name="dapri_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, iunit

      integer iii,illa,ilma,inoa,inva,ioa,iout,ipoa,j,k
      integer(8) i, ii
!     ***************************
!       Frank
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
!
      if(ina.lt.1.or.ina.gt.nda) then
        print*,'ERROR IN DAPRI, INA = ',ina
!        X = SQRT(-ONE)
!        PRINT*,X
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')                       &
     &daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,            &
     &'***********'//'**********************************'
!
      iout = 0
      ioa = 0

      if(inva.eq.0) then
        write(iunit,'(A)')                                              &
     &  '    I  VALUE  '
        do i = ipoa,ipoa+illa-1
          write(iunit,'(I6,2X,G20.14)') i-ipoa, cc(i)
        enddo
      elseif(nomax.eq.1) then
        if(illa.ne.0) write(iunit,'(A)')                                &
     &  '    I  COEFFICIENT          ORDER   EXPONENTS'
        if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
        do i=1,illa
          do k=1,inva
            j(k)=0
          enddo
          iout=iout+1
          if(i.ne.1) then
            j(i-1)=1
            ioa=1
          endif
          if(abs(cc(ipoa+i-1)).gt.eps) then
            write(iunit,'(I6,2X,E21.14,I5,4X,18(2I2,1X))')              &
     &            iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
!            write(iunit,*) cc(ipoa+i-1)
          endif
        enddo
      else
        if(illa.ne.0) write(iunit,'(A)')                                &
     &  '    I  COEFFICIENT          ORDER   EXPONENTS'
        if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
        do ioa = 0,inoa
          do ii=ipoa,ipoa+illa-1
            if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
            call dancd(i1(ii),i2(ii),j)
!ETIENNE
            if(abs(cc(ii)).gt.eps) then
!ETIENNE
              iout = iout+1
              write(iunit,'(I6,2X,E21.14,I5,4X,18(2I2,1X))')            &
     &        iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!ETIENNE
!              write(iunit,*) cc(ii)
            endif
!ETIENNE
!
 100        continue
          enddo
        enddo
!
      endif

      write(iunit,'(A)') '                                      '
!
      return
      end subroutine

      subroutine dapriold(ina,iunit)
      implicit none
!     ***************************
!       Frank
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------
!
      integer iii,illa,ilma,ina,inoa,inva,ioa,iout,iunit,j,k
      integer(8) i, ii, ipoa
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      if(ina.lt.1.or.ina.gt.nda) then
         print*,'ERROR IN DAPRI, INA = ',ina
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)')
     &     daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,
     &     '*********************************************'
!
      iout = 0
      ioa = 0

      if(inva.eq.0) then
         write(iunit,'(A)') '    I  VALUE  '
         do i = ipoa,ipoa+illa-1
            write(iunit,'(I6,2X,G20.14)') i-ipoa, cc(i)
         enddo
      elseif(nomax.eq.1) then
         if(illa.ne.0) write(iunit,'(A)')
     &        '    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
         do i=1,illa
            do k=1,inva
               j(k)=0
            enddo
            iout=iout+1
            if(i.ne.1) then
               j(i-1)=1
               ioa=1
            endif
            write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     &           iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
            write(iunit,*) cc(ipoa+i-1)
         enddo
      else
         if(illa.ne.0) write(iunit,'(A)')
     &        '    I  COEFFICIENT          ORDER   EXPONENTS'
         if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
         do ioa = 0,inoa
            do ii=ipoa,ipoa+illa-1
               if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
               call dancd(i1(ii),i2(ii),j)
!ETIENNE
               if(abs(cc(ii)).gt.eps) then
!ETIENNE
                  iout = iout+1
                  write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')
     &                 iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!ETIENNE
                  write(iunit,*) cc(ii)
               endif
!ETIENNE
!
 100           continue
            enddo
         enddo
!
      endif

      write(iunit,'(A)') '                                      '
!
      return
      end subroutine

      subroutine dapri77(ina,iunit)
      implicit none
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j
      character c10*10,k10*10
!     ***************************
!       Etienne
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      if(iunit.eq.0) return
!
      if(ina.lt.1.or.ina.gt.nda) then
        write(6,*)'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
!      WRITE(IUNIT,*) INA, ' in dapri ', DANAME(INA)
!      WRITE(6,*) INA, ' in dapri ', DANAME(INA)
! 611  WRITE(6,*) ' MORE '
!        READ(5,*) MORE
!        IF(MORE.GT.0) THEN
!        WRITE(6,*) MORE,' ',DANAME(MORE)
!        GOTO 611
!        ENDIF
      write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)')                  &
     &daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,            &
     &'***********'//'**********************************'
!
      if(illa.ne.0) write(iunit,'(A)')                                  &
     &'    I  COEFFICIENT          ORDER   EXPONENTS'
      if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
!
      c10='      NO ='
      k10='      NV ='

      write(iunit,'(A10,I6,A10,I6)') c10,inoa,k10,inva

      iout = 0
!
!      DO 100 IOA = 0,INOA
      do ioa = 0,nocut
        do ii=ipoa,ipoa+illa-1
          if(nomax.ne.1) then
            if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          endif
!ETIENNE
          if(dabs(cc(ii)).gt.eps) then
!ETIENNE

            if(nomax.ne.1) then
              call dancd(i1(ii),i2(ii),j)
              iout = iout+1
            else
              if(ii.eq.ipoa.and.ioa.eq.1) goto 100
              if(ii.gt.ipoa.and.ioa.eq.0) goto 100
              do i=1,lnv
                j(i)=0
              enddo
              if(ii.ne.ipoa) j(ii-ipoa)=1
              iout = iout+1
            endif
!

!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
            if(dabs(cc(ii)).gt.eps) then
              if(eps.gt.1.e-37) then
                write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
              else
                write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
              endif
            endif
 501        format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503        format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502        format(' ', i5,1x,g23.16,1x,100(1x,i2))

          endif
!ETIENNE
!
 100      continue
        enddo
      enddo
!
      do i=1,lnv
        j(i)=0
      enddo

      if(iout.eq.0) iout=1

      write(iunit,502) -iout,0.d0,(j(i),i=1,inva)
!
      return
      end subroutine

      subroutine dashift(ina,inc,ishift)
      implicit none
      integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,j
      double precision c
!       Frank
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      integer inb,ishift,ich,ik,jd(lnv),inc
!-----------------------------------------------------------------------------3
!
!

      inb=0
      if(ina.lt.1.or.ina.gt.nda) then
        write(6,*)'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
      call daall1(inb,'$$DAJUNK$$',inoa,inva)

!      WRITE(IUNIT,*) INA, ' in dapri ', DANAME(INA)
!      WRITE(6,*) INA, ' in dapri ', DANAME(INA)
! 611  WRITE(6,*) ' MORE '
!        READ(5,*) MORE
!        IF(MORE.GT.0) THEN
!        WRITE(6,*) MORE,' ',DANAME(MORE)
!        GOTO 611
!        ENDIF
      iout = 0
!
!      DO 100 IOA = 0,INOA
      do ioa = 0,nocut
        do ii=ipoa,ipoa+illa-1
          if(nomax.ne.1) then
            if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
          endif
!ETIENNE
          if(dabs(cc(ii)).gt.eps) then
!ETIENNE

            if(nomax.ne.1) then
              call dancd(i1(ii),i2(ii),j)
              iout = iout+1
            else
              if(ii.eq.ipoa.and.ioa.eq.1) goto 100
              if(ii.gt.ipoa.and.ioa.eq.0) goto 100
              do i=1,lnv
                j(i)=0
              enddo
              if(ii.ne.ipoa) j(ii-ipoa)=1
              iout = iout+1
            endif
!

!      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
            if(dabs(cc(ii)).gt.eps) then
              if(eps.gt.1.e-37) then
!       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
                ich=1
                do ik=1,ishift
                  if(j(ik).ne.0) ich=0
                enddo
                if(ich.eq.1) then
                  do ik=1,lnv
                    jd(ik)=0
                  enddo
                  do ik=ishift+1,lnv
                    jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                  enddo
                endif
                call dapok(inb,jd,cc(ii))
              else
!       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
                ich=1
                do ik=1,ishift
                  if(j(ik).ne.0) ich=0
                enddo
                if(ich.eq.1) then
                  do ik=1,lnv
                    jd(ik)=0
                  enddo
                  do ik=ishift+1,lnv
                    jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                  enddo
                endif
                call dapok(inb,jd,cc(ii))
              endif
            endif
 501        format(' ', i3,1x,g23.16,1x,100(1x,i2))
 503        format(' ', i3,1x,g23.16,1x,100(1x,i2))
 502        format(' ', i5,1x,g23.16,1x,100(1x,i2))

          endif
!ETIENNE
!
 100      continue
        enddo
      enddo
!
      do i=1,lnv
        j(i)=0
      enddo

      call dacop(inb,inc)
      call dadal1(inb)
!
      return
      end subroutine

      subroutine darea(ina,iunit) bind(C, name="darea_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, iunit

      integer i,iche,ii,ii1,ii2,iin,illa,ilma,inoa,inva,io,io1,         &
     &iwarin,iwarno,iwarnv,j,nno
      integer(8) ic, ipoa
      double precision c
!       Frank
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      character c10*10
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
        print*,'ERROR IN DAREA, INA = ',ina
!        X = SQRT(-ONE)
!        PRINT*,X
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do i=1,lnv
        j(i) = 0
      enddo
!
      call daclr(1)
!
      ic = 0
!
      iwarno = 0
      iwarnv = 0
      iwarin = 0
!
      read(iunit,'(A10)') c10
      read(iunit,'(18X,I4)') nno
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10

!
!
      iin = 0
!
  10  continue
      iin = iin + 1
      read(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                     &
     &ii,c,io,(j(i),i=1,inva)
!
      if(ii.eq.0) goto 20
!ETIENNE
      read(iunit,*) c
!ETIENNE
      if(ii.ne.iin) then
        iwarin = 1
      endif
      io1 = 0
      do i=1,inva
        io1 = io1 + j(i)
      enddo
!
      if(io1.ne.io) then
        iwarnv = 1
        goto 10
      endif
      if(io.gt.inoa) then
!        IF(IWARNO.EQ.0) PRINT*,'WARNING IN DAREA, FILE ',
!    *              'CONTAINS HIGHER ORDERS THAN VECTOR '
        iwarno = 1
        goto 10
      endif
!
      if(nomax.ne.1) then
        ic = ic + 1
        call dadcd(j,ii1,ii2)
        ic = ia1(ii1) + ia2(ii2)
        cc(ic) = c
        goto 10
      else
        iche=0
        do i=1,inva
          if(j(i).eq.1) iche=i
        enddo
        cc(ipoa+iche)=c
        goto 10
      endif

!
  20  continue
!
      if(nomax.ne.1) call dapac(ina)
!
      return
      end subroutine
!FF
!
      subroutine darea77(ina,iunit) bind(C, name="darea77_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) ina, iunit

      integer i,ic,iche,ii,ii1,ii2,iin,illa,ilma,inoa,inva,ipoa,        &
     &j,k,nojoh,nvjoh
      double precision c
!     ***************************
!     Etienne
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      character c10*10,k10*10
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
        write(6,*)'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do i=1,lnv
        j(i) = 0
      enddo
!
      call daclr(1)
!
!
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10,I6,A10,I6)') c10,nojoh,k10,nvjoh
!
      iin = 0
!
  10  continue
      iin = iin + 1
      read(iunit,*) ii,c,(j(k),k=1,nvjoh)
      if(ii.lt.0) goto 20

      do i=inva+1,nvjoh
        if(j(i).ne.0) goto 10
      enddo
      iche=0
      do i=1,inva
        iche=iche+j(i)
      enddo
      if(iche.gt.nomax) goto 10
      if(nomax.ne.1) then
        call dadcd(j,ii1,ii2)
        ic = ia1(ii1) + ia2(ii2)
        cc(ic) = c
      else
        iche=0
        do i=1,inva
          if(j(i).eq.1) iche=i
        enddo
        cc(ipoa+iche)=c
      endif
      goto 10
!
  20  continue
!
      if(nomax.ne.1) call dapac(ina)
!
      return
      end subroutine

      subroutine dadeb(iunit,c,istop)
      implicit none
      integer istop,iunit
!     *******************************
!
!     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL. IT PRINTS ALL
!     NONZERO INFORMATION IN THE COMMON BLOCKS AND ALL DA  VECTORS.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      character c*10
!
!etienne

      write(6,*) '  ',c
      call flush(6)
      stop
!      write(6,*) daname(lda+1)
      end subroutine
!
!
!
      subroutine danum(no,nv,numda)
      implicit none
      integer i,mm,no,numda,nv
!     *****************************
!
!     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
!     ORDER NO AND NUMBER OF VARIABLES NV
!
      numda = 1
      mm = max(nv,no)
!
      do i=1,min(nv,no)
        numda = (numda*(mm+i))/i
      enddo
!
      return
      end subroutine
!
      subroutine dainf(inc,inoc,invc,ipoc,ilmc,illc)
      implicit none
      integer illc,ilmc,inc,inoc,invc
      integer(8) ipoc
!     **********************************************
!
!     THIS SUBROUTINE SEARCHES THE NUMBER OF DA VECTOR C
!     AND RETURS THE INFORMATION IN COMMON DA
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      if(inc.ge.1.and.inc.le.nda) then
        inoc = idano(inc)
        invc = idanv(inc)
        ipoc = idapo(inc)
        ilmc = idalm(inc)
        illc = idall(inc)
        return
      endif
!
      write(6,*) 'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
      call dadeb(31,'ERR DAINF ',1)
!
      return
      end subroutine
!
      subroutine dapac(inc)
      implicit none
      integer illc,ilmc,inc,inoc,invc
      integer(8) i, ic, ipoc
      double precision ccc
!     ************************
!
!     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR 1
!     INTO THE VECTOR INC.
!     INVERSE IS DAUNP.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      ipoc = idapo(inc)

!
      ic = ipoc - 1
!
      do i=1,nmmax
        ccc = cc(i)
        if(dabs(ccc).lt.eps) goto 100
        ic = ic + 1
        cc(ic) = ccc
        i1(ic) = ie1(i)
        i2(ic) = ie2(i)
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DAPAC '
        call dadeb(31,'ERR DAPAC ',1)
      endif
!
      return
      end subroutine
!
!
      subroutine dachk(ina,inoa,inva, inb,inob,invb, inc,inoc,invc)
      implicit none
      integer ierr,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,invsum,lsw
!     *************************************************************
!
!     THIS SUBROUTINE CHECKS IF THE VECTORS A, B AND C
!     HAVE COMPATIBLE ATTRIBUTES
!
      parameter(lsw=1)
!
      if(lsw.eq.1) return
!
      ierr = 0
!
!     CASE OF A UNARY OPERATION
!     *************************
!
      if(inob.eq.-1.and.invb.eq.-1) then
        invsum = inva + invc
        if(invsum.eq.0) then
          if(inoa.gt.inoc) ierr = 1
        elseif(invsum.eq.1) then
          ierr = 1
        else
          if(inoa.gt.inoc.or.inva.gt.invc) ierr = 1
        endif
        if(ierr.eq.1) then
          write(6,*)'ERROR IN DACHK, ',ina,' AND ',inc,                 &
     &    ' ARE INCOMPATIBLE',inoa,inva,inoc,invc
          call dadeb(31,'ERR DACHK1',1)
        endif
!
!     CASE OF A BINARY OPERATION
!     **************************
!
      else
        invsum = inva + invb + invc
        if(invsum.eq.0) then
          if(inoa.gt.inoc.or.inob.gt.inoc) ierr = 1
        elseif(invsum.eq.1.or.invsum.eq.2) then
          ierr = 1
        else
          if(inoa.gt.inoc.or.inob.gt.inoc.or.                           &
     &    inva.gt.invc.or.invb.gt.invc) ierr = 1
        endif
        if(ierr.eq.1) then
          write(6,*)'ERROR IN DACHK, ',ina,',',inb,' AND ',inc,         &
     &    ' ARE INCOMPATIBLE'
          call dadeb(31,'ERR DACHK2',1)
        endif
      endif
!
      return
      end subroutine
!
      subroutine damch(iaa,ia)
      implicit none
      integer i,ia,iaa,illa,ilma,ino1,inoi,inv1,invi
      integer(8) ipoa
!     ************************
!
!     THIS SUBROUTINE CHECKS IF THE IA VECTORS IN THE MATRIX IA HAVE
!     IDENTICAL ATTRIBUTES.
!
      dimension iaa(*)
!
      call dainf(iaa(1),ino1,inv1,ipoa,ilma,illa)
!
      do i=2,ia
        call dainf(iaa(i),inoi,invi,ipoa,ilma,illa)
        if(ino1.ne.inoi.or.inv1.ne.invi) then
          write(6,*)'ERROR IN DAMCH, VECTORS ',iaa(1),' AND ',iaa(i),   &
     &    ' ARE INCOMPATIBLE '
          stop
        endif
      enddo
!
      return
      end subroutine
!
      subroutine dadcd(jj,ic1,ic2)
      implicit none
      integer i,ibase,ic1,ic2,isplit
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv)
      ibase = nomax + 1
      isplit = (nvmax+1)/2
      ic1 = 0
      ic2 = 0
!
      do i=nvmax,isplit+1,-1
        ic2 = ic2*ibase + jj(i)
      enddo
!
      do i=isplit,1,-1
        ic1 = ic1*ibase + jj(i)
      enddo
!
      return
      end subroutine
!
      subroutine dancd(ic1,ic2,jj)
      implicit none
      integer i,ibase,ic,ic1,ic2,isplit
      double precision x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(*)
      ibase = nomax + 1
      isplit = (nvmax+1)/2
!
      ic = ic1
      do i=1,isplit
        x  = dble(ic)/dble(ibase)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      enddo
!
      ic = ic2
      do i=isplit+1,nvmax
        x  = dble(ic)/dble(ibase)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      enddo
!
      do i=nvmax+1,lnv
        jj(i) = 0
      enddo
!
      return
      end subroutine

!ETIENNE
      subroutine datra(idif,ina,inc) bind(C, name="datra_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) idif, ina, inc

      integer ibase,ider1,ider1s,ider2,ider2s,iee,ifac,illa,            &
     &illc,ilma,ilmc,inoa,inoc,inva,invc,jj
      integer(8) i, ic, ipoa, ipoc
      double precision x,xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE PSEUDO DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!     dx^n/dx= x^(n-1)
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
!       CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!

      if(nomax.eq.1) then
        call dader(idif,ina,inc)
        return
      endif
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
        ider1  = 0
        ider1s = 0
        ider2  = idif-(nvmax+1)/2
        ider2s = 1
        do jj=1,ider2-1
          ider2s = ider2s*ibase
        enddo
        xdivi  = ider2s*ibase
      else
        ider1  = idif
        ider1s = 1
        do jj=1,ider1-1
          ider1s = ider1s*ibase
        enddo
        ider2  = 0
        ider2s = 0
        xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do i=ipoa,ipoa+illa-1
!
        if(ider1.eq.0) then
          iee = i2(i)
        else
          iee = i1(i)
        endif
!
        x = iee/xdivi
        ifac = int(ibase*(x-int(x+epsmac)+epsmac))
!
        if(ifac.eq.0) goto 100
!
!etienne      IFAC = INT(IBASE*(X-INT(X)+1.D-8))
!
!etienne      IF(IFAC.EQ.0) GOTO 100
!
        ic = ic + 1
        cc(ic) = cc(i)
        i1(ic) = i1(i) - ider1s
        i2(ic) = i2(i) - ider2s
!
 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DADTRA'
        call dadeb(111,'ERR DADTRA',1)
      endif
!
      return
      end subroutine

      subroutine etred(no1,nv1,ic1,ic2,no2,nv2,i11,i21)
      implicit none
      integer i,i11,i21,ic,ic1,ic2,no1,no2,nv1,nv2
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jj(lnv)

      if(nv1.gt.lnv.or.nv2.gt.lnv) then
        write(6,*) ' ERROR IN RECODING '
        stop123
      endif
      call dehash(no1,nv1,ic1,ic2,jj)
      ic=0

      do i=nv1+1,nv2
        if(jj(i).gt.0) then
          i11=-1
          i21=0
          return
        endif
      enddo

      do i=1,nv2
        ic=ic+jj(i)
      enddo

      if(ic.gt.no2) then
        i11=-1
        i21=0
        return
      endif

      call hash(no2,nv2,jj,i11,i21)

!
      return
      end subroutine

      subroutine hash(no1,nv1,jj,ic1,ic2) bind(C, name="hash_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) no1, nv1, jj(*), ic1, ic2

      integer i,ibase,isplit
!     ****************************
!
!     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"

      ibase = no1 + 1
      isplit = (nv1+1)/2
      ic1 = 0
      ic2 = 0
!
      do i=nv1,isplit+1,-1
        ic2 = ic2*ibase + jj(i)
      enddo
!
      do i=isplit,1,-1
        ic1 = ic1*ibase + jj(i)
      enddo
!
      return
      end subroutine
!
      subroutine dehash(no1,nv1,ic1,ic2,jj) bind(C, name="dehash_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      integer(c_long) no1, nv1, ic1, ic2, jj(*)

      integer i,ibase,ic,isplit
      double precision x
!     ****************************
!
!     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES I1,I2.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!

      epsmac=1.e-7
      ibase = no1 + 1
      isplit = (nv1+1)/2
!
      ic = ic1
      do i=1,isplit
        x  = dble(ic)/dble(ibase)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      enddo
!
      ic = ic2
      do i=isplit+1,nv1
        x  = dble(ic)/dble(ibase)
        ic = int(x+epsmac)
        jj(i) = nint(ibase*(x-ic))
      enddo
!
      return
      end subroutine

      subroutine daswap(j1,j2,inb)
      implicit none
      integer ic1,ic2,illb,ilmb,inb,inob,invb,j1,j2,jj,k1,k2
      integer(8) ia, ic, ipob
!     *************************
!
!     SWAP A DA VECTOR
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension jj(lnv)

      call dainf(inb,inob,invb,ipob,ilmb,illb)

      call daclr(1)
!

      do ia = ipob,ipob+illb-1

        call dehash(nomax,nvmax,i1(ia),i2(ia),jj)
        k1=jj(j1)
        k2=jj(j2)
        jj(j1)=k2
        jj(j2)=k1
        call hash(nomax,nvmax,jj,ic1,ic2)

        ic=ia1(ic1)+ia2(ic2)

        cc(ic) = cc(ia)
        i1(ic) = ic1
        i2(ic) = ic2
!
      enddo
!
!
      call dapac(inb)
      return
      end subroutine

      subroutine dagauss(ina,inb,nd2,anorm)
      implicit none
      integer i,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,        &
     &ja,jb,nd2
      integer(8) ia, ib, ipoa, ipob
      double precision anorm,gau
!     ***************************
!
!     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension ja(lnv),jb(lnv)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
      anorm = 0.d0

      do ia=ipoa,ipoa+illa-1
        do ib=ipob,ipob+illb-1
          call dancd(i1(ia),i2(ia),ja)
          call dancd(i1(ib),i2(ib),jb)
          gau=1.d0
          do i=1,nd2
            gau= facint(ja(i)+jb(i))*gau
          enddo
          anorm = anorm + cc(ia)*cc(ib)*gau
        enddo
      enddo

!
      return
      end subroutine

      subroutine daran(ina,cm,xran)
      implicit none
      integer illa,ilma,ina,inoa,inva
      integer(8) i, ipoa
      double precision bran,cm,xran
!     ************************
!
!     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
!     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
!     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
!     ABS(CM) IS THE FILLING FACTOR
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
!
      if(inva.eq.0.or.nomax.eq.1) then
        do i=ipoa,ipoa+ilma-1
          if(cm.gt.0d0) then
            cc(i) = bran(xran)
            if(cc(i).gt.cm) cc(i) = 0d0
          elseif(cm.lt.0d0) then
            cc(i) = int(1+10*bran(xran))
            if(cc(i).gt.-1d1*cm) cc(i) = 0d0
          endif
        enddo
        idall(ina) = idalm(ina)
        return
      endif
!
      if(inoa.ne.nomax.or.inva.ne.nvmax) then
        write(6,*)'ERROR IN DARAN, ONLY VECTORS WITH NO = NOMAX AND'    &
     &  //' NV = NVMAX ALLOWED'
        call dadeb(31,'ERR DARAN1',1)
      endif
!
      call daclr(1)
!
      do i=1,nmmax
        if(cm.gt.0.d0) then
          cc(i) = bran(xran)
          if(cc(i).gt.cm) cc(i) = 0.d0
        elseif(cm.lt.0.d0) then
          cc(i) = int(1+10*bran(xran))
          if(cc(i).gt.-10.d0*cm) cc(i) = 0.d0
        else
          write(6,*)'ERROR IN ROUTINE DARAN'
          call dadeb(31,'ERR DARAN2',1)
        endif
      enddo
!
      call dapac(ina)
!
      return
      end subroutine
!
      real(c_double) function bran(xran) bind(C, name="bran_")
      use iso_c_binding, only: c_char, c_long, c_double
      implicit none
      real(c_double) xran

!     ************************************
!
!     VERY SIMPLE RANDOM NUMBER GENERATOR
!
      xran = xran + 10.d0
      if(xran.gt.1.d4) xran = xran - 9999.12345
      bran = dabs(sin(xran))
      bran = 10*bran
      bran = bran - int(bran)
!      IF(BRAN.LT. .1D0) BRAN = BRAN + .1D0
!
      return
      end function

      subroutine danorm2(ina,inc)
      implicit none
      integer illc,ilmc,ina,inc,incc,inoc,invc
      integer(8) ipoc
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call danorm2t(ina,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call danorm2t(ina,inc)
      endif

      return
      end subroutine

      subroutine danorm2t(ina,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8) ia, ib, ipoa, ipob
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
!
      ib = ipob - 1
!
      do ia=ipoa,ipoa+illa-1
!
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
        ib = ib + 1
        cc(ib) = cc(ia)**2
        i1(ib) = i1(ia)
        i2(ib) = i2(ia)
!
 100    continue
      enddo
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
        write(6,*)'ERROR IN DANORM'
        call dadeb(31,'ERR DANOR1',1)
      endif
!
      return
      end subroutine

      subroutine danormr(ina,inc)
      implicit none
      integer illc,ilmc,ina,inc,incc,inoc,invc
      integer(8) ipoc
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!

      if(ina.eq.inc) then
        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
        incc=0
        call daall1(incc,'$$DAJUNK$$',inoc,invc)
        call danormrt(ina,incc)
        call dacop(incc,inc)
        call dadal1(incc)
      else
        call danormrt(ina,inc)
      endif

      return
      end subroutine

      subroutine danormrt(ina,inb)
      implicit none
      integer illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb
      integer(8) ia, ib, ipoa, ipob
!     ******************************
!
!     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
!     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
!     THE DA VECTOR DENOTED WITH THE INTEGER E.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inb,inob,invb,ipob,ilmb,illb)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
!
!
      ib = ipob - 1
!
      do ia=ipoa,ipoa+illa-1
!
        if(ieo(ia1(i1(ia))+ia2(i2(ia))).gt.nocut) goto 100
        ib = ib + 1
        cc(ib) = dsqrt(cc(ia))
        i1(ib) = i1(ia)
        i2(ib) = i2(ia)
!
 100    continue
      enddo
!
      idall(inb) = ib-ipob+1
      if(idall(inb).gt.idalm(inb)) then
        write(6,*)'ERROR IN DANORM '
        call dadeb(31,'ERR DANOR2',1)
      endif
!
      return
      end subroutine
      subroutine dakey(c)
      implicit none
      character c*(*)
!

      return
!
      end subroutine
! ANFANG UNTERPROGRAMM
      subroutine dapri6(ina,result,ien,i56)
      implicit none
      integer i,i56,ien,ihp,illa,ilma,ina,inoa,inva,ioa,iout,           &
     &j
      integer(8) ipoa, ii
      double precision result
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     RESULT CONTAINS THE (IEN-1) TH DERIVATIVE TO ENERGY AFTER EXECUTION
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
      dimension j(lnv)
      result=0.
      if(ina.lt.1.or.ina.gt.nda) then
        print*,'ERROR IN DAPRI6, INA = ',ina
        stop
      endif
      inoa=idano(ina)
      inva=idanv(ina)
      ipoa=idapo(ina)
      ilma=idalm(ina)
      illa=idall(ina)
      iout=0
      if(nomax.eq.1) then
        do i=1,illa
          if(ien.eq.1) then
            if(i-1.ne.0) goto 90
            result=cc(ipoa+i-1)
            return
          endif
          if(ien.eq.2) then
            if(i-1.ne.i56) goto 90
            result=cc(ipoa+i-1)
            return
          endif
 90       continue
        enddo
      else
        do ioa=0,inoa
          do ii=ipoa,ipoa+illa-1
            if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
            iout=iout+1
            call dancd(i1(ii),i2(ii),j)
            if(i56.eq.6) then
              do ihp=1,5
                if(j(ihp).ne.0) goto 100
              enddo
              if(j(6).eq.(ien-1)) then
                result=cc(ii)
                return
              endif
            else if(i56.eq.5) then
              do ihp=1,4
                if(j(ihp).ne.0) goto 100
              enddo
              if(j(5).eq.(ien-1)) then
                result=cc(ii)
                return
              endif
            else if(i56.eq.4) then
              do ihp=1,3
                if(j(ihp).ne.0) goto 100
              enddo
              if(j(4).eq.(ien-1)) then
                result=cc(ii)
                return
              endif
            endif
 100        continue
          enddo
        enddo
      endif
      return
      end subroutine
! ANFANG UNTERPROGRAMM

      subroutine darea6(ina,zfeld,i56)
      implicit none
      integer i,i56,ii1,ii2,iin,illa,ilma,ina,inoa,inva,io,io1,ip,      &
     &iwarin,iwarno,iwarnv,j
      integer(8) ipoa, ic
      double precision zfeld
!     *************************************
!
!     THIS SUBROUTINE IS FOR REDUCED STORAGE DA VERSION JULY 91
!     I56 SAYS WHETHER THE 5TH OR THE 6TH COORDINATE IS THE ENERGY
!     AND MUST HAVE THE VALUE 5 OR 6 ACCORDINGLY
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      character daname(lda)*10
      common / daname / daname
      dimension zfeld(100)
!-----------------------------------------------------------------------------3
      dimension j(lnv)
      if(ina.lt.1.or.ina.gt.nda) then
        print*,'ERROR IN DAREA6, INA = ',ina
        stop
      endif
      inoa=idano(ina)
      inva=idanv(ina)
      ipoa=idapo(ina)
      ilma=idalm(ina)
      illa=idall(ina)
      do i=1,lnv
        j(i)=0
      enddo
      call daclr(1)
      ic=0
      iwarno=0
      iwarnv=0
      iwarin=0
      iin=0
      if(nomax.eq.1) then
        do i=1,illa
          if (i-1.eq.0) then
            cc(ipoa+i-1)=zfeld(1)
          else if (i-1.eq.i56) then
            cc(ipoa+i-1)=zfeld(2)
          endif
        enddo
        return
      endif
      do ip=1,inva
        j(ip)=0
      enddo
      io=0
  10  continue
      iin=iin+1
      io1=0
      do i=1,inva
        io1=io1+j(i)
      enddo
      if(io1.ne.io) then
        if(iwarnv.eq.0) print*,'WARNING IN DAREA6, FILE ',              &
     &  'CONTAINS MORE VARIABLES THAN VECTOR'
        iwarnv = 1
        goto 10
      endif
      if(io.gt.inoa) then
        iwarno = 1
        goto 10
      endif
      ic = ic + 1
      call dadcd(j,ii1,ii2)
      ic = ia1(ii1) + ia2(ii2)
      cc(ic) = zfeld(io+1)
      j(i56)=j(i56)+1
      io=io+1
      if (io.gt.inoa) goto 20
      goto 10
  20  continue
      call dapac(ina)
      return
      end subroutine
! ANFANG FUNKTION
      double precision function dare(ina)
      implicit none
      integer ii,illa,ilma,ina,inoa,inva,ioa,j,jj
      integer(8) ipoa
!     ***********************************
!     NEW VERSION OF DARE, AUGUST 1992
!     SUPPOSED TO TREAT THE 0TH COMPONENT ACCURATELY
!
!     30.10 1997 E.Mcintosh & F.Schmidt
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      dimension j(lnv)
!-----------------------------------------------------------------------------9
!
!      CALL DAINF(INA,INOA,INVA,IPOA,ILMA,ILLA)
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)

!FRS 30.10.1997
      if(nomax.eq.1) then
        dare = cc(ipoa)
        return
      endif
!FRS 30.10.1997
!FRS March 1997
!      IF(NOMAX.EQ.1) goto 110
!FRS March 1997

      ioa = 0
      do ii=ipoa,ipoa+illa-1
        if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
        call dancd(i1(ii),i2(ii),j)
        do jj=1,inva
          if(j(jj).ne.0) goto 100
        enddo
! 110    continue
        dare = cc(ipoa)
        return
 100    continue
      enddo
      dare = 0d0
      return
      end function
! ANFANG UNTERPROGRAMM

      subroutine daprimax(ina,iunit)
      implicit none
      integer iii,illa,ilma,ina,inoa,inva,ioa,iout,iunit,j
      integer(8) i, ii, ipoa
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!-----------------------------------------------------------------------------9
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      dimension j(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
        print*,'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      iout = 0
      ioa = 0

      if(inva.eq.0) then
        do i = ipoa,ipoa+illa-1
          write(iunit,'(I6,2X,G20.14)') i-ipoa, cc(i)
        enddo
      elseif(nomax.eq.1) then
        do i=1,illa
          iout=iout+1
          if(i.ne.1) then
            j(i-1)=1
            ioa=1
          endif
          write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                &
     &    iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
          write(iunit,*) cc(ipoa+i-1)
        enddo
      else
        iout = 0
        do ioa = 0,inoa
          do ii=ipoa,ipoa+illa-1
            if(ieo(ia1(i1(ii))+ia2(i2(ii))).ne.ioa) goto 100
            call dancd(i1(ii),i2(ii),j)
!ETIENNE
            if(abs(cc(ii)).gt.eps) then
!ETIENNE
              iout = iout+1
              write(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')            &
     &        iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
!ETIENNE
              write(iunit,*) cc(ii)
            endif
!ETIENNE
!
 100        continue
          enddo
        enddo
      endif
!

!     WRITE(IUNIT,'(A)') '                                      '
!
      return
      end subroutine
!FF

!  unknown stuff
      subroutine damono(ina,jd,cfac,istart,inc)
      implicit none
      integer illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,          &
     &istart,jd
      integer(8) ia, ic, ipoa, ipoc
      double precision cfac
!     *****************************
!
!     THIS SUBROUTINE RETURNS THE MONOMIALS ONE BY ONE
!     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
!     RESULT IN C
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      dimension jd(*)
!
      if(ina.eq.inc) then
        write(6,*) ' USE DIFFERENT POWER SERIES IN DAMONO '
        stop999
      endif
      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
      if(istart.eq.0) then
        istart=illa
        return
      endif
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ic = ipoc - 1
!
      ia=ipoa+istart-1
!
      call dancd(i1(ia),i2(ia),jd)

!
      ic = ic + 1
      cc(ic) = cc(ia)
      cfac=cc(ia)
      i1(ic) = i1(ia)
      i2(ic) = i2(ia)
!
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DAMONO'
        call dadeb(31,'ERR DAMONO',1)
      endif
!
      return
      end subroutine
!
!

      subroutine dacycle(ina,ipresent,value,j,illa)
      implicit none
      integer i,illa,ilma,ina,inoa,inva,iout,ipoa,ipresent,j
      integer(8) ii
      double precision value
!     ***************************
!
!     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      dimension j(lnv)
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
!
      if(ina.lt.1.or.ina.gt.nda) then
        write(6,*)'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      if(ina.eq.0) then
        value=0.d0
        illa=0
        do i=1,lnv
          j(i)=0
        enddo
        return
      endif

      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
      iout = 0
      ipresent=1+ipresent
      if(ipresent.gt.illa) then
        ipresent=1
      endif
      ii=ipresent+ipoa-1
      call dancd(i1(ii),i2(ii),j)
      value=cc(ii)
      return

      end subroutine
      subroutine daorder(ina,iunit,jx,invo,nchop)
      implicit none
      integer i,ii,ii1,ii2,iin,illa,ilma,ina,inoa,inva,invo,io,io1,     &
     &ipoa,iunit,iwarin,iwarno,iwarnv,j,jh,jt,jx,nchop
      integer(8) ic
      double precision c
!     ***************************
!
!     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
      character daname(lda)*10
      common / daname / daname
!-----------------------------------------------------------------------------3
!
      character c10*10
      dimension j(lnv),jx(lnv),jt(lnv)
!
      if(ina.lt.1.or.ina.gt.nda) then
        write(6,*)'ERROR IN DAPRI, INA = ',ina
        stop
      endif
!
      inoa = idano(ina)
      inva = idanv(ina)
      ipoa = idapo(ina)
      ilma = idalm(ina)
      illa = idall(ina)
!
      do i=1,lnv
        jt(i)=0
        j(i) = 0
      enddo
!
      call daclr(1)
!
      ic = 0
!
      iwarno = 0
      iwarnv = 0
      iwarin = 0
!
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
      read(iunit,'(A10)') c10
!
      iin = 0
!
  10  continue
      iin = iin + 1
      read(iunit,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')                     &
     &ii,c,io,(jt(i),i=1,invo)
!
      if(ii.eq.0) goto 20
!etienne
      read(iunit,*) c
      do jh=1,invo
        j(jh)=jt(jx(jh))
      enddo
      do jh=nchop+1,inva
        j(jh)=0
      enddo
!etienne
      io1 = 0
      do i=1,inva
        io1 = io1 + j(i)
      enddo
!
      ic = ic + 1
      call dadcd(j,ii1,ii2)
      ic = ia1(ii1) + ia2(ii2)
      cc(ic) = c
      goto 10
!
  20  continue
!
      call dapac(ina)
!
      return
      end subroutine
!
!ETIENNE
      subroutine datrash(idif,ina,inc)
      implicit none
      integer ibase,ider1,ider1s,ider2,ider2s,idif,ikil1,ikil2,         &
     &illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,invc,jj
      integer(8) i, ic, ipoa, ipoc
      double precision xdivi
!     ******************************
!
!     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
!     OF THE VECTOR A AND STORES THE RESULT IN C.
!
!-----------------------------------------------------------------------------1
      include "TPSALib_prm.f"
!
      integer jx(lnv)

!      call daclr(1)
!      call dacop(ina,1)
!      call dapac(ina)

      call dainf(ina,inoa,inva,ipoa,ilma,illa)
      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
!
!
!      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
!
      ibase = nomax + 1
!
      if(idif.gt.(nvmax+1)/2) then
        ider1  = 0
        ider1s = 0
        ider2  = idif-(nvmax+1)/2
        ider2s = 1
        do jj=1,ider2-1
          ider2s = ider2s*ibase
        enddo
        xdivi  = ider2s*ibase
      else
        ider1  = idif
        ider1s = 1
        do jj=1,ider1-1
          ider1s = ider1s*ibase
        enddo
        ider2  = 0
        ider2s = 0
        xdivi  = ider1s*ibase
      endif
!
      ibase = nomax+1
!
      ic = ipoc-1
!
      do i=ipoa,ipoa+illa-1
!
        call dancd(i1(i),i2(i),jx)

        ikil1=0
        ikil2=0
        if(idif.gt.(nvmax+1)/2) then
          ikil2=jx(idif)
        else
          ikil1=jx(idif)
        endif
!
!      X = IEE/XDIVI
!etienne      IFAC = INT(IBASE*(X-INT(X)+1.D-8))
!
!etienne      IF(IFAC.EQ.0) GOTO 100
!
        ic = ic + 1
        cc(ic) = cc(i)
        i1(ic) = i1(i) - ikil1*ider1s
        i2(ic) = i2(i) - ikil2*ider2s
!
! 100    continue
      enddo
!
      idall(inc) = ic - ipoc + 1
      if(idall(inc).gt.idalm(inc)) then
        write(6,*)'ERROR IN DATRASH '
        call dadeb(111,'ERR DATRAS',1)
      endif
!
      return
      end subroutine

      integer function mypause(i)
      implicit none
! Replaces obsolescent feature pause
      integer i
!      write (*,'(A,i6)',ADVANCE='NO') ' PAUSE:',i
      write (*,'(A,i6)') ' PAUSE: ',i
      read(*,*)
      mypause=i
      end function mypause

      end module tpsa
