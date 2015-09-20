!     PARAMETERS:
!
!     LDA: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
!     LST: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
!     LEA: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO, NV
!     LIA: DIMENSION OF IA1, IA2;           CAN BE INCREASED FOR LARGE NO, NV
!     LNO: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
!     LNV: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
!
!-----------------------------------------------------------------------------1
      integer  i1, i2, ia1, ia2, idall, idalm, idano, idanv, idapo,     &
     &         ie1, ie2, ieo, ifi, lfi, nda, ndamaxi, nmmax, nocut,     &
     &         nomax, nst, nvmax
      double precision  cc, eps, epsmac, facint
!      parameter(lda=30000, lst=132000000, lea=500000, lia=50000,        &
      parameter(lda=30000, lst=7000000, lea=500000, lia=50000,          &
!     RHIC
!      parameter(lda=15000, lst=130000000, lea=500000, lia=50000,        &
     &          lno=15, lnv=8)
      common /fordes/ nda, ndamaxi
      common /da/ cc(lst), eps, epsmac
      common /dai/ i1(lst), i2(lst), ie1(lea), ie2(lea), ieo(lea),      &
     &             ia1(0:lia), ia2(0:lia), ifi(lea), idano(lda),        &
     &             idanv(lda), idapo(lda), idalm(lda), idall(lda),      &
     &             nst, nomax, nvmax, nmmax, nocut, lfi
      common /factor/ facint(0:lno)
!-----------------------------------------------------------------------------9
