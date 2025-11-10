      subroutine gfwigg(nstep, eps, len, lambda, rhoi, kx, x)

c     Symplectic integrator for wiggler based on generating function
c     exact in x and y but expanded to second order in h, px and py.
c     The implicit transformation is inverted by a Newton search.
c     Initially by G. Wuestefeld and J. Bahrdt at BessyII.

      implicit none

      integer	miter
      real*8	eps
      parameter	(miter = 20)

      integer	nstep, i, j
      real*8 	len, lambda, rhoi, kx, x(*)
      real*8 	rho, rhoi2, k, ky, kx2, ky2, k2, ky4
      real*8 	dx, dy, dpxi, dpyi, pxipxf, pyipyf, pi, h
      real*8	qxi, qyi, pxi, pyi, cx, sx, shy, chy, cx2, sx2, shy2, chy2
      real*8	pxf, pyf, qx, qy, z, cz, sz, cz2, sz2, k4
      real*8	ans1, ans2, ans3, ans4

      pi = 4d0*datan(1d0)

      if (rhoi .ne. 0d0) then
        rho = 1d0/rhoi
      else
	rho = 0d0
      endif
      rhoi2 = rhoi**2

      k = 2d0*pi/lambda
      ky = dsqrt(k**2+kx**2)

      kx2 = kx**2
      ky2 = ky**2
      ky4 = ky2**2
      k2 = k**2
      k4 = k2**2

      qxi = x(1)
      pxi = x(2)
      qyi = x(3)
      pyi = x(4)

      z = len/nstep

      if (dabs(nint(z/lambda)-z/lambda) .gt. 1d-5) write(*, *)
     +  '** gfwigg: step size has to be multiple of period length'

      do 10 i=1, nstep
        cx = dcos(kx*qxi)
        sx = dsin(kx*qxi)
        shy = dsinh(ky*qyi)
        chy = dcosh(ky*qyi)
        cx2=cx**2
        sx2=sx**2
        shy2=shy**2
        chy2=chy**2

c	Get starting values

      PXF=-1./2.*((KX2*SHY2-KY2*CHY2)*Z*KX*CX*SX-2.*PXI*
     . RHO**2*KY2*K2)*RHO**(-2)*KY2**(-1)*K2**(-1)

      PYF=-1./2.*((KX2-KY2)*Z*SHY*CHY*SX2+Z*KY2*SHY*CHY-
     . 2.*KY*PYI*RHO**2*K2)*KY**(-1)*RHO**(-2)*K2**(-1)

c	Newton search
	do 20 j=1, miter

      DPXI=PXI - (PXF*((KX2*SHY2-KY2*CHY2)*(CX2-SX2)*Z**2*KX2*
     . RHOI2+4.*KY2*K2)+2.*((KX2*SHY2-KY2*CHY2)*CX*RHOI2-
     . 2.*(KX2+KY2)*KY*PXF*PYF*SHY*RHOI+(KX2-KY2)*Z*KY*PYF*
     . CX*SHY*CHY*RHOI2)*Z*KX*SX-4.*(PXF**2-PYF**2)*Z*KX2*
     . KY2*CX*CHY*RHOI)/(4.*KY2*K2)

      PXIPXF=(-4.*PYF*(KX2+KY2)*Z*KX*KY*SX*SHY*RHOI+((KX2*
     . SHY2-KY2*CHY2)*Z*CX2*RHOI2-8.*PXF*KY2*CX*CHY*RHOI)*
     . Z*KX2-(KX2*SHY2-KY2*CHY2)*Z**2*KX2*SX2*RHOI2+4.*KY2
     . *K2)/(4.*KY2*K2)

      DPYI=PYI - (-(2.*(2.*(PXF**2-PYF**2)*KY2*RHOI-(KX2-KY2)*Z*
     . PXF*CX*CHY*RHOI2)*Z*KX*SX*SHY-2.*(2.*(KX2+KY2)*KY*
     . PXF*PYF*CX*RHOI+KY2*SHY*RHOI2)*Z*CHY-(2.*(KX2-KY2)*
     . SHY*CHY+(SHY2+CHY2)*Z*KY*PYF*KX2)*Z*SX2*RHOI2-(SHY2
     . +CHY2)*Z**2*KY*PYF*KY2*CX2*RHOI2-4.*KY*PYF*K2))/(4.
     . *KY*K2)

      PYIPYF=(4.*PXF*(KX2+KY2)*Z*CX*CHY*RHOI+((SHY2+CHY2)*
     . Z*KY2*CX2*RHOI2+8.*KX*KY*PYF*SX*SHY*RHOI)*Z+(SHY2+
     . CHY2)*Z**2*KX2*SX2*RHOI2+4.*K2)/(4.*K2)

	  dx = dpxi/pxipxf
	  dy = dpyi/pyipyf
          pxf = pxf + dx
          pyf = pyf + dy

          if ((dabs(dy)+dabs(dx)) .lt. eps) goto 30
20	continue

30      continue

      QX=(4.*PYF*(KX2+KY2)*Z*KY*CX*SHY*RHOI+((KX2*SHY2-
     . KY2*CHY2)*Z*CX*RHOI2-8.*PXF*KY2*CHY*RHOI)*Z*KX*SX+
     . 4.*(Z*PXF+QXI)*KY2*K2)/(4.*KY2*K2)

      QY=(4.*PXF*(KX2+KY2)*Z*CX*SHY*RHOI+8.*PYF*Z*KX*KY*
     . SX*CHY*RHOI+4.*(Z*PYF+QYI)*KY*K2+(KX2*SX2+KY2*CX2)*
     . Z**2*SHY*CHY*RHOI2)/(4.*KY*K2)

        qxi = qx
        pxi = pxf
        qyi = qy
        pyi = pyf
10    continue

      x(1) = qx  
      x(2) = pxf
      x(3) = qy
      x(4) = pyf

      return
      end
