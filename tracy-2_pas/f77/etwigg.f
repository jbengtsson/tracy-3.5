      subroutine vecpot(z, BoBrho, kx, ky, kz, x, AxoBrho, AyoBrho,
     +			radia)

      implicit none

      integer	nharm
      parameter	(nharm = 1)

      integer	radia, i
      real*8	z, BoBrho, kx, ky, kz, x(*), AxoBrho(0:3), AyoBrho(0:3)
      real*8	cx, sx, chy, shy, sz, cz

      cx = dcos(kx*x(1))
      sx = dsin(kx*x(1))

      do 10 i=0, 3
	AxoBrho(i) = 0d0
	AyoBrho(i) = 0d0
10    continue
      do 20 i=1, nharm
c	the sum over harmonics assumes kx zero
        chy = dcosh(i*ky*x(3))
        shy = dsinh(i*ky*x(3))
        sz = dsin(i*kz*z)

        AxoBrho(0) = AxoBrho(0) + BoBrho/kz*cx*chy*sz
        AyoBrho(0) = AyoBrho(0) + BoBrho*kx/(ky*kz)*sx*shy*sz

c       derivatives with respect to x
        AxoBrho(1) = AxoBrho(1) - BoBrho*kx/kz*sx*chy*sz
        AyoBrho(1) = AyoBrho(1) + BoBrho*kx**2/(ky*kz)*cx*shy*sz

c       derivatives with respect to y
        AxoBrho(2) = AxoBrho(2) + BoBrho*ky/kz*cx*shy*sz
        AyoBrho(2) = AyoBrho(2) + BoBrho*kx/kz*sx*chy*sz

        if (radia .eq. 1) then
          cz = dcos(kz*z)
c	  derivatives with respect to z
          AxoBrho(3) = AxoBrho(3) + BoBrho*cx*chy*cz
          AyoBrho(3) = AyoBrho(3) + BoBrho*kx/ky*sx*shy*cz
        endif
20    continue

      return
      end


      subroutine etwigg(nstep, len, lambda, BoBrho, kx, x,
     +			crad, pthlen, radia)

c     first order symplectic integrator for wiggler using expanded Hamiltonian

      implicit none

      integer	nstep, pthlen, radia, i
      real*8	len, lambda, BoBrho, kx, x(*), crad
      real*8	pi, ky, kz, AxoBrho(0:3), AyoBrho(0:3), B2, x2
      real*8	dp, h, hodp, det, a11, a12, a21, a22, xf(6)
      real*8	c11, c12, c21, c22, z, d1, d2, xp, yp, B(3)
      real*8	B2perp

      pi = 4d0*datan(1d0)

      kz = 2d0*pi/lambda
      ky = dsqrt(kz**2+kx**2)

      h = len/nstep

      z = 0d0
      do 10 i=1, nstep
        call vecpot(z, BoBrho, kx, ky, kz, x, AxoBrho, AyoBrho, radia)

        dp = 1d0 + x(5)
        hodp = h/dp

        a11 = hodp*AxoBrho(1)
        a12 = hodp*AyoBrho(1)
        a21 = hodp*AxoBrho(2)
        a22 = hodp*AyoBrho(2)
        det = 1d0 - a11 - a22 + a11*a22 - a12*a21

        d1 = hodp*AxoBrho(0)*AxoBrho(1)
        d2 = hodp*AxoBrho(0)*AxoBrho(2)

        c11 = (1d0-a22)/det
        c12 = a12/det
        c21 = a21/det
        c22 = (1d0-a11)/det

        x2 = c11*(x(2)-d1) + c12*(x(4)-d2)
        x(4) = c21*(x(2)-d1) + c22*(x(4)-d2)
        x(2) = x2
        x(1) = x(1) + hodp*(x(2)-AxoBrho(0))
        x(3) = x(3) + hodp*x(4)

        x(6) = x(6) + h*(((x(2)-AxoBrho(0))/dp)**2
     +	     +((x(4)-AyoBrho(0))/dp)**2)/2d0
        if (pthlen .eq. 1) x(6) = x(6) + h

        if (radia .eq. 1) then
          xp = x(2)/dp
	  yp = x(4)/dp
	  B(1) = -AyoBrho(3)
	  B(2) =  AxoBrho(3)
	  B(3) = AyoBrho(1)-AxoBrho(2)
	  B2 = B2perp(0d0, B, x, xp, yp)

          xf(5) = - crad*dp**2*B2*(1d0+(xp**2+yp**2)/2d0)

c         good in large machine: conservation of the dx/dz and dy/dz
          xf(2) = xp*xf(5)
          xf(4) = yp*xf(5)

          x(2) = x(2) + h*xf(2)
          x(4) = x(4) + h*xf(4)
          x(5) = x(5) + h*xf(5)
        endif
        z = z + h
10    continue

      return
      end
