 
      subroutine ruptdirct2024 ( mag, U, T, Ry0, Smax1, Smax2, Ztor, 
     1                 Rake, specT, Version, lnfd, PhiRed) 

      implicit none
      include 'pfrisk.h'
C   input            
      real mag, U, T, Ry0, Smax1, Smax2, Ztor, Rake     
      real specT
      integer Version
C   output
      real lnfd, PhiRed 
C   internal
      integer nper, iflag, count1, count2, i
      real pi, RakeRad, PhiPer(13), e1(13), Amax, k, SigG, C
      real Tpeak, x, A, S, Srake, D, S2, fs2
      real theta, ftheta, Rt, R1, AR, Rmax, fdist
      real fztor, fG, L1, L2, fGbar, fGprime, fD
      
C***************************************************************************************************************************      
C      REQUIRED FUNCTION INPUT
C      
C      
C      FUNCTION OUTPUT
C      
C      tobe written
C***************************************************************************************************************************             
       
      pi = acos(-1.0)
      RakeRad=Rake*pi/180
      data PhiPer/ 0.01  , 0.3 , 0.4  , 0.5 , 0.75  , 1.0  , 1.5  , 2.0, 
     1             3.0   , 4.0  , 5.0 , 7.5 , 10. /      
      if (Version. eq. 1 ) then
C          -- simulation based model
       data e1/ 0.00, 0.000, 0.0003, 0.011, 0.038, 0.072, 0.107,0.143, 
     1             0.172,0.189,0.195,0.206,0.200/   
       Amax=0.54
       k=1.58
       SigG=0.38
      elseif (Version. eq. 2 ) then
C          -- data based model
       data e1/ 0.00, 0.000, 0.0024, 0.0074, 0.024, 0.041, 0.064,0.076, 
     1             0.091,0.110,0.124,0.145,0.157/   
       Amax=0.34
       k=1.58
       SigG=0.26
      endif
      if (specT.lt.0.01) specT=0.01 ! lower limit
      if (specT.gt.10.) specT=10. ! upper limit

C     -- Interpolate to get PhiRed
      nper=13
      C=0.0
      do i=1, nper-1
         if (specT .ge. PhiPer(i) .and. specT .le. PhiPer(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo
      if (specT .gt. PhiPer(nper) ) then
            C=e1(nper)
      endif
      goto 1002
1001     call interp (PhiPer(count1),PhiPer(count2),e1(count1),
     1                e1(count2),specT,C,iflag)
1002  continue
            
C     -- Continue with median adjustment
C     Coeffs
      Tpeak = 10**(-2.15+0.404*mag)
      x = log10(specT/Tpeak)
      A = Amax*exp(-x**2 /(2*SigG**2) )

C     Get S and theta
      S = U
      if (S .lt. 0) then
       if (abs(S) .gt. abs(Smax1)) then
        S = abs(Smax1)
       else
        S = abs(U)
       endif
      else
       if (S .gt. Smax2) then
        S = Smax2
       endif
      endif
      Srake = S*cos(RakeRad) 
      D = 3
      S2 = sqrt(D**2 + Srake**2)
      fs2 = log(S2)
      if (fs2 .gt. 6.142) then
        fs2 = 6.142
      endif
      theta = abs(atan(T/U))
      if (T .eq. 0 .and. U .eq. 0) then
       theta = 0
      endif
      ftheta = abs(cos(2*theta))

C     Distance Taper
      Rt = sqrt(T**2 + Ry0**2 + Ztor**2)
      if (Rt .le. 0.1) then
       Rt = 0.1
      endif
      R1 = sqrt(T**2 + Ry0**2 )
      if (mag .lt. 5) then
       Rmax=40
      elseif (mag .gt. 7) then  
       Rmax=80
      else
       Rmax=-60+20*mag
      endif
      AR = -4*Rmax
      fdist = 1-exp(AR/Rt - AR/Rmax)
      if (Rt .ge. Rmax) then
       fdist = 0
      endif

C     Ztor Taper
      if (abs(Ztor) .lt. 20) then
       fztor=1-abs(Ztor)/20
      else
       fztor=0
      endif

C     Directivity predictor fG
      fG=fs2*ftheta
      L2=abs(Smax1)
      L1=Smax2
         call fg_center_sub (L1,L2,R1,D,Rake,fGbar)
      fGprime=fG-fGbar

C     Calculate Directivity adjustment fD
      fD=A*(2/(1+exp(-k*fGprime*fdist*fztor))-1)

C     Final parameters, median and variability adj
      lnfd = fD
      if (Rt .le. Rmax) then
       PhiRed = -C
      else
       PhiRed = 0
      endif
             
      return
      END

C-----------------------------------------------------------------------

      subroutine fg_center_sub (L1,L2,R1,D,Rake,fGbar) 
     
      implicit none
      
      real L1, L2, R1, D, Rake, fGbar
      real pi, dx, rk, ysum
      real x, xr, s, l, lr, r
      real y1, y2, y3, y4
      integer n1, n2, n3, n4, i1

      pi = acos(-1.0)

      if (R1 .lt. 0.1) then
       R1 = 0.1
      endif

      dx=0.1;
      rk=Rake*pi/180
      ysum=0

C between the ends of the fault, strike direction
      n1=int(L1/dx+1)
      
      do i1=1,n1
       x=(i1-1)*dx
       xr=x*cos(rk)
       s=sqrt(xr**2 + D**2)
       y1=log(s)*abs(cos(2*atan(R1/x)))
       ysum=ysum+y1
      enddo

C between the ends of the fault, anti-strike direction
      n2=int(L2/dx)
      do i1=1,n2
       x=(i1)*dx
       xr=x*cos(rk)
       s=sqrt(xr**2 + D**2)
       y2=log(s)*abs(cos(2*atan(R1/x)))
       ysum=ysum+y2
      enddo

C off the end of the fault, strike direction
      l=L1
      lr=l*cos(rk)
      s=sqrt(lr**2 + D**2)
      n3=int(R1/dx)
      do i1=1,n3
       x=l+(i1)*dx
       r=sqrt(R1**2 - (x-l)**2)
       y3=log(s)*abs(cos(2*atan(r/x)))
       ysum=ysum+y3
      enddo

C off the end, anti-strike dir
      l=L2
      lr=l*cos(rk)
      s=sqrt(lr**2 + D**2)
      n4=int(R1/dx)
      do i1=1,n4
       x=l+(i1)*dx
       r=sqrt(R1**2 - (x-l)**2)
       y4=log(s)*abs(cos(2*atan(r/x)))
       ysum=ysum+y4
      enddo

C   total
      fGbar=ysum/(n1+n2+n3+n4)
        
      return
      end     