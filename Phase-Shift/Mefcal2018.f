      program mefcal2018
c
c     	21.08.1972   merle mit fehler rms-radius                                                                               
c	mefitb? from mainz may 1982
c	modifications for dec10, jvdl
c       included input via rhogr, july 2016, hpb
c       modified constants and use AMASS instead of IA, jan. 2018, hpb
c
      implicit double precision (a-h,o-z)
      double precision lr1,lr2,fcoul
      dimension anu(25),rho(0:150)

      common/wiwq/ e,wimi,wima,dwc,ifalt,dphi,dthe,wqka(3),
     1             tardik,strl,witar,deng,ewidth,disp
      common/index/ npw1,nq,nc,lm
      common/hilf/ gam,amass,pi,wf,hc,hund,xk,gammak,betak,fmot
      common/lcom/par(32),ad(7),ipab(32),rd(13),b(35),bi(35),cm(35),
     1            cmi(35),rms,drms(26)
      common phd(35),test(35),wqab(3),va(4),dum(720,2)
      integer unit1, unit2
      parameter(unit1=1, unit2=2)

      data alfa,fmo,eps/7.297352d-3,.719987813d0,1.d-11/
      data drk,en,rmax,xfma/4*0.d0/ 
      data idn,iz,j,kb,mtwq,nd,nk,npot,nueb/9*0/
      data one/1.d0/
c     data pi/3.1415926535/

      call iofile(unit1,unit2)

      write(unit2,400)
      write(6,400)
  400 format(' *****  Mefcal version Nikhef 2018 *****',/)
c
c     einlesen und vorbesetzung der parameter
c
	drk=0.07
	nc=25
      read(unit1,*) iz,amass
      zalf=iz/137.0
c	radiation length computed 
c	830406 jvdl
      goto (31,32,33,34) iz
      lr1=alog(184.15*iz**(-1.0/3.0))
      lr2=alog(1194.0*iz**(-2.0/3.0))
      goto 35
 31   lr1=5.31
      lr2=6.144
      goto 35
 32   lr1=4.79
      lr2=5.621
      goto 35
 33   lr1=4.74
      lr2=5.805
      goto 35   
 34   lr1=4.71
      lr2=5.924
 35   zasqr=zalf**2
      fcoul=zasqr*(1.202-zasqr*(1.0369-zasqr*1.008/(1+zasqr)))
      strl=716.405*amass/((lr1-fcoul)*iz**2+lr2*iz)
c
c     coulomb-koeffizienten (lenz)
c
      gam=iz*alfa
      call dirfck(gam,nc)
c
c     parameter
c

      call rhogr(iz,rmax,anu,nq,unit1,unit2)

      npw1=nq
      nq1=nq+1
      par(32)=rmax
      par(31)=rmax
      rhosum=0.0
      nrmax=rmax*10
      dqf=pi/rmax
      do 52 j=0,nrmax
      qr=j*dqf/10.0
      rho(j)=0.0
      do 53 i=1,nq
      if (j .eq. 0) then
      	rho(j)=rho(j)+anu(i)
      else
      	rho(j)=rho(j)+anu(i)*sin(qr*i)/(qr*i)
      end if
 53   continue
      rhosum=rhosum+rho(j)*(j/10.0)**2
 52   continue
      rhosum = rhosum * 0.4 * pi
      if (abs((rhosum-iz)/iz) .ge. 1e-4) then
        write(6,1001) rhosum
        write(unit2,1001) rhosum
      else
        write(unit2,1002) rhosum
      	rhosum = iz
      endif
      do 54 j=1,nq
 54   par(j)=anu(j)/rhosum*4*rmax**3/(pi*j**2)*(-1)**(j+1)
      par(nq+1)=1.0
      par(nq+2)=1.0
c
c     energien,winkel
c
c
      call foumom(par,nq,rmax,2,rms,drms)
      write(unit2,500) iz,amass,rmax,rms
      write(6,510) iz,amass
      write(unit2,1000)
      write(unit2,1005) (anu(j),j=1,nq)
      write(unit2,1006)
      write(unit2,1007) (rho(j),j=0,nrmax)
      ifalt=0
      tardik=0.
      witar=0.
      deng=0.
      ewidth=0.
      disp=0.
      dthe=0.
      dphi=0.
   60 read(unit1,*) e, ewidth
      if (e.le.0) goto 100
      write(unit2,1020) e
      if (ewidth.lt.0) goto 70
      read(unit1,*) deng, disp
      read(unit1,*) tardik, dthe, dphi, witar
      ifalt=1
      write(unit2,1040) ewidth,deng,disp,dthe,dphi,tardik,strl
      if(witar .eq. -2.) write(unit2,1044) tardik
      if(witar .eq. -1.) write(unit2,1045)
      if(witar .eq. 0.) write(unit2,1046)
      if(witar .gt. 0.) write(unit2,1047) witar
c
c     input parameters used for correction for solid angle
c     and energy-width contributions
c
c       tardik = target thickness (g/cm2)
c       witar  = target angle (deg)
c              = -2 effective thickness, no correction
c              = -1 target in reflection
c              = 0 target in transmission
c              > 0 target angle fixed to witar
c       ewidth = beam energy width (mev)
c       deng   = energy acceptation on target (mev)
c       disp   = energy dispersion on target (mrad/mev)
c       dthe   = horizontal acceptation of spectrometer (deg)
c       dphi   = vertical acceptation of spectrometer (deg)
c
   70 read(unit1,*) wimi,wima,dwc
      tardik=tardik/1000.0
      deng=deng*e/hund
      ewidth=ewidth*e/hund
      disp=disp*hund/e/111.5
      dthe=dthe*0.18/pi
      dphi=dphi*0.18/pi
      ifdr=0
c
c     maximalradius,schrittanzahl
c     potential und rms-radius
c
      npot=rmax*e/(200.*drk)
      npot=(npot/36+1)*36
      if(npot.gt.720) npot=720
      call potenf(npot,rmax,gam,npw1)
c
c     integration der diracgleichung
c
      fmot=fmo*float(iz)/e
      betak=one/(one+amass*931.47/e)
      gammak=one/dsqrt(one-betak*betak)
      en=gammak*e*(one-betak)
      xk=en/hc
c
c     bestimmung der schrittweite
c
      nd=xk*rmax/drk
      nk=npot/nd
      if(nk.lt.1) nk=1
      if(nk.gt.3) nk=3
      nd=npot/nk
      call digraf(xk,idn,nd,nk,rmax,gam,lm,npw1,eps*hund,eps)
      call sigma
      write(6,520)e
      goto 60
  100 continue
      close (unit=unit1)
      close (unit=unit2)
c
 1000 format(' f.b.coefficients:',/)
 1001 format(' integral over rho:',f5.2,
     +  ' F.B.coefficients normalized',/)
 1002 format(' integral over rho: ',f5.2,
     + ', no normalizing performed'/)
 1005 format(6(e15.5))
 1006 format(/,' rho(r) in steps of 0.1 fm, starting at r=0: ',/)
 1007 format(10(f10.5))
 1020 format(//,' ENERGY :',f8.3,' MeV',/, 21(1h=),/)
 1040 format('  folding correction parameters: ',/,
     +  ' E-width          =',f6.2,' %',/,
     +  ' E-slit           =',f6.2,' %',/,
     +  ' E-disp           =',f6.2,' cm/%',/,
     +  ' theta            =',f6.2,' mrad',/,
     +  ' phi              =',f6.2,' mrad',/,
     +  ' target-thickness =',f6.2,' mg/cm**2',/,
     +  ' rad.length       =',f6.2)
 1044  format(' fixed effective target thickness =',f6.2,' mg/cm**2',/)
 1045  format(' target in reflection',/)
 1046  format(' target in transmission',/)
 1047  format(' fixed target-angle =',f6.2,' deg',/)
  500 format(' z:',i8,/,' a:',f8.4,/,' rmax:',f8.2,' fm',/,
     +  ' rms: ',f8.4,' fm',/)
  510 format(' rho and FB calculated for Z=',i3,' A=',f8.4)
  520 format(' calculation for energy',f7.2,' done')
      end


      SUBROUTINE RHOGR(IZ,RMAX,ANU,NNQ,UNIT1,UNIT2)
      DOUBLE PRECISION RMAX,ANU
      COMMON /GRD/ Z,RK,GNU(16),MODEL,NQ
      DIMENSION RHO(150),B(16),Q(16),ANU(16)
      integer unit1, unit2
      DATA PI /3.1415926535/
      DO 5 N=1,16
    5 GNU(N)=0.0
      Z=FLOAT(IZ)
      READ(UNIT1, *) RK,MODEL,NQ
      WRITE(UNIT2,200)
      WRITE(UNIT2,205) Z,RK,NQ
      IF(MODEL.EQ.0) READ(UNIT1,*) (GNU(I),I=1,NQ)
      IF(MODEL.EQ.0) GOTO 45
      IF(MODEL.EQ.1) READ(UNIT1,*) (GNU(I),I=1,3)
      IF(MODEL.EQ.1) WRITE(UNIT2,210) (GNU(I),I=1,3)
      IMAX=10*INT(RK)
      IF(MODEL.EQ.2) READ(UNIT1,*) (RHO(I),I=1,IMAX) 
C  calculate FB coefficients for 3pF and numerical model
      DO 10 N=1,NQ
      B(N)=0.
   10 Q(N)=FLOAT(N)*PI/RK
      DELX=RK/500.
      IF(MODEL.EQ.2) DELX=0.1
      IF(MODEL.EQ.1) IMAX=500
      DO 30 I=1,IMAx
      X=DELX*FLOAT(I)
      IF(MODEL.eq.1) CALL CHARG0(X,C0)
      IF(MODEL.eq.2) C0=RHO(I)
      C1=C0*X
      DO 20 N=1,NQ
   20 B(N)=B(N)+C1*SIN(Q(N)*X)
   30 CONTINUE
      DO 40 N=1,NQ
   40 GNU(N)=B(N)*Q(N)*2.0*DELX/RK
   45 CONTINUE
C normalize to total charge Z
      S=1.
      SUM=0.
      DO 50 N=1,NQ
      SUM=SUM+GNU(N)*S/FLOAT(N)**2
      S=-S
   50 CONTINUE
      SUM=4.0*SUM*RK**3/PI
      DO 60 N=1,NQ
      GNU(N)=GNU(N)*Z/SUM
   60 CONTINUE
   70 CONTINUE
C   calculate and print rho
      DO 80 I=1,150
  80  RHO(I)=0.0
      ANM=0.
      RMS=0.0
      MODEL=0
C    calculate rho by using the FB coefficients
      DO 85 I=1,150
      X=0.1*FLOAT(I)
      IF(X.GT.RK) GOTO 85
      CALL CHARG0(X,C0)
      A=FLOAT(1+I-2*(I/2))*X**2*0.1/1.5
      RHO(I)=C0
      ANM=ANM+C0*A
      RMS=RMS+C0*A*X**2
      IMAX=I
   85 CONTINUE
      RMS=SQRT(RMS/ANM)
      ANM=Z/(4.*PI*ANM)
      DO 90 I=1,IMAX
   90 RHO(I)=RHO(I)*ANM
C     WRITE(UNIT2,100)
C     WRITE(UNIT2,110) (RHO(I),I=1,IMAX)
C     WRITE(UNIT2,120) RMS
C     WRITE(UNIT2,125)
C     WRITE(UNIT2,130) GNU
      RMAX=RK
      NNQ=NQ
      DO 99 N=1,16
   99 ANU(N)=DBLE(GNU(N))       
C     WRITE(UNIT2,130) ANU
      RETURN

C 100 FORMAT(10X,'GROUND STATE CHARGE DISTRIBUTION:')
C 110 FORMAT(5X,10F10.5)
C 120 FORMAT(10X,'GROUND STATE RMS-RADIUS =',F8.3,' FM')
C 125 FORMAT(10X,'FB-COEFFICIENTS FOR GROUND STATE:')
C 130 FORMAT(8F12.7)
  200 FORMAT(' GROUND STATE CHARGE CALCULATION'
     2,1X,'AND FOURIER BESSEL EXPANSION')
  205 FORMAT(' NUCLEAR CHARGE:',F5.1,
     + ', CUT-OFF RADIUS:',F5.1,' FM, NUMBER OF FB COEFFICIENTS:',I3)
  210 FORMAT(' GROUND STATE IS 3-PARAMETER FERMI DISTR. WITH ',
     2  'C =',F7.4,' FM, Z =',F7.4,' FM AND W =',F7.4,/)
      END

      SUBROUTINE CHARG0(R,C0)
C---------------------------------------------------------------------------
C     Subroutine calculates ground state density C0 for models 0 and 1
C     other models not yet implemented; model 2 does not use this subroutine
C---------------------------------------------------------------------------
      COMMON /GRD/ Z,RK,GNU(16),MODEL,NQ
      EQUIVALENCE (C,GNU(1)),(A,GNU(2)),(W,GNU(3)) 
      DATA PI /3.1415926535/
      C0=0.0
      IF(MODEL.NE.0) GOTO 10
      DO 5 J=1,NQ
      B=FLOAT(J)*PI*R/RK
      C0=C0+GNU(J)*SIN(B)/B
    5 CONTINUE
      RETURN
   10 IF(MODEL.NE.1) GOTO 20
      WX=W/C**2
      C0=(1.0+WX*R*R)/(1.0+EXP((R-C)/A))
      RETURN
   20 CONTINUE
C     other models not yet implemented; model 2 does not use this subroutine
      RETURN
      END



c	function agauss
c
c	purpose
c	 evaluate integral of gaussian probability function
c
c	description of parameters
c	 x	- limit for integral
c	 integration range is +/- abs(x)
c
      function agauss(x)

      double precision x, z, y2, term, sum, denom, agauss

      z = dabs(x)
      agauss = 0.
      if (z) 42, 42, 21
21    term = 0.7071067812 * z
      sum = term
      y2 = (z**2) / 2.
      denom =1.
c
c	accumulate sum of terms
c
31    denom = denom + 2.
      term = term * (y2 + 2. / denom)
      sum = sum + term
      if (term / sum - 1.e-10) 41, 41, 31
41    agauss = 1.128379167 * sum * dexp(-y2)
42    return
      end
      block data
      implicit double precision (a-h,o-z)
      common/lcom/par(32),ad(7),ipab(32),rd(13),b(35),bi(35),cm(35),
     *cmi(35),rms,drms(26)
      common/hilf/ gam,amass,pi,wf,hc,hund,xk,gammak,betak,fmot
      data pi,hc/3.1415926536d0,197.3285d0/
      data wf,hund/1.745329252d-2,100.d0/
      data ad(1),ad(2),ad(3),ad(4),ad(5),ad(6),ad(7)/41.d0,216.d0,27.d0,
     *272.d0,27.d0,216.d0,41.d0/
      end
      subroutine cogam (xarg,yarg,nln,xgam,ygam)
      implicit double precision (a-h,o-z)
      dimension a(7)
      data pi,half,zero,xln/3.1415926536d0,0.5d0,0.d0,
     *0.9189385332d0/
      data a(1),a(2),a(3),a(4),a(5),a(6),a(7)/2.269488974d0,
     *1.517473649d0,1.011523068d0,.525606469d0,.2523809524d0,
     *.33333333333d-1,.8333333333d-1/
      x = xarg
      y = yarg
      xsum=zero
      ysum=zero
      n = 0
      m=0
    5 if(x.ge.10.) go to 10
      xsum = dlog(x*x+y*y)+xsum
      ysum = datan(y/x)+ysum
      if(x.lt.zero) m=m+1
      n = n+1
      x=xarg+n
      go to 5
   10 if(m-m/2*2.eq.1) ysum=ysum+pi
      xgam=zero
      ygam=zero
      do 20 i = 1,7
      xz=xgam+x
      yz=ygam+y
      ygam=a(i)/(xz*xz+yz*yz)
      xgam=xz*ygam
   20 ygam=-yz*ygam
      xz=half*dlog(x*x+y*y)
      yz = datan(y/x)
      if(x.lt.zero) yz=yz+pi
      xgam=(x-half)*xz-y*yz-x+xln+xgam-half*xsum
      ygam=(x-half)*yz+y*xz-y+ygam-ysum
      if(nln.eq.1) return
      ygam=dsin(ygam)/dcos(ygam)
      xgam=dexp(xgam)/dsqrt(1.+ygam*ygam)
      ygam=xgam*ygam
      return
      end
      subroutine comdiv (xzaehl,yzaehl,xnenn,ynenn,xquot,yquot)
      implicit double precision (a-h,o-z)
      an=xnenn*xnenn+ynenn*ynenn
      xquot=(xnenn*xzaehl+ynenn*yzaehl)/an
      yquot=(xnenn*yzaehl-xzaehl*ynenn)/an
      return
      end
      subroutine comult(x1,y1,x2,y2,xprod,yprod)
      implicit double precision (a-h,o-z)
      xprod = x1*x2-y1*y2
      yprod = x1*y2+x2*y1
      return
      end
      subroutine coulfm(x,xl,gam,gi,fi,gr,fr,eps,test)
      implicit double precision (a-h,o-z)
      dimension h(3)
c     double precision ad,xld,sjd,dsqrt
      data half,pih,z0,z1,z2/.5d0,1.5707963269d0,0.d0,1.d0,2.d0/
      sm=z1
      ad=gam/xl
      xld=xl
      sjd=xld*dsqrt(1.d0-ad*ad)
      xd=xld-sjd
      sj=sjd
      sjs=xl+sj
      sjd=-xd
c     diese subroutine berechnet die regulaere und irreg. coulomb-
c     funktion (loesung der dirac-gleichung fuer coulomb-potential)
c     an der stelle x (sj positiv=regulaer,sj negativ=irreg.)
c     dazu wird das von buehring angegebene verfahren der fortge-
c     setzten potenzreihen entwicklung benutzt.
  100 sj2=z2*sj
      call cogam(sj,gam,1,xgam,ygam)
      call cogam(sj2+z1,z0,1,xres,yres)
c     bestimmung des anfangswertes
      if(sm) 10,10,11
   11 ro=1.5*xl/dsqrt(xl+gam*gam*half)
      go to 12
   10 roi=0.735759*sj*sj/xl
      if(roi.gt.ro) ro=roi
   12 xa=ro
      if(xa.gt.x) xa=x
      xnorm=half*dlog(xl*sjs*half)+pih*gam+xgam-xres
      xs=xa*z2
      xnorm=xnorm+sj*dlog(xs)
      xnorm=dexp(xnorm)*sm
      merk=0
      f=gam/sjs
      g=z1
c     berechnung der potenzreihen , um 0 fuer merk=0,umy fuer merk=1
  200 a=f
      b=g
      sf=z0
      sg=z0
      m=0
      a2=z0
      b2=z0
   50 m=m+1
      xm=m
      if(merk) 1,1,2
    1 fm=xa/(xm*(xm+sj2))
      a1=fm*(-gam*a+(xm+sjd)*b)
      b=fm*(-gam*b-(xm+sjs)*a)
      go to 3
    2 xm1=m-1
      fm=hm/(xm*y)
      a1=fm*(-a*(xl+xm1)+b*d+hm*b2)
      b2=b
      b=fm*(b*(xl-xm1)-a*d-hm*a2)
      a2=a
    3 a=a1
      sf=sf+a
      sg=sg+b
      f1=f+sf
      g1=g+sg
      if(abs(a).gt.eps*abs(f1).or.abs(b).gt.eps*abs(g1)) go to 50
      f=f1
      g=g1
      y=xa
      if(xa.ge.x) go to 300
      merk=1
c     bestimmung der schrittweite hm
      h(1)=0.3*y
      h(2)=y/xl
      h(3)=x-y
      hm=1.
      do 20 i=1,3
   20 if(h(i).lt.hm) hm=h(i)
      xa=y+hm
      d=y+gam
      go to 200
  300 if(sm) 4,4,5
    5 fr=f
      gr=g
      xsr=xs
      xnr=xnorm
      sj=-sj
      sjd=-sjs
      sjs=xd
      sm=-z1
      go to 100
    4 sj=abs(sj)
      fakt=(xsr/xs)**sj
      test=(gr*f-fr*g)*fakt*gam/(z2*sj)-z1
      fi=f*xnorm
      gi=g*xnorm
      fr=fr*xnr
      gr=gr*xnr
      return
       end
      double precision function etac(xl,gam)
      implicit double precision (a-h,o-z)
      data pihalb,half,one/1.5707963269d0,.5d0,1.d0/
      ad=gam/xl
      xld=xl
      rod=xld*dsqrt(one-ad*ad)
      ro=rod
      xd=xld-rod
      nln=1
      call cogam(ro,-gam,nln,x,y)
      etac=half*datan(-gam/ro)+y+pihalb*xd
      return
      end
      double precision function etad(xl,gam)
      implicit double precision (a-h,o-z)
      data pi,one/3.1415926536d0,1.d0/
      ad=gam/xl
      xld=xl
      rod=xld*dsqrt(one-ad*ad)
      ro=rod
      xd=xld-rod
      xd=xd*pi
      a=dexp(pi*gam)
      coth=(a+one/a)/(a-one/a)
      etad=datan(-coth*dsin(xd)/dcos(xd))-xd
      return
      end
c     function determ
c   
c     purpose
c       calculate the determinant of a square matrix
c   
c     comments
c   
c       this subprogram destroys the input matrix array
c       dimension statement valid for norder up to 10
c       based on routine of bevington
c   
c     description of parameters
c   
c       array  -  matrix
c       norder -  order of determinant ( order of matrix )
c
c
      function determ(array,norder)

      double precision save, array
      dimension array(10,10)

      determ = 1.
      do 50 k = 1, norder
c
c      interchange columns if diagonal element is zero
c
      if(array(k,k)) 41, 21, 41
21    do 23 j = k, norder
      if(array(k,j)) 31, 23, 31
23    continue
      determ = 0.
      goto 60
31    do 34 i = k, norder
      save = array(i,j)
      array(i,j) = array(i,k)
34    array(i,k) = save
      determ = -determ
c
c      subtract row k from lower rows to get diagonal matrix
c
41    determ = determ  * array(k,k)
      if(k - norder) 43, 50, 50
43    k1 = k + 1
      do 46 i = k1, norder
      do 46 j = k1, norder
46    array(i,j) = array(i,j) - array(i,k) * array(k,j) / array(k,k)
50    continue
60    return
      end
      subroutine digraf(xk,idn,nd,nk,xfit,gam,lm,npw1,
     *epsd,eps)
      implicit double precision (a-h,o-z)
      dimension f1(7),g1(7),ed(35),fr(35),gr(35)
c     fuer aisdruck der wellenfunktionen
c     dimension fd(3),gd(3),nz(3)
      dimension fi(35),gi(35)
      dimension a(5),b(5),ak(33),ah(6),il(6),ih(6)
      common/lcom/par(32),ad(7),ip(32),dul(180)
      common phd(35),te(35),wqab(3),va(4),v(720),fg2(720)
      data il(1),ih(1),ah(1),ak(1),ak(2),ak(3)/2,3,12.d0,
     *5.d0,8.d0,-1.d0/
      data il(2),ih(2),ah(2),ak(4),ak(5),ak(6),ak(7)/
     *5,7,24.d0,9.d0,19.d0,-5.d0,1.d0/
      data il(3),ih(3),ah(3),ak(8),ak(9),ak(10),ak(11),ak(12)/
     *9,12,720.d0,251.d0,646.d0,-264.d0,106.d0,-19.d0/
      data il(4),ih(4),ah(4),ak(13),ak(14),ak(15),ak(16),
     *ak(17),ak(18)/14,18,1440.d0,475.d0,1427.d0,-798.d0,482.d0,
     *-173.d0,27.d0/
      data il(5),ih(5),ah(5),ak(19),ak(20),ak(21),ak(22),
     *ak(23),ak(24),ak(25)/20,25,60480.d0,19087.d0,65112.d0,-46461.d0,
     *37504.d0,-20211.d0,6312.d0,-863.d0/
      data il(6),ih(6),ah(6),ak(26),ak(27),ak(28),ak(29),
     *ak(30),ak(31),ak(32),ak(33)/27,33,120960.d0,36799.d0,
     *139849.d0,-121797.d0,123133.d0,-88547.d0,41499.d0,-11351.d0,
     *1375.d0/
      data z0,z1,z2,z3,z4,z5/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0/
      data is,lphm/0,35/
      data rv/0.d0/
      lm=0
c
c     lesen der coul.-wellenfunktionen vom disk
c
      icou=15
	if(is.eq.1) go to 122
      do 123 i=1,6
      ilw=il(i)-1
      ihw=ih(i)
      do 123 lp=ilw,ihw
  123 ak(lp)=ak(lp)/ah(i)
      is=1
  122 dx=xfit/nd
      fai=dx/140.d0
      l=0
  100 l=l+1
      xl=l
c     if(iskip(n).eq.1) write(6,121)l
c 121 format(3h1l=,i3)
      l2=2*l
      xl2=l2
c     besetzung der anfangswerte
      f1(1)=(xk-va(1))/(xl2+z1)
      g1(1)=z0
      j=0
      do 150 i=1,nd
      j=j+nk
      x=i*dx
      loop=i-i/7*7+1
      vd=xk-v(j)
      xl2x=xl2/x
      if(i.ge.2) go to 44
      fl=l2+1
      vd1=xk-va(1)
      a(1)=vd1/fl
      b(2)=-vd1*a(1)/z2
      a(3)=(vd1*b(2)-va(3))/(fl+z2)
      a(4)=-va(4)/(fl+z3)
      b(4)=(-vd1*a(3)+va(3)*a(1))/z4
      a(5)=(vd1*b(4)-va(3)*b(2))/(fl+z4)
      b(5)=(vd1*a(4)+va(4)*a(1))/z5
      f=x*(a(1)+x*x*(a(3)+x*(a(4)+x*a(5))))
      g=z1+x*x*(b(2)+x*x*(b(4)+x*b(5)))
      go to 43
   44 if(i.gt.7) go to 46
      lp=i-1
      ilw=il(lp)
      ihw=ih(lp)
      f3=dx*ak(ilw-1)
   46 xf=z0
      xg=z0
      do 11 lo=ilw,ihw
      lob=ihw-lo+ilw
      lp=loop-ihw+lo-1
      if(lp.lt.1) lp=lp+7
      xf=xf+ak(lob)*f1(lp)
      xg=xg+ak(lob)*g1(lp)
   11 continue
      xf=f+dx*xf
      xg=g+dx*xg
      a(1)=z1+f3*xl2x
      a(2)=-vd*f3
c     a(3)=-a(2)
c     a(4)=1.
      det=a(1)+a(2)*a(2)
      f=(xf-xg*a(2))/det
      g=(xf*a(2)+xg*a(1))/det
   43 f1(loop)=vd*g-xl2x*f
      g1(loop)=-vd*f
c
c     j=i-1
c     ld=j-j/3*3+1
c     fd(ld)=f
c     gd(ld)=g
c     nz(ld)=i
c     if(ld.eq.3.and.iskip(n).eq.1) write(6,120)
c     *(nz(j),fd(j),gd(j),j=1,3)
c 120 format(3(i5,1x,2e16.7,1x))
c
      fn = (x/xfit)**l
      fg2(i)=fn*(f*f+g*g)*fn
  150 continue
c     write(6,250) (fg2(i),i=1,nd)
  250 format(10e13.4)
  101 xs=x*xk
      if(icou.eq.0.and.l.le.lm) go to 131
      icou=15
      call coulfm(xs,xl,gam,gi(l),fi(l),gr(l),fr(l),eps,te(l))
      ed(l)=etad(xl,gam)
  131 fg=(f*gr(l)-g*fr(l))/(fi(l)*g-gi(l)*f)
      fl=dcos(ed(l))
      phd(l)=datan(fg*dsin(ed(l))/(z1+fg*fl))
      if(abs(f)-abs(g)) 85,85,86
   85 fn=g/(gr(l)+gi(l)*fg)
      go to 87
   86 fn=f/(fr(l)+fi(l)*fg)
   87 fn=z1/(fn*fn*(fg*(fg+z2*fl)+z1))
c
c     ableitungen nach parametern
c
      fn=fn*fai
c     write(6,250) fn,fai
      do 90 k=1,npw1
      sm=0.
      lp=ip(k)
      if(abs(phd(l)).gt.epsd.and.l.lt.lphm) go to 100
      lm=l
      if(icou.eq.0) return
      rv=xfit
   90	continue
      return
       end
      subroutine dirfak(lm,nc,gam)
      implicit double precision (a-h,o-z)
c     29.11.1972  merle
      dimension phd(35),cf(35),cfi(35)
      common/lcom/dum(39),idu(32),du(13),b(35),bi(35),cm(35),cmi(35)
     *,dul(27)
      common cl(35),cli(35),duz(1447)
      equivalence (cl(1),phd(1))
c     dirfak berechnet koeffizienten fuer entwicklung nach legendre-
c     polynomen
c     b(l),bi(l) = koff. von (fc-fquer) (berchnet in dirfck nach lenz)
c     phd(l)=phasendifferenz ausged.-punktkern
c     cm(l),cmi(l)=real und imag.teil von xl*dexp(2.*i*etac(xl,gam))
c     ausgabe
c     cl(l),cli(l)= koeff. von fa-fquer
c     cf(l),cfi(l)= koeff. von abl von fa-fquer
c
      data two/2.d0/
      lm1=lm
      if(nc.gt.lm) lm=nc
      do 80 l=1,lm
      xl=l
      eta=0.
      if(l.gt.lm1) go to 81
      eta=phd(l)
   81 m=0
      b1=dsin(eta)
      b1=-two*b1*b1
      eta=two*eta
      c1=dsin(eta)
   83 ex=cm(l)*b1-cmi(l)*c1
      ey=cmi(l)*b1+cm(l)*c1
      if(m.eq.1) go to 84
      cl(l)=ex+b(l)
      cli(l)=ey+bi(l)
      b1=dcos(eta)
      m=1
      go to 83
   84 cf(l)=-two*ey
   80 cfi(l)=two*ex
      return
       end
      subroutine dirfck(gam,nc)
      implicit double precision (a-h,o-z)
c     25.07.1972   merle
      common/lcom/dum(39),idu(32),f(6),fi(6),gm,
     =b(35),bi(35),cm(35),cmi(35),dul(27)
      dimension g(3)
c     fbl berechnet die koeffizienten der diff.reihe fc-fquer
c     bei entwicklung nach legendre-polynomen,ausserdem gibt die
c     subroutine die hilfsfunktionen f(bzw.fi) zur berechnung von
c     fquer aus.
c     eingabe
c     gam=z*alpha  nmax=anzahl der zu berechnenden koeffizienten
c     ausgabe
c     b(l),bi(l)=real u. imag-teil der koeffizienten der diff.reihe
c     f(i),fi(i)=real u. imag.-tteil der hilfsf. fuer fquer(in fqwq)
c     g(i)=potenzen von gam
c     cm(l),cmi(l)=xl*dcos(2.*etac(xl,gam))),=xl*dsin(2.*etac(xl,gam))
c
      data pi,z0,z1,z2,z4,z8,z9,ha,z3,z5,z10,z12,z15/3.1415926536d0,
     *0.d0,1.d0,2.d0,4.d0,8.d0,9.d0,.5d0,3.d0,5.d0,10.d0,12.d0,15.d0/
c     berechnung der f,fi und potenzen von gam
c
      pi2=pi*pi
      pig=pi*gam
      pig2=pig*pig
      g(1)=gam
      g(2)=g(1)*g(1)
      g(3)=g(2)*g(1)
      f(1)=g(2)/z2*(z1-pig2/z4)
      f(2)=pi*g(3)*(z3-z4*g(2))
      f(3)=g(2)/z4*(-z4+g(2)*(19.d0+pi2+g(2)*(-z4+pi2*(-z2+pig2/z4/z12
     *))))
      f(4)=pi*g(3)*(-37.5d0+g(2)/z3*(337.d0+g(2)*(-116.d0+pi2*(-7.5d0+
     *z2*g(2)))))
      fi(1)=g(3)
      fi(2)=pi*g(2)*(-z1+g(2)*(z5-pig2/z3/z2))
      fi(3)=g(3)/z3*(-z8-z3+g(2)*(z10+pi2*(z1+ha-z3/z4*g(2))))
      fi(4)=pi*g(2)*(z9+g(2)*(-82.5d0+g(2)/z3*(258.d0+z5*pi2+g(2)*(
     *-z12*z2+pi2*(-5.5d0+pig2/z4/z10)))))
      l=0
c     berechnung der b(l) und bi(l)
c
      call cogam(z1,-g(1),1,x,wi1)
      call cogam(ha,-gam,1,x,wi2)
      phi1=wi1
      phi2=wi2
      ex1=z0
      ey1=z0
  100 l=l+1
      if(l.gt.35) go to 50
      if(l.gt.nc) go to 49
      xl=l-1
      xlf=z2*xl+z1
      cphi=dcos(z2*phi1)
      sphi=dsin(z2*phi1)
      xl2=xl*xl
      b1=xl2+xl+g(2)
      a1=b1-z2
      c1=z4*(xl2+g(2))-z1
      d1=c1+z4*xl-z4-z10
      call comdiv(f(1),fi(1),a1,-z3*g(1),xq,yq)
      x=z2*(a1*a1-z15*g(2)-z4*a1)
      y=-z8*g(1)*(z2*a1-z3)
      call comdiv(f(3),fi(3),x,y,xq1,yq1)
      xs=g(2)+(z1+g(2))*xq+(z2+g(2))*xq1
      ys=(z1+g(2))*yq+(z2+g(2))*yq1
      call comdiv(xs,ys,b1,-g(1),xq2,yq2)
      x=z1+xq+xq1-xq2
      y=yq+yq1-yq2
      b(l)=x*cphi-y*sphi
      bi(l)=y*cphi+x*sphi
c
c     berechnung des zweiten gliedes
c
      call comdiv(f(2),fi(2),d1,-16.d0*g(1),xq,yq)
      x=d1*d1-384.d0*g(2)-z2*z10*d1
      y=-z2*z10*g(1)*(z2*d1-z2*z8)
      call comdiv(f(4),fi(4),x,y,xq1,yq1)
      ad=z4*(xl-g(2))-z1
      xs=(ad-z4)*xq+(ad-z8)*xq1
      ys=pi*g(2)*ad+yq*(ad-z4)+yq1*(ad-z8)
      call comdiv(xs,ys,c1,-z4*g(1),xq2,yq2)
      x=xq+xq1+xq2
      y=pi*g(2)+yq+yq1+yq2
      call comdiv(x,y,xlf+z2,z2*g(1),xs,ys)
      cphi=dcos(phi2*z2)
      sphi=dsin(phi2*z2)
      b(l)=xlf*(b(l)+xs*cphi-ys*sphi)
      bi(l)=xlf*(bi(l)+ys*cphi+xs*sphi)
   49 xl=l
      eta=etac(xl,g(1))*z2
      ex=xl*dcos(eta)
      ey=xl*dsin(eta)
      cm(l)=ex
      cmi(l)=ey
      if(l.gt.nc) go to 52
      ex1=b(l)-ex1
      ey1=bi(l)-ey1
      b(l)=ex-ex1
      bi(l)=ey-ey1
      go to 51
   52 b(l)=z0
      bi(l)=z0
   51 continue
c
c     rekursion fuer gamma-funktion
c
      phi1=phi1-datan(g(1)/xl)
      phi2=phi2-datan(g(1)/(xl-ha))
      go to 100
   50 gm=gam
c
c     endgueltige berechnung der f und fi fuer dithet
c
      phi1=wi1
      phi2=wi2
      cphi=dcos(z2*phi1)
      sphi=dsin(z2*phi1)
      cphi2=dcos(z2*phi2)
      sphi2=dsin(z2*phi2)
      call comdiv(f(2),fi(2),z2-z8*g(2),z8*g(1),xq,yq)
      call comdiv(f(3),fi(3),-z4*g(2),z2*(g(1)-g(3)),xq1,yq1)
      call comult(z2-z8*g(2),z8*g(1),z9-z4*g(2),z12*g(1),
     *xs,ys)
      call comdiv(f(4),fi(4),xs,ys,xq2,yq2)
      xs=f(1)
      ys=fi(1)
      f(1)=-g(1)*sphi
      f(2)=-pi*g(2)*sphi2/z2
      f(3)=-(xs*sphi+ys*cphi)/g(1)
      f(4)=-xq*cphi2+yq*sphi2
      f(5)=xq1*cphi-yq1*sphi
      f(6)=xq2*cphi2-yq2*sphi2
      fi(1)=g(1)*cphi
      fi(2)=pi*g(2)*cphi2/z2
      fi(3)=(xs*cphi-ys*sphi)/g(1)
      fi(4)=-xq*sphi2-yq*cphi2
      fi(5)=xq1*sphi+yq1*cphi
      fi(6)=xq2*sphi2+yq2*cphi2
      return
      end
      subroutine dithet(ifalt,xk,gk,bk,lm,fa,xc,yc)
      implicit double precision (a-h,o-z)
c
c     30.8.1973  merle
c     berechnet wq und 1. sowie 2. ableitung nach winkel th   (wqab )
c     benutzt fuer punktstreu-amplitude  das von lenz angegebene ver-
c     verfahren zur erzielung schneller konvergenz. diese koeffizienten
c     werden in dirfck berechnet und in dirfak mit koeffizienten der
c     differenzreihe f-fc zusammengezogen
c     fa(i) - wq mit 1. u. 2. winkelableitung
c     fa(i) - punkt-wq mit entspr. abl.
c     fa(i) j=3,np  abl. nach parametern und energie mit wi-abl.
c
      dimension fa(3),fsr(3),fsi(3)
      common/lcom/dum(39),ipab(32),f(6),fi(6),gam,b(35),bi(35),dul(97)
      common cl(35),cli(35),fai(3),va(4),dphd(720,2)
      data wf,pi,zero,one,two,z3/1.745329252d-2,3.1415926536d0,
     *0.d0,1.d0,2.d0,3.d0/
      data lphm/35/
c
c     transformation in cms-system und ableitung des tr.-faktors
c     nach laborwinkel
c
c     xc=dcos(theta/2.)
c     yc=dsin(theta/2.)
      cx=xc*xc-yc*yc
      sx=xc*yc*two
      a1=one-bk*cx
      wi=datan(sx/(gk*(cx-bk)))
      if(wi.lt.zero) wi=wi+pi
      dwi=one/(gk*a1)
      dwi2=-bk/gk*sx/(a1*a1)
      xk2=one/(xk*xk*two*two)
c
c     reihe fuer lenz-hilfsfunktion fquer und ableitungen
c     nach cms-winkel
c
      if=1
      if(ifalt.ne.0) if=3
      x=dcos(wi/two)
      y=dsin(wi/two)
      x2=x*x
      y2=y*y
      sx=x*y*two
      cx=dcos(wi)
      cotx=x/y
      tx=y/x
      do 10 i=1,if
      fsr(i)=zero
      fsi(i)=zero
      fa(i)=zero
      fai(i)=zero
   10 continue
      ey=one/y2
      do 11 i=1,6
      ex=ey
      do 12 j=1,if
      fsr(j)=fsr(j)+ex*f(i)
      fsi(j)=fsi(j)+ex*fi(i)
   12 ex=(i-j-2)*ex
   11 ey=ey*y
c
c     multiplikation mit dexp(2*i*gam*dlog(dsin(wi/2.))) mit
c     ableitungen
c
      gl=gam*dlog(y)*two
      gc=gam*cotx
      a1=cotx/two
      xl=1.
      ex=dcos(gl)
      ey=dsin(gl)
      do 13 i=1,if
      xm=fsr(i)
      fsr(i)=(ex*fsr(i)-ey*fsi(i))*xl
      fsi(i)=(ex*fsi(i)+ey*xm)*xl
   13 xl=xl*a1
      if(if.eq.1) go to 14
      ey=gam/(two*y2)
      fsr(3)=-gc*gc*fsr(1)+ey*fsi(1)-two*gc*fsi(2)+fsr(3)
      fsi(3)=-gc*gc*fsi(1)-ey*fsr(1)+two*gc*fsr(2)+fsi(3)
c     fehlender term von summen-ableitung
      fsr(3)=fsr(3)-fsr(2)*tx/two
      fsi(3)=fsi(3)-fsi(2)*tx/two
      fsr(2)=-gc*fsi(1)+fsr(2)
      fsi(2)=gc*fsr(1)+fsi(2)
c
c     multiplikation mit dcos(wi/2.)**2  mit abl.
c
   14 do 15 j=1,if
      i=if+1-j
      fsr(i)=fsr(i)*x2
      fsi(i)=fsi(i)*x2
      if(i.ne.3) go to 16
      fsr(i)=fsr(i)-(two*sx*fsr(2)+cx*fsr(1))/two
      fsi(i)=fsi(i)-(two*sx*fsi(2)+cx*fsi(1))/two
   16 if(i.ne.2) go to 15
      fsr(i)=fsr(i)-sx*fsr(1)/two
      fsi(i)=fsi(i)-sx*fsi(1)/two
   15 continue
c
c     legendre-polynom-reihen mit ableitungen
c
      pl1=zero
      pl2=one
      do 41 l=1,lm
      pl3=((2*l-1)*cx*pl2-(l-1)*pl1)/l
      pl1=pl2
      pl2=pl3
      k=l
c
c     cf(l) durch cl(k) aufgerufen
c
      xl=1.
      xm=pl1
      do 43 i=1,if
      a1=xl*(pl2+xm)
      fa(i)=fa(i)+cl(k)*a1
      fai(i)=fai(i)+cli(k)*a1
      xm=-xm
   43 xl=xl*l
   41 continue
      if(if.eq.1) go to 45
      fa(2)=cotx*fa(2)
      fai(2)=cotx*fai(2)
      fa(3)=-fa(3)-fa(2)/sx
      fai(3)=-fai(3)-fai(2)/sx
c
c     zusammenfassung von fquer mit der legendre-pol.-reihe
c
   45 do 46 i=1,if
      fa(i)=fa(i)+fsr(i)
   46 fai(i)=fai(i)+fsi(i)
c
c     winkel-vorfaktor mit ableitungen
c
      gc=gk/x
      fsr(1)=(one+bk*cx)*gc
      if(if.eq.1) go to 47
      fsr(2)=-bk*sx*gc
      fsr(3)=-bk*cx*gc
      fsr(3)=two*fsr(2)*fsr(2)+fsr(1)*two*(two*fsr(2)*tx+fsr(3)
     *+fsr(1)*(one+z3*tx*tx)/two/two)
      fsr(2)=fsr(1)*(two*fsr(2)+fsr(1)*tx)
   47 fsr(1)=fsr(1)*fsr(1)
c
c     betrags-quadrat der streu-amplitude sowie ableitungen nach
c     parametern mit jeweiligen ableitungen nach cms-winkel
c
      jb=1
      do 51 i1=1,if
      i=if+1-i1
      xm=zero
      do 52 k=1,i
      m=i+1-k
      ex=fa(k)*fa(m)+fai(k)*fai(m)
c     abpruefung auf k=m , i=3
      if(k.eq.m.and.i.eq.3) ex=ex*two
   52 xm=xm+ex
c     faktor 1./(4.*xk*xk)
   51 fa(i)=xm*xk2
c
c     multiplikation mit winkel-vorfaktor und ableitungen
c     ruecktransformation in laborsystem
c
      do 53 i1=1,if
      i=if+1-i1
      xm=zero
      do 54 k=1,i
      m=i+1-k
      ex=fsr(k)*fa(m)
      if(k.eq.m.and.i.eq.3) ex=ex*two
   54 xm=xm+ex
   53 fa(i)=xm
c
c     ableitungen nach laborwinkel
c
      if(if.eq.1) go to 50
      fa(3)=fa(3)*dwi*dwi+fa(2)*dwi2
      fa(2)=fa(2)*dwi
   50 continue
      return
      end
c
c	function factor
c
c	purpose
c	 calculates factorial function for integers
c
c	description of parameters
c	 n  - integer argument
c
      function factor (n)

      double precision fi, sum

      factor = 1.
      if (n-1) 40, 40, 13
13    if (n-10) 21, 21, 31
c
c	n less than 11
c
21    do  23 i = 2, n
      fi = i
23    factor = factor * fi
      goto 40
c
c	n greater than 10
c
31    sum = 0.
      do 34 i = 11, n
      fi = i
34    sum = sum + dlog(fi)
      factor = 3628800. * dexp(sum)
40    return
      end
      subroutine foumom(any,nq,rm,k,xmo,dxmo)
      implicit double precision (a-h,o-z)
      dimension rmo(2),any(20),dxmo(20)
      common/lcom/p(39),ip(32),dul(180)
c     20.8.1973  merle  berechnung der momente fuer
c     ladungsverteilung nach fourierentwicklung
      data pi/3.1415926536d0/
      kr=k+3-(k+2)/2*2
      fny=-1.
      xmo=0.
      do 11 ny=1,nq
      fny=-fny
      piny=pi*ny
      pin2=piny*piny
      rmo(1)=fny*piny*si(piny)
      rmo(2)=1.+fny
      if(k.lt.0) go to 12
      kd=k+1
      do 13 kl=kr,kd,2
      if(kl.eq.1) rmo(kr)=-rmo(kr)/pin2
      if(kl.eq.3) rmo(kr)=1.
      if(kl.gt.1) rmo(kr)=1.-kl*(kl-1)*rmo(kr)/pin2
   13 continue
   12 xmo=xmo+rmo(kr)*any(ny)
      dxmo(ny)=rmo(kr)
   11 continue
      if(k.ne.0) xmo=xmo**(1./float(k))
      if(k.eq.0) xmo=dexp(xmo)
      xmo=rm*xmo
      fk=xmo*(rm/xmo)**k
      if(k.ne.0) fk=fk/abs(float(k))
      do 14 ny=1,nq
   14 dxmo(ny)=fk*dxmo(ny)
      return
      end
      double precision function si(x)
      implicit double precision (a-h,o-z)
c     diese subroutine berechnet den integraldsinus fuer
c     werte x groesser 1 (abramowitz 5.2.8 seite 232 und
c     5.2.38 bzw. 5.2.39 seite 233 )
      dimension f(2)
      dimension a(4,2),b(4,2)
      data a(1,1),a(2,1),a(3,1),a(4,1),a(1,2),a(2,2),a(3,2),a(4,2)
     *,b(1,1),b(2,1),b(3,1),b(4,1),b(1,2),b(2,2),b(3,2),b(4,2)/
     *38.027264d0,265.187033d0,335.677320d0,38.102495d0,42.242855d0,
     *302.75786d0,352.018498d0,21.821899d0,40.021433d0,322.624911d0,
     *570.236280d0,157.105423d0,48.196927d0,482.485984d0,
     *449.690326d0,1114.978885d0/
      y=x*x
      rx=1.
      do 2 j=1,2
      oben=1.
      unten=1.
      rx=rx*x
      do 1 i=1,4
      oben=oben*y+a(i,j)
    1 unten=unten*y+b(i,j)
    2 f(j)=oben/(unten*rx)
      si=1.570796d0-f(1)*dcos(x)-f(2)*dsin(x)
      return
       end
c
c	function gamma
c
c	purpose
c	 calculate the gamma function for integers and half-integers
c
c	description of parameters
c	 x  - integer or half-integer
c
c	subroutines and function subprograms required
c	 factor(n)
c	  calculates n factorial for integers
c
      function gamma (x)

      double precision prod, sum, fi
c
c	integrize argument
c
      n = x - .25
      xn = n
      if (x - xn - .75) 31, 31, 21
c
c	argument is integer
c
21    gamma = factor(n)
      goto 60
c
c	argument is half-integer
c
31    prod = 1.77245385
      if (n) 44, 44, 33
33    if (n-10) 41, 41, 51
41    do 43 i = 1, n
      fi = i
43    prod = prod * (fi - .5)
44    gamma =prod
      goto 60
51    sum = 0.
      do 54 i = 11, n
      fi = i
54    sum = sum + dlog(fi - .5)
      gamma = prod * 639383.8623 * dexp(sum)
60    return
      end
	subroutine iofile(kan1,kan2)

	character*1 ifil(80), ofil(84)
	character*80 ifilnam
	character*84 ofilnam
	equivalence (ifil,ifilnam), (ofil,ofilnam)

	if (iargc() .eq. 1) then
	  call getarg(1, ifilnam)
	else
          write(6,1100) kan1
	  read(5,*) ifilnam
	end if
	do 10 i = 1, 80
	iend = i - 1
	if (ifil(i) .eq. ' ') goto 20
	ofil(i) = ifil(i)
10	continue
20	ofil(iend+1) = '.'
	ofil(iend+2) = 'o'
	ofil(iend+3) = 'u'
	ofil(iend+4) = 't'
	do 30 i = iend+5, 84
	ofil(i) =' '
30	continue

	open(unit=kan1,file=ifilnam,access='sequential',status='old')
	open(unit=kan2,file=ofilnam,access='sequential')
	rewind(unit=kan1)
	rewind(unit=kan2)
        return
1100	format(1x,'kanaal=',i2,'  inputfile : ',$)

	end
c
c	function pchisq
c
c	purpose
c	 evaluate integral of x**2 * exp(-x**2/2)
c	 between chisqr and infinity
c
c	description of parameters
c	 chisqr - lower limit of integral
c
c	subroutines and function subprograms required
c	 gamma(x)
c	  calculates gamma function for integers and half - integers
c
      function pchisq (chisqr,nfree)

      double precision chisqr, z , term, sum, pchisq

      if (nfree) 12, 12, 14
12    pchisq = 0.
      goto 60
14    free = nfree
      z = chisqr * free / 2.
      neven = 2 * (nfree / 2)
      if (nfree - neven) 21, 21, 41
c
c	number of degrees of freedom is even
c
21    imax = nfree / 2
      term = 1.
      sum = 0.
      do 34 i = 1, imax
      fi = i
      sum = sum + term
34    term = term * z / fi
      pchisq = sum * dexp(-z)
      go to 60
c
c	number of degrees of freedom is odd
c
41    if (z - 25) 44, 44, 42
42    z = chisqr * (free - 1.) / 2.
      goto 21
44    pwr = free / 2.
      term = 1.
      sum = term / pwr
      do 56 i = 1, 1000
      fi = i
      term = -term * z / fi
      sum = sum + term / (pwr + fi)
      if ( dabs(term/sum) - .00001) 57, 57, 56
56    continue
57    pchisq = 1. - (z**pwr) * sum / gamma(pwr)
60    return
      end
c
c      subroutine polfit
c
c      purpose
c       make a least-squares fit to data with a polynomial curve
c          y = a(1) + a(2)*x + a(3)*x**2 + ......
c
c      description of parameters
c        x          - array of data points for independent variable
c        y          - array of data points for dependent variable
c        npts       - number of pairs of data points
c        nterms     - number of coefficients (degree of polynomial + 1)
c        a          - array of coefficients of polynomial
c
c      subroutines and function subprograms required
c        determ(array,norder)
c           evaluates the determinant of a symmetric two-dimensional
c           matrix of order norder
c
c      comments
c        dimension statement valid for nterms up to 10
c

      subroutine polfit(x, y, npts, nterms, a)
      double precision x, y, a, sumx, sumy, xterm, yterm, array
      dimension x(*), y(*), a(*), sumx(19), sumy(10), array(10,10)

c
c       accumulate weigthed sums
c
      nmax = 2 * nterms - 1
      do 13 n = 1,  nmax
   13 sumx(n) = 0.
      do 15 j = 1, nterms
   15 sumy(j) = 0.
      do 50 i = 1, npts
      xi = x(i)
      yi = y(i)
      weight = 1.
      xterm = weight
      do 44 n = 1, nmax
      sumx(n) = sumx(n) + xterm
   44 xterm = xterm * xi
      yterm = weight * yi
      do 48 n = 1, nterms
      sumy(n) = sumy(n) + yterm
   48 yterm = yterm * xi
   50 continue
c
c       construct matrices and calculate coefficients
c
      do 54 j = 1, nterms
      do 54 k = 1, nterms
      n = j + k -1
   54 array(j,k) = sumx(n)
      delta = determ(array,nterms)
      if (delta) 61,57,61
   57 do 59 j = 1, nterms
   59 a(j) = 0.
      goto 80
   61 do 70 l = 1, nterms
      do 66 j = 1, nterms
      do 65 k = 1, nterms
      n = j + k - 1
   65 array(j,k) = sumx(n)
   66 array(j,l) = sumy(j)
   70 a(l) = determ(array,nterms) / delta
   80 return
      end
      subroutine potenf(n,rmax,gam,np)
      implicit double precision (a-h,o-z)
      common/lcom/par(39),ip(32),xdu(153),rms,drms(26)
      common dum(70),dum1(3),va(4),v(720),x(720)
c
c     aufsummation der einzelanteile des potentials dvp(k,i) mit
c     faktoren par(i) zum potential v fuer digrad
c
      data pi,half,one,z6/3.1415926536d0,.5d0,1.d0,6.d0/
c
c     vorbesetzung
c
      dx=n
      dx=rmax/dx
      do 15 i=1,4
   15 va(i)=0.
      fg=-gam/rmax
      fm=-one
      do 10 k=1,n
      x(k)=k*dx
   10 v(k)=fg
      rms=0.d0
      do 11 i=1,np
      fm=-fm
      va(1)=va(1)+(fm+1.)*par(i)*fg
      qny=pi/rmax*i
      drms(i)=rmax*rmax*(one-z6/(pi*pi*i*i))
      rms=rms+drms(i)*par(i)
      va(3)=va(3)-fm*fg*qny*qny/z6*par(i)
      fq=fm/qny*fg
      j=ip(i)
      do 12 k=1,n
      dv=fq/x(k)*dsin(x(k)*qny)
   12 v(k)=v(k)+par(i)*dv
   11 continue
      rms=dsqrt(rms)
      do 13 i=1,np
c  13 drms(i)=0.5*(drms(i)/rms-rms)
   13 drms(i)=half*drms(i)/rms
      return
      end
      subroutine sigma
      implicit double precision (a-h,o-z)
      common/wiwq/ e,wimi,wima,dwc,ifalt,dphi,dthe,wqka(3),
     1             tardik,strl,witar,deng,ewidth,disp
      common/index/ npw1,nq,nc,lm
      common/hilf/ gam,amass,pi,wf,hc,hund,xk,gammak,betak,fmot
      common/lcom/ par(31),rmax,ad(7),ipab(32),dul(180)
      common dum(70),wqab(3),va(4),dum2(720,2)
      dimension thx(9),yr(9),crka(5),temp(3)
      integer unit1, unit2
      parameter(unit1=1,unit2=2)

      data z1,z2/1.d0,2.d0/
      fab=wf*wf/6.d0
c
c     koeffizienten der legendre-polynom-reihe
c
      call dirfak(lm,nc,gam)
c
c     verschiedene winkel
c
      if (dwc.ne.0.or.wima.le.wimi) go to 427
      dwc=(wima-wimi)/10
      write(6,428) wimi,wima,dwc
  427 if(ifalt .ne. 1) goto 20
      x=0.
      fact=0.
      if (ewidth .eq. 0. .and. deng .eq. 0.) x=0.
      if (ewidth .eq. 0. .and. deng .ne. 0.) x=3.
      if (ewidth .ne. 0.) x=dabs(deng/ewidth)
      if (x .le. 0.3) fact=deng**2/24.
      if (x .gt. 0.3 .and. x .lt. 3.0)
     +   fact=(1.-pchisq(x,3))/agauss(x)/8.*ewidth**2
      if (x .ge. 3.0) fact=ewidth**2/8.
   20 nwq=(wima-wimi)/dwc+1
	if (wima.le.wimi) nwq=1
      write(unit2,14)
      witarf = wf * witar
      do 413 n=1,nwq
      wi1=wimi+float(n-1)*dwc
      wif = wf * wi1
      if (ifalt .ne. 1) goto 10
c
c target thickness dependent : independent,transmission,reflection or fixed
c
      if (witar .eq. -2) then
	thick = tardik
      else
	if (witar .eq. -1) a1 = wf * (180.0 - wi1/2.0)
	if (witar .eq. 0)  a1 = wf * (180.0 - wi1)/2.0
	if (witar .gt. 0)  a1 = wf * witar
	a2 = wf * wi1 + a1
	thick = tardik * 0.5 * (dabs(1./dsin(a1)) + dabs(1./dsin(a2)))
      endif
      eps2 = (21.5/e)**2 * thick / (4.*strl*fab)
   10 swi=dsin(wif/z2)
      cwi=dcos(wif/z2)
      fwi=swi/cwi
      cotw=(z1/fwi-fwi)/z2
      am=amass
      rueck=z1/(z1+z2*e/(am*931.47)*swi*swi)
      q=z2*swi*e/hc*dsqrt(rueck)
      xkv1=e/(va(1)*hc)
      xkv=z1-z1/xkv1
      qkor=q*xkv
      call dithet(1,xk,gammak,betak,lm,wqka,cwi,swi)
      wqmot=z1/(swi*fwi)*fmot
      wqmot=wqmot*wqmot*rueck
      dwqt=wqka(2)/wqka(1)
      dwt2=wqka(3)/wqka(1)
      if (ifalt .ne. 1) goto 40
      thfalt=dthe*wf+0.5*cotw*(wf*dphi)**2
      npts=9
      nterms=5
      dth=thfalt/float(npts)
      thmi=wif-dth*(float(npts)/2.0-0.5)
      do 30 i=1,npts
      th=thmi+(i-1)*dth
      sth=dsin(th/z2)
      cth=dcos(th/z2)
      call dithet(0,xk,gammak,betak,lm,temp,cth,sth)
      thx(i)=i-0.5-float(npts)/2.0
   30 yr(i)=temp(1)/wqka(1)
      call polfit(thx,yr,npts,nterms,crka)
      facu=1.
      do 50 i=2,nterms
      crka(i)=facu*dth**(1-i)/crka(1)*crka(i)
50    facu=facu*float(i)
c
c     berechnung d. 1.und 2.energie-abl.-jeweils dividiert durch wqka(1,
c     bei dwqe faltung ueber sp.-oeffnung beruecksichtigt
c
   40 wif=z2*fwi/(e*xkv)
      ge=xkv*xkv/(fwi*fwi)
      aa=z2/(e*(xkv1-z1))
      be=-(z1/ge+z1)/fwi
      dwqe=wif*(dwqt-be)+aa
      dcre=wif*(crka(2)-be)+aa
      de=z1/(z2*ge)+z2+1.5*ge
      ca1=aa*(z2*dwqe-aa*(z1/z2+xkv1))
      ca2=aa*(z2*dcre-aa*(z1/z2+xkv1))
      dwe2 = wif*wif*(dwt2-de+dwqe*(e*xkv/z2/z2-z2*be/wif))+ca1
      dce2 = wif*wif*(crka(3)-de+dcre*(e*xkv/z2/z2-z2*be/wif))+ca2
      dwqt=dwqt*wf
      form=wqka(1)/wqmot
      dwqe=dwqe*hund
      dwqt=dwqt*hund
      wqp=wqka(1)
      if (ifalt.ne.1) goto 420
      falt1=1.0-fab*(dthe*dthe+dphi*dphi)+
     + fab/wqka(1)*(wqka(2)*(eps2+dphi*dphi)*cotw+
     + wqka(3)*(eps2+dthe*dthe))+
     + fact*(dwe2+(disp)**2*(3.+cotw*wqka(2)/wqka(1)))
      falt2=1.0+fab*(dthe*dthe+dphi*dphi)+
     + fab*(crka(2)*(eps2+dphi*dphi)*cotw+
     + crka(3)*(eps2+dthe*dthe)+
     + crka(4)/6.*cotw*fab*(dthe*dphi)**2+
     + crka(5)*3./10.*fab*(dthe)**4)+
     + fact*(dce2+(disp)**2*(3.+cotw*crka(2)))
      wqp=falt1*wqka(1)
      write(unit2,104) wi1,q,qkor,wqmot,form,dwqe,dwqt,wqp,falt1,falt2
      goto 430
  420 write(unit2,104) wi1,q,qkor,wqmot,form,dwqe,dwqt,wqp
  430 continue
  413 continue
	write(unit2,16)
      return
c
  14  format('  angle   q   q_eff  sigma_mott      sigma/     ',
     +	'dsigma/   dsigma/    sigma      falt.coeff.      testfalt',/,
     +	'                          ',
     +	'         sigma_mott      de     dtheta',/,
     +	' [deg]  [fm-1][fm-1] [fm^2/sr]                  [%/mev]',
     +	'  [%/deg]   [fm^2/sr]',/)
   16 format(7h -100.0)
   18 format(2x//)
  103 format(1x,f6.2)
  104 format(f7.3,2f6.3,2e13.5,2f9.2,e14.5,f14.5,f14.5)
  428 format(' theta(min)=',f5.2,' theta(max)=',f5.2,' no stepsize given
     +       !, set to ',f5.2,/)
      end
