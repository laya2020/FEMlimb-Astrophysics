MODULE lensFs2
USE nrtype
USE nrutil
USE Events
IMPLICIT NONE

	INTEGER(I1B) :: Nd1, Nd2
	INTEGER	::  Ncnst
	REAL(DP), ALLOCATABLE, DIMENSION(:)	:: Inodes, rnodes, Al1, tLC 
	REAL(DP)    ::  c1, c2		
	!-------------------------------------------------------------------------!
	!## below are variables for meshing procedure:
	!##	Nt is number of elements in l(t) & Nndt is number of Nods in l(t) 	##!
	!##	Nr is number of elements in r and Nndr is number of Nods in r		##!
	!-------------------------------------------------------------------------!	
	INTEGER(I2B)    :: Nt, Nndt, Nnd, Nr,  Nndr, near
	
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: l, dell, rv, delr										  										
	INTEGER(I2B), ALLOCATABLE, DIMENSION(:,:)	:: nglobr, nglobl			
	REAL(DP), PARAMETER	:: rmin=eps, eps2=1.d4*eps
	REAL(DP), PARAMETER	:: invsqrt2pi=1.d0/dsqrt(2.d0*PI_d), accuracy=1.d-7, &
							h1=1.d-4, hmin=1.d1*eps, invPI=1/PI_D	
	
CONTAINS

FUNCTION mFLUX(rmin, ll, Nn)
IMPLICIT NONE
			
		INTEGER(I2B), INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		REAL(DP), INTENT(IN)  :: rmin
		!REAL(DP), INTENT(OUT) :: FBL
		REAL(DP), DIMENSION(1:Nn) :: mFLUX
		INTEGER(I1B)	:: nf
		INTEGER(I2B)	:: il
		REAL(DP)  :: li, r, dr, rc
		REAL(DP), DIMENSION(1:1)	:: ystart
		
		!print*, 'FUNCTION mFLUX, c1, c2:',c1,c2
	 	nf=4_I1B ! => Integrand: r*I_Nsqrt(r)*A(l,r)
		dr=1.d0	!won't use in case of nf=4,16
		rc=1.d0	!won't use in case of nf=4,16
	 	r=rmin
	DO il=1, Nn
		li=ll(il)
		ystart=0
		CALL odeint(nf,li,dr,rc,  ystart, r, 1.d0)	
		mFLUX(il)=ystart(1)
	END DO	
								
END FUNCTION mFLUX
!======================!
FUNCTION FLUX_step(rmin, ll, Nn)
IMPLICIT NONE
			
		INTEGER(I2B), INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		REAL(DP), INTENT(IN)  :: rmin
		!REAL(DP), INTENT(OUT) :: FBL
		REAL(DP), DIMENSION(1:Nn) :: FLUX_step
		INTEGER(I1B)	:: nf
		INTEGER(I2B)	:: il
		REAL(DP)  :: li, r, wi, dr, rc
		REAL(DP), DIMENSION(1:1)	:: ystart
 		
		!print*, 'FUNCTION FLUX_step, c1, c2:',c1,c2
	  		!WRITE(*,*) 'mFLUX: Curve Νο.',ic
	 	nf=16_I1B ! => Integrand: r*NLD_step(r)*A(l,r)
		dr=1.d0	!won't use in case of nf=4,16
		rc=1.d0	!won't use in case of nf=4,16
	 	r=rmin
	DO il=1, Nn
		li=ll(il)
		ystart=0
		CALL odeint(nf,li,dr,rc,  ystart, r, 1.d0)!, accuracy, h1, hmin)	
		FLUX_step(il)=ystart(1)
		!print*,'FUNCTION FLUX_step:',il,Nn,li,FLUX_step(il)
	END DO	
								
END FUNCTION FLUX_step
!======================!
FUNCTION FLUX_fringe(rmin, ll, Nn)
IMPLICIT NONE
			
		INTEGER(I2B), INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		REAL(DP), INTENT(IN)  :: rmin
		!REAL(DP), INTENT(OUT) :: FBL
		REAL(DP), DIMENSION(1:Nn) :: FLUX_fringe
		INTEGER(I1B)	:: nf
		INTEGER(I2B)	:: il
		REAL(DP)  :: li, r, wi, dr, rc
		REAL(DP), DIMENSION(1:1)	:: ystart
 		
		print*, 'FUNCTION FLUX_step, c1, c2:',c1,c2
	  		!WRITE(*,*) 'mFLUX: Curve Νο.',ic
	 	nf=17_I1B ! => Integrand: r*NLD_fringe(r)*A(l,r)
		dr=1.d0	!won't use in case of nf=4,16
		rc=1.d0	!won't use in case of nf=4,16
	 	r=rmin
	DO il=1, Nn
		li=ll(il)
		ystart=0
		CALL odeint(nf,li,dr,rc,  ystart, r, 1.d0)!, accuracy, h1, hmin)	
		FLUX_fringe(il)=ystart(1)
	END DO	
								
END FUNCTION FLUX_fringe
!======================!
FUNCTION FLUX_tanh(rmin, ll, Nn)
IMPLICIT NONE
			
		INTEGER(I2B), INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		REAL(DP), INTENT(IN)  :: rmin
		!REAL(DP), INTENT(OUT) :: FBL
		REAL(DP), DIMENSION(1:Nn) :: FLUX_tanh
		INTEGER(I1B)	:: nf
		INTEGER(I2B)	:: il
		REAL(DP)  :: li, r, wi, dr, rc
		REAL(DP), DIMENSION(1:1)	:: ystart
 		
		print*, 'FUNCTION FLUX_step, c1, c2:',c1,c2
	  		!WRITE(*,*) 'mFLUX: Curve Νο.',ic
	 	nf=18_I1B ! => Integrand: r*tanhLD(r)*A(l,r)
		dr=1.d0	!won't use in case of nf=4,16
		rc=1.d0	!won't use in case of nf=4,16
	 	r=rmin
	DO il=1, Nn
		li=ll(il)
		ystart=0
		CALL odeint(nf,li,dr,rc,  ystart, r, 1.d0)!, accuracy, h1, hmin)	
		FLUX_tanh(il)=ystart(1)
	END DO	
								
END FUNCTION FLUX_tanh
!======================!

FUNCTION FLUX_FEM(rmin, ll, Nn)
IMPLICIT NONE
		!REAL(DP), INTENT(IN)	:: u1, u2	
		INTEGER(I2B), INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		REAL(DP), INTENT(IN)  :: rmin
		!REAL(DP), INTENT(OUT) :: FBL
		REAL(DP), DIMENSION(1:Nn) :: FLUX_FEM
		REAL(DP), DIMENSION(1:10)	:: x, w		
		INTEGER	  :: il, ig
		INTEGER(I1B)	:: nf
		REAL(DP)  :: li, r, wi, dr, rc, sumg!, accuracy, h1, hmin
		REAL(DP), DIMENSION(1:1)	:: ystart
 		
	 	nf=15_I1B ! => Integrand: r*I_FEM_2Nodes(r)*A(l,r)
		dr=1.d0	!won't use in case of nf=4
		rc=1.d0	!won't use in case of nf=4
	 	r=rmin
	DO il=1, Nn
		li=ll(il)
		ystart=0
		CALL odeint(nf,li,dr,rc,  ystart, r, 1.d0)!, accuracy, h1, hmin)	
		FLUX_FEM(il)=ystart(1)
	END DO	!it: time
								
END FUNCTION FLUX_FEM
	
!======================!
	FUNCTION FRAND0(Flux, Nn, sigma, nsigma, idum1)
		INTEGER(I2B), INTENT(IN)	:: Nn
		REAL(SP), INTENT(IN)	:: nsigma
		INTEGER, INTENT(INOUT) :: idum1
		INTEGER :: i ,j!, idum2
		REAL(DP), DIMENSION(1:Nn), INTENT(IN)	:: Flux, sigma
		!REAL(DP), INTENT(IN)	:: sigma
		!REAL(DP), DIMENSION(1:Nndt,1:NRan)	::	FRAND
		REAL(DP), DIMENSION(1:Nn)	::	FRAND0
		REAL(DP), DIMENSION(1:Nn)	:: maxG!, sigma2
		REAL(DP)	::	randF, randG, df, GaussF, df2d2!, sumran
		
		
	  	!WRITE(*,*) 'FRAND0: Curve Νο.',ic
	  	maxG=0.d0
		maxG=invsqrt2pi/sigma!1.d0/(sigma*sqrt2_d*dsqrt(PI_d))
		
		j=0
		i=0
			!PRINT*,'Random func. seeds:',idum1, idum2
			!write(*,*) 'lensFs2.f95, SUB. FRAND): (dF)/sigmaj, GaussF, RandG	, MaxG'				
		DO
		i=i+1
			randG=ran2(idum1)*maxG(j+1)
			randF=ran2(idum1)
			df=(2*nsigma*randF-nsigma)!*sigma(j+1)
			df2d2=df**2/2
			GaussF=maxG(j+1)*exp(-df2d2)!/sigma2(j+1))
						!write(*,'(I10,4(f12.4,x))')i, sigma2(j+1), GaussF, randG, maxG(j+1)
		IF (randG<GaussF) THEN
			j=j+1
			FRAND0(j)=df*sigma(j)+Flux(j)			

			IF (j==Nn) EXIT			
		END IF

		
		END DO										
		!WRITE(*,*) i, j	, sumran/Nn
	END FUNCTION FRAND0
!======================!
	FUNCTION FI0(rmin, ll, Nn)

		INTEGER, INTENT(IN)  :: Nn
		REAL(DP), DIMENSION(1:Nn), INTENT(IN) :: ll
		!REAL(DP), DIMENSION(1:Nt), INTENT(IN)   :: dell
		!~INTEGER(I2B), INTENT(OUT)  :: Nndt
		REAL(DP), INTENT(IN)  :: rmin
		REAL(DP), DIMENSION(1:10) :: x,w
		REAL(DP), DIMENSION(1:Nn) :: FI0
		REAL(DP), DIMENSION(1:1)	:: ystart
		INTEGER	  :: il, ig
		INTEGER(I1B)	:: nf
		REAL(DP)  :: li, r, dr,rc,sumg
		
		dr=1.d0
		rc=1.d0
		nf=5_I1B	! => Integrand: r*A(l,r)
	 	r=rmin
		WRITE(*,*) 'FI0 FUNC. curve No.',ic
		DO il=1, Nn
			li=ll(il)
			ystart=0
			CALL odeint(nf,li,dr,rc, ystart, r, 1.d0)!, accuracy, h1, hmin)	
			FI0(il)=ystart(1)
		END DO	
			
	END FUNCTION FI0	
!======================!
 FUNCTION FI0l(nf,rmin, li)
	
		REAL(DP), INTENT(IN)  :: rmin, li
		REAL(DP), DIMENSION(1:10) :: x,w
		REAL(DP) :: FI0l
		REAL(DP), DIMENSION(1:1)	:: ystart
		INTEGER	  ::  ig
		INTEGER(I1B), INTENT(IN)	:: nf
		REAL(DP)  :: r,dr,rc,sumg
		

		sumg=0
		dr=1.d0
		rc=1.d0
	 	r=rmin
		ystart=0
			CALL odeint(nf,li,dr,rc, ystart, r, 1.d0)!, accuracy, h1, hmin)	
		FI0l=ystart(1)
		
	END FUNCTION FI0l		
!======================!
	FUNCTION IntegA(li, r)
	USE Elliptic
		REAL(DP), INTENT(IN) :: r, li
		REAL(DP)     :: IntegA, mu, rml, l2pE2, l2, &
						lpr, lmr2, lpr2, en, ak, abslmr

		l2  =li**2
		rml =r-li
		lmr2=rml**2
		abslmr=dabs(rml)
		lpr =li+r
		lpr2=lpr**2
		l2pE2=l2+Event(ic)%El2

	IF (abslmr>1.d-6) THEN
		
		IF (r>1.d-15.and.li>1.d-15) THEN
			mu=4.d0/(lpr*dsqrt(lmr2+4*Event(ic)%El2))
			en=-4*li*r/lpr2
			ak=4*Event(ic)%El/lpr*dsqrt(li*r/(lmr2+4*Event(ic)%El2))
			IntegA=mu*(2*Event(ic)%El2*ellf_d(PIO2_d, ak)+lmr2*ellpi_d(PIO2_d, en, ak))
		ELSE IF (r<=1.d-15) THEN
			IntegA=4/(li*dsqrt(l2+4*Event(ic)%El2))*PIO2_d*(2*Event(ic)%El2+l2)
		ELSE 
			IntegA=4/(r*dsqrt(r**2+4*Event(ic)%El2))*PIO2_d*(2*Event(ic)%El2+r**2)		
		END IF
	ELSE
		
		IF (abslmr<=1.d-7) THEN  !because for very small abslmr  elliptic integration is expensive
			rml=dsign(1.d-7,rml)
			abslmr=1.d-7
		END IF	
		mu=Event(ic)%El*(l2+l2pE2)/(l2*l2pE2)
		IntegA=2*Event(ic)%El/li*(1-rml/(2*li))*dlog( 8*Event(ic)%El*li/(abslmr*dsqrt(l2pE2)) )+ &
				4*datan(li/Event(ic)%El)+mu*rml
!IF (abslmr<=1.d-8) write(*,*) rml, IntegA		
	END IF

		
	END FUNCTION IntegA
!======================!
	FUNCTION testMAG(rho, t)
		REAL(DP), INTENT(IN) :: rho, t
		REAL(DP) ::  testMAG

		testMAG=(t+rho)/dsqrt(t**2+rho**2)
	
	END FUNCTION testMAG	
!========================!
	FUNCTION QLD(r, I0, a2, a4)
		REAL(DP), INTENT(IN) :: r, I0, a2, a4
		REAL(DP)	:: QLD, mu, mu2
		
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			mu2=(1.d0-r)*(1.d0+r)
			QLD=I0*( 1.d0-a2*(1.d0-mu)-a4*(1.d0+mu2-2.d0*mu) )	
			
	END FUNCTION QLD
!========================!
	FUNCTION NQLD(r, a2, a4)
		REAL(DP), INTENT(IN) :: r, a2, a4
		REAL(DP)	:: NQLD, mu, mu2, nrml
		
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			mu2=(1.d0-r)*(1.d0+r)
			nrml=PI_D*(1.d0-a2/3.d0-a4/2.d0)
			NQLD=nrml*( 1.d0-a2*(1.d0-mu)-a4*(1.d0+mu2-2.d0*mu) )	
			
	END FUNCTION NQLD
!========================!
	FUNCTION NSQrtLD(r, cLD, dLD)
		REAL(DP), INTENT(IN) :: r, cLD, dLD
		REAL(DP)	:: NSQrtLD, mu, sqrtMu, nrml
		
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			sqrtMu=dsqrt(mu)
			nrml=1.d0/(PI_D*(1.d0-cLD/3.d0-dLD*0.2d0))
			NSQrtLD=nrml*( 1.d0-cLD*(1.d0-mu)-dLD*(1.d0-sqrtMu) )	
			
	END FUNCTION NSQrtLD
!========================!
	FUNCTION NSQrtmu(mu, cLD, dLD)
		REAL(DP), INTENT(IN) :: mu, cLD, dLD
		REAL(DP)	:: NSQrtmu, sqrtMu, nrml
		
			sqrtMu=dsqrt(mu)
			nrml=1.d0/(PI_D*(1.d0-cLD/3.d0-dLD*0.2d0))
			NSQrtmu=nrml*( 1.d0-cLD*(1.d0-mu)-dLD*(1.d0-sqrtMu) )	
			
	END FUNCTION NSQrtmu
!========================!
	FUNCTION tanhLD(r)
	IMPLICIT NONE
	REAL(DP), INTENT(IN)	:: r
	REAL(DP)	::	tanhLD
	REAL(DP)	:: mu, term1
	
	mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
	
	term1=4.d0*(1.d0-mu)-2.5d0
	tanhLD=(1.d0-dtanh(term1))/5.098156152d0
		
	END FUNCTION tanhLD
!========================!
	FUNCTION SLD(r, I0, u_SLD)
		REAL(DP), INTENT(IN) :: r, I0, u_SLD
		REAL(DP)	:: SLD, mu
		
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			SLD=I0*((1.d0-u_sLD)+u_SLD*mu)	
			
	END FUNCTION SLD
!========================!
	FUNCTION NLLD(r, u1, u2)
		REAL(DP), INTENT(IN) :: r, u1, u2
		REAL(DP)	:: NLLD, Gamma, mu
		
		
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			!term2=dsqrt(mu)
			Gamma=2.d0*u1/(3.d0-u1)
			NLLD=( (1.d0-Gamma)+1.5d0*Gamma*mu )/PI_d	
			
	END FUNCTION NLLD	
!========================!
	FUNCTION NLD_step(r, u1, u2)
		REAL(DP), INTENT(IN) :: r, u1, u2
		REAL(DP), parameter	::	r_s=0.8d0
		REAL(DP)	:: NLD_step, Gamma, mu,&
					 mu_s, Nrml_s
				
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			IF (r>r_s) mu=0.d0
			
			mu_s=dsqrt(1.d0-r_s)*dsqrt(1.d0+r_s)
			Gamma=2.d0*u1/(3.d0-u1)
			Nrml_s=1.d0+Gamma*mu_s**3
			NLD_step=( (1.d0-Gamma)+1.5d0*Gamma*mu )/(PI_d*Nrml_s)	
			
	END FUNCTION NLD_step
!========================!
	FUNCTION NLD_fringe(r, u1, u2)
		REAL(DP), INTENT(IN) :: r, u1, u2
		REAL(DP), parameter	::	r_f1=0.45d0, r_f2=0.55d0
		REAL(DP)	:: NLD_fringe, Gamma, mu,term1, term2, &
					 mu1, mu2, Nrml_s,r_m,mu_m, I_m
			
			r_m=1.d0!(1.d0+r_f2)/2.d0
			mu_m=dsqrt(1.d0-r_m)*dsqrt(1.d0+r_m)
			I_m=(1.d0-Gamma)+1.5d0*Gamma*mu_m
				
			mu=dsqrt(1.d0-r)*dsqrt(1.d0+r)
			IF (r>r_f1.and.r<r_f2) mu=mu_m
			
			mu1=dsqrt(1.d0-r_f1)*dsqrt(1.d0+r_f1)
			mu2=dsqrt(1.d0-r_f2)*dsqrt(1.d0+r_f2)
						
			Gamma=2.d0*u1/(3.d0-u1)
			term1=(Gamma-1.d0+I_m)*(mu1-mu2)*(mu1+mu2)
			term2=Gamma*(mu2**3-mu1**3)
			Nrml_s=1.d0+term1+term2
			NLD_fringe=( (1.d0-Gamma)+1.5d0*Gamma*mu )/(PI_d*Nrml_s)	
			
	END FUNCTION NLD_fringe			
!=======================!
	FUNCTION I_sin(r, c1)
		REAL(DP), INTENT(IN) :: r, c1
		REAL(DP)	:: I_sin, mu

			mu=1-c1
			I_sin=mu+mu*dsin(4*r)/2
			
	END FUNCTION I_sin
!=======================!	
	SUBROUTINE gauleg(x1,x2,x,w,n)

		REAL(DP), INTENT(IN) :: x1,x2
		INTEGER(I2B), INTENT(IN) ::  n 
		REAL(DP), DIMENSION(1:n), INTENT(OUT) :: x,w
		INTEGER :: its,j,m,i
		INTEGER, PARAMETER :: MAXIT=10
		REAL(DP) :: xl,xm
		REAL(DP), DIMENSION(1:(N+1)/2) :: p1,p2,p3,pp,z,z1
		LOGICAL, DIMENSION(1:(n+1)/2) :: unfinished
		m=(n+1)/2
		xm=0.5D0*(x2+x1)
		xl=0.5D0*(x2-x1)
	DO i=1, m
    	z(i)=dcos(PI_D*( dfloat(i)-0.25D0 )/(n+0.5D0) )
	ENDDO
! WRITE(*,'(5(f10.5,x))')(z(i), i=1,5)
		unfinished=.true.
	DO its=1,MAXIT
    	where (unfinished)
        	p1=1.0
        	p2=0.0
    	end where
    	DO j=1,n
        	where (unfinished)
            	p3=p2
            	p2=p1
            	p1=((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
        	end where
    	ENDDO
    	where (unfinished)
        	pp=n*(z*p1-p2)/(z*z-1.0D0)
        	z1=z
        	z=z1-p1/pp
        	unfinished=(abs(z-z1) > EPS)
    	end where
    	IF (.not. any(unfinished)) exit
	ENDDO
	IF (its == MAXIT+1) THEN
    	WRITE(*,*)'too many iterations in gauleg'
    	STOP 'program terminated by gauleg'
	ENDIF
		x(1:m)=xm-xl*z
		x(n:n-m+1:-1)=xm+xl*z
		w(1:m)=2.0D0*xl/((1.0D0-z**2)*pp**2)
		w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauleg
!=======================!
SUBROUTINE testRAN
	INTEGER	:: i, j, idum
	INTEGER, PARAMETER	:: Nran=20
	REAL(DP), DIMENSION(1:Nran,1:10)	:: Fi
	
	Fi=0.d0
DO i=1, 10
	idum=i	
	DO j=1, Nran
		Fi(j,i)=ran2(idum)
	END DO
END DO

	OPEN(21, file='testRAN.txt')
DO i=1, Nran
	WRITE(21,'(I4,10(f11.8,x))') i,(Fi(i,idum), idum=1,10)
END DO
	CLOSE(21)	
	
		
END SUBROUTINE testRAN
!=======================!
FUNCTION ran2(idum)

	INTEGER, INTENT(INOUT) ::	idum
	INTEGER, PARAMETER	:: IM1=2147483563,IM2=2147483399,IMM1=IM1-1, &
					IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
					IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB 
	REAL(DP) ::	ran2
	REAL(DP), PARAMETER ::	AM=1.d0/IM1, RNMX=1.d0-EPS
!Long period (> 2 × 10^18) random number generator of L’Ecuyer with Bays-Durham shuffle and added safeguards.
! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). 
!Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. 
!RNMX should approximate the largest floating value that is less than 1.
	INTEGER, DIMENSION(1:NTAB)	:: iv
	INTEGER :: idum2,j,k,iy
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/

if (idum.le.0) then 
	idum=max(-idum,1) 
	idum2=idum
	do  j=NTAB+8,1,-1
		k=idum/IQ1
		!Initialize.
		!Be sure to prevent idum = 0.
		!Load the shuffle table (after 8 warm-ups).
		idum=IA1*(idum-k*IQ1)-k*IR1 
		if (idum.lt.0) idum=idum+IM1 
		if (j.le.NTAB) iv(j)=idum
	enddo 
	iy=iv(1) 
endif

	k=idum/IQ1 
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1 
	k=idum2/IQ2 
	idum2=IA2*(idum2-k*IQ2)-k*IR2 
	if (idum2.lt.0) idum2=idum2+IM2 
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum 
	if(iy.lt.1)iy=iy+IMM1 
	ran2=min(AM*iy,RNMX)
	
END FUNCTION ran2
!=======================!	
SUBROUTINE odeint(nf,li,dr,rc, ystart,x1,x2)!,epsme,h1,hmin) 
USE nrutil, ONLY : nrerror,reallocate
!USE ode_path
IMPLICIT NONE
	INTEGER(I1B), INTENT(IN)	:: nf
	REAL(DP), INTENT(IN)	:: li,dr,rc
	REAL(DP), DIMENSION(1:1), INTENT(INOUT) :: ystart 
	REAL(DP), INTENT(IN) :: x1,x2!,epsme,h1,hmin 
	REAL(DP), PARAMETER :: TINY=1.0e-30_dp 
	INTEGER(I4B), PARAMETER :: MAXSTP=10000
	INTEGER(I4B) :: nstp,	nok,nbad,kount
	REAL(DP) :: h,epsme,hdid,hnext,x,xsav, dxsav
	REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal 
	!REAL(DP), DIMENSION(:), POINTER :: xp 
	!REAL(DP), DIMENSION(:,:), POINTER :: yp 
	x=x1
	h=sign(h1,x2-x1)
	epsme=accuracy
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
!if (save_stepsme) then
	!xsav=x-2.0_dp*dxsav
	!nullify(xp,yp)
	!allocate(xp(256)) 
	!allocate(yp(size(ystart),size(xp)))
!end if
do nstp=1,MAXSTP
    do	!me
		call derivs(nf,li,dr,rc, x,y,dydx)
		if (dabs(dydx(1))>eps) exit	!me
	 	x=x+h!epsme	!me
	 	!if (nstp>1) write(*,*)nstp,'dydx(1),x',dydx(1),x	!me
	end do !me	
	!if (nf==1) write(*,*) 'odient:	',nstp, x, y,  dydx
	yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
!	if (save_stepsme .and. (abs(x-xsav) > abs(dxsav))) & !Store intermediate results.
!	call save_a_step
	if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x 
	!If stepsmeize can overshoot, decrease.
	!IF (dabs(dydx)<eps) h=2*h
	call rkqs(nf,li,dr,rc, y,dydx,x,h,epsme,yscal,hdid,hnext)
	if (hdid == h) then
		nok=nok+1 
	else
		nbad=nbad+1 
	end if
	if ((x-x2)*(x2-x1) >= 0.0) then 
		ystart(:)=y(:)
		!if (save_stepsme) call save_a_step
		RETURN 
	end if
!Are we done?
	if (abs(hnext) < hmin) then
			write(*,*) hnext, ystart(1)
			call nrerror('stepsize smaller than minimum in odeint')
	end if		
    h=hnext
end do
	call nrerror('too many steps in odeint')
	
!CONTAINS
	!SUBROUTINE save_a_step 
		!kount=kount+1
	!if (kount > size(xp)) then
		!xp=>reallocate(xp,2*size(xp))
		!yp=>reallocate(yp,size(yp,1),size(xp)) 
	!end if
		!xp(kount)=x 
		!yp(:,kount)=y(:)
		!xsav=x
	!END SUBROUTINE save_a_step 
	
END SUBROUTINE odeint
!=======================!
SUBROUTINE derivs(nf,li,dr,rc, x,y,dydx)
	USE nrtype
	IMPLICIT NONE
		REAL(DP), INTENT(IN)	:: li,dr,rc
		INTEGER(I1B), INTENT(IN)	::	nf
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(1:1), INTENT(IN) :: y 
		REAL(DP), DIMENSION(1:1), INTENT(OUT) :: dydx 
		REAL(DP), DIMENSION(1:1)	:: ystart
		REAL(DP)	::	r, rb, BF2, BF1, IAntg			
		SELECT CASE(nf)
		CASE(0)
			r=dr*x+rc
			!rb=(x-rc)/dr
			BF2=1.d0!SHAPE1D2(x,nf)
			IAntg=IntegA(li,r)
			dydx=r*BF2*IAntg		
		CASE(1,2,3)
			r=dr*x+rc
			!rb=(x-rc)/dr
			BF2=SHAPE1D2(x,nf)
			IAntg=IntegA(li,r)
			dydx=r*BF2*IAntg
		CASE(4)		! I=I(r)
			IAntg=IntegA(li,x)	
			!dydx=x*IAntg*NLLD(x,c1,c2)
			dydx=x*IAntg*NSQrtLD(x,c1,c2)
		CASE(5)		! I=cte
			IAntg=IntegA(li,x)
			dydx=x*IAntg
		CASE(6)		! A=1	
			!dydx=x*NLLD(x,c1,c2)
			dydx=x*NSQrtLD(x,c1,c2)
		CASE(12)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*I_FEM_modified(x)							
		CASE(14)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*I_FEM(x)
		CASE(15)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*I_FEM_2Nodes(x)
		CASE(16)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*NLD_step(x, c1, c2)
		CASE(17)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*NLD_fringe(x, c1, c2)
		CASE(18)
			IAntg=IntegA(li,x)
			dydx=x*IAntg*tanhLD(x)
	
		!CASE(16,17,18)
			!r=dr*x+rc
			!BF1=SHAPE1D1(x,nf-15_I1B)
			!ystart=0
			!CALL odient(nf-15_I1B,r,dr,rc, ystart,-1,1)
			!dydx=BF1*ystart(1)		
			!write(*,*) li, x, I_FEM(x), IAntg, dydx		
		END SELECT	
		!IF (dydx(1)<1.d-14.and.nf>2.and.nf<5) write(*,*) 'l,	r:	', li, x
		!IF (dydx(1)<1.d-14.and.nf==1) write(*,*) nf, 'l,	rb,	V:	',li, x, BF2
		
END SUBROUTINE derivs
!=======================!
FUNCTION I_FEM_modified(r)
	IMPLICIT NONE
	REAL(DP)	::	I_FEM_modified
	REAL(DP), INTENT(IN) :: r
	
	REAL(DP)	::	r1, r2, rb, dr
	INTEGER	::	i, dim1, inode

	dim1=size(rnodes)

	DO i=1, dim1
		IF (r<rnodes(i)) EXIT		
	ENDDO	
	inode=i
	
	IF (inode>1) THEN
		r2=rnodes(inode)
		r1=rnodes(inode-1)
		dr=r2-r1
		rb=2*(r-r1)/dr-1
		I_FEM_modified=SHAPE1D2(rb,1_I1B)*Inodes(inode-1)+SHAPE1D2(rb,2_I1B)*Inodes(inode)
		!print*, INODE, r, I_FEM_modified
	ELSE	!happens just when r=(r1~0)
		I_FEM_modified=Inodes(1)
		print*, 'inode<2',INODE, r, I_FEM_modified
	ENDIF
					   
END FUNCTION I_FEM_modified
!=======================!
FUNCTION I_FEM(r)
!REAL(DP), DIMENSION(:), INTENT(IN) :: Inodes
REAL(DP), INTENT(IN) :: r
REAL(DP)	:: I_FEM 
REAL(DP)	:: rb, sum1, r1, r2, term
INTEGER(I1B) :: id
INTEGER(I2B) :: i, j, Nel, n

n=size(Inodes)
!PRINT*, Nr
!PRINT*, Nglobr(Nr,1_I1B)
DO i=1, Nr
!PRINT*, Nglobr(i, Nd2)
	IF (Nd2>1) THEN
   		r2=rv(Nglobr(i,Nd2))
   	ELSE
   		r2=rv(i+1)
   	END IF	
  
  IF (r<=r2) THEN
  	Nel=i
  	EXIT
  END IF	
  	
END DO

IF (Nd2>1) THEN
	r1=rv(Nglobr(Nel,1_I1B))
	rb=2*(r-r1)/delr(Nel)-1
!WRITE(*,*) 'dr:',delr(Nel), '	r1, r2:',r1, r2,'	r:',r ,'I_FEM function:	rb:',rb
	sum1=0
	DO id=1, Nd2
		i=Nglobr(Nel,id)
		term=SHAPE1D2(rb,id)*Inodes(i)
	!WRITE(*,*) Nel, i, term
		sum1=sum1+term
	END DO
	I_FEM=sum1
ELSE
	I_FEM=Inodes(i)
END IF	

END FUNCTION I_FEM
!=======================!
FUNCTION I_FEM_2nodes(r)
IMPLICIT NONE
	REAL(DP), INTENT(IN) :: r
	REAL(DP)	:: I_FEM_2Nodes
	REAL(DP)	:: rb, r1, r2, V1, V2
	INTEGER(I1B) :: id
	INTEGER(I2B) :: i, Nel, n

	n=size(Inodes)
	
	DO i=1, n-1  
  		IF (r<=rnodes(i+1)) THEN
  			Nel=i
  			EXIT
  		END IF	  	
	END DO

	IF (Nel<n) THEN	
		r1=rnodes(Nel)
		r2=rnodes(Nel+1)
		rb=2.d0*(r-r1)/(r2-r1)-1.d0	
	
		V1=SHAPE1D1(rb,1_I1B)
		V2=SHAPE1D1(rb,2_I1B)
		I_FEM_2Nodes = Inodes(Nel)*V1+Inodes(Nel+1)*V2
	ELSE
		I_FEM_2Nodes=INodes(n)	
	ENDIF
	
END FUNCTION I_FEM_2nodes
!==================================!
	FUNCTION SHAPE1D1(xb, id)
		REAL(DP), INTENT(IN)     :: xb
		INTEGER(I1B), INTENT(IN) :: id
		REAL(DP)	:: SHAPE1D1
	SELECT CASE(Nd1)
	CASE(2)
		SELECT CASE(id)
			CASE(1)
				SHAPE1D1=(1.d0-xb)/2.d0
			CASE(2)
				SHAPE1D1=(1.d0+xb)/2.d0
		END SELECT
	CASE(3)
		SELECT CASE(id)
			CASE(1)
				SHAPE1D1=xb*(xb-1.d0)/2.d0
			CASE(2)
				SHAPE1D1=(1.d0-xb)*(1.d0+xb)
			CASE(3)
				SHAPE1D1=xb*(1.d0+xb)/2.d0
		END SELECT		
	END SELECT		
		
	END FUNCTION SHAPE1D1
!==================================!
	FUNCTION SHAPE1D2(xb, id)
		REAL(DP), INTENT(IN)     :: xb
		INTEGER(I1B), INTENT(IN) :: id
		REAL(DP)	:: SHAPE1D2
	SELECT CASE(Nd2)
	CASE(2)
		SELECT CASE(id)
			CASE(1)
				SHAPE1D2=(1.d0-xb)/2.d0
			CASE(2)
				SHAPE1D2=(1.d0+xb)/2.d0
		END SELECT
	CASE(3)
		SELECT CASE(id)
			CASE(1)
				SHAPE1D2=xb*(xb-1.d0)/2.d0
			CASE(2)
				SHAPE1D2=(1.d0-xb)*(1.d0+xb)
			CASE(3)
				SHAPE1D2=xb*(1.d0+xb)/2.d0
		END SELECT		
	END SELECT		
		
	END FUNCTION SHAPE1D2
!==================================!
SUBROUTINE rkqs(nf,li,dr,rc, y,dydx,x,htry,epsme,yscal,hdid,hnext) 
USE nrutil, ONLY : assert_eq,nrerror
!USE nr, ONLY : rkck
IMPLICIT NONE
	INTEGER(I1B), INTENT(IN)	:: nf
	REAL(DP),	INTENT(IN)	:: li,dr,rc
	REAL(DP), DIMENSION(1:1), INTENT(INOUT) :: y 
	REAL(DP), DIMENSION(1:1), INTENT(IN) :: dydx,yscal 
	REAL(DP), INTENT(INOUT) :: x
	REAL(DP), INTENT(IN) :: htry,epsme
	REAL(DP), INTENT(OUT) :: hdid,hnext
	INTEGER(I4B) :: ndum
	REAL(DP) :: errmax,h,htemp,xnew,hme
	REAL(DP), DIMENSION(size(y)) :: yerr,ytemp
	REAL(DP), PARAMETER :: SAFETY=0.9_dp,PGROW=-0.2_dp,PSHRNK=-0.25_dp,&
						ERRCON=1.89e-4

	ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
	h=htry !Set stepsmeize to the initial trial value. 
	hme=htry
	do
		call rkck(nf,li,dr,rc, y,dydx,x,h,ytemp,yerr) !Take a step.
		 
		errmax=maxval(abs(yerr(:)/yscal(:)))/epsme !Evaluate accuracy.
		
		if (errmax <= 1.0) exit !Step succeeded. 
		
		htemp=SAFETY*h*(errmax**PSHRNK) !Truncation error too large, reduce stepsize. 		
		h=sign(max(abs(htemp),0.1_dp*dabs(h)),h) !No more than a factor of 10.
		xnew=x+h
		
		if (dabs(xnew - x)<eps) then
			write(*,*) 'dydx, x: ',dydx, x
			call nrerror('stepsmeize underflow in rkqs')
		end if	
		
	end	do
	
	if (errmax > ERRCON) then
		hnext=SAFETY*h*(errmax**PGROW) 
	else
		hnext=5.0_dp*h 
	end if
	hdid=h
	x=x+h
	y(:)=ytemp(:)
END SUBROUTINE rkqs
!-----------------------!
SUBROUTINE rkck(nf,li,dr,rc, y,dydx,x,h,yout,yerr) 
USE nrutil, ONLY : assert_eq 
IMPLICIT NONE
	REAL(DP),	INTENT(IN)	:: li,dr,rc
	INTEGER(I1B), INTENT(IN)	:: nf
	REAL(DP), DIMENSION(1:1), INTENT(IN) :: y,dydx 
	REAL(DP), INTENT(IN) :: x,h
	REAL(DP), DIMENSION(1:1), INTENT(OUT) :: yout,yerr 
	INTEGER(I4B) :: ndum
	REAL(DP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp 
	REAL(DP), PARAMETER :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
				A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,& 
				B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,& 
				B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,& 
				B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,& 
				B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,& 
				B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,& 
				C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,& 
				C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,& 
				DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,& 
				DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp
				
	ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
	
		ytemp=y+B21*h*dydx
		
	call derivs(nf,li,dr,rc, x+A2*h,ytemp,ak2)
		ytemp=y+h*(B31*dydx+B32*ak2)
		
	call derivs(nf,li,dr,rc, x+A3*h,ytemp,ak3)
		ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
		
	call derivs(nf,li,dr,rc, x+A4*h,ytemp,ak4)
		ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
		
	call derivs(nf,li,dr,rc, x+A5*h,ytemp,ak5) 
		ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
		
	call derivs(nf,li,dr,rc, x+A6*h,ytemp,ak6) !Sixth step. 
		yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6) !Accumulate increments with proper weights. 
		
		yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
		
!Estimate error as difference between fourth and fifth order methods.
END SUBROUTINE rkck

!============================!
	SUBROUTINE MblFit(lmin, tmax, max_sm, minMbl, dMbl_n, dMbl_p, minb, db_n, db_p, minChi2, DoF)
	IMPLICIT NONE	
		REAL(DP), INTENT(IN)	:: 	lmin, tmax, max_sm
		REAL(DP), INTENT(OUT)	:: 	minMbl, dMbl_n, dMbl_p, minb, db_n, db_p, minChi2
		INTEGER, INTENT(OUT)	::	DoF
		INTEGER(I1B)	:: ierror
		INTEGER			:: i,j,k,nb, nmb
		INTEGER(I2B)	:: nData, np
		REAL(DP)		:: mt, dm, db, sm, mins, maxs, res0, mbl1, Chi2 &
						  , At, At0, delm, maxl, minl,b0,t, delb, &
						  dChi2, ntE, mint, maxt, lData, ave_sm
		REAL(DP), ALLOCATABLE, DIMENSION(:)	:: mtLc, lLc, sigma, res, &
										       FIl0, mI0, Ap
		CHARACTER(LEN=130)	:: FNAM
	
		!WRITE(*,*)'	Sub. MblFit:'

		!DO ic=1, Ncurve
		 
	 	WRITE(FNAM,'(a,a,a)') trim(DATA_PATH),trim(event(ic)%event_fldr), &
	 	 					   trim( event(ic)%tel_file(nlc) )
	 	 					  
	  	OPEN(301, file=FNAM)
	  	!print *,  FNAM,'file opened to obtain blending and m_bl.'
	  	j=0
	  	mins=1.d1
	  	maxs=0.d0
	  	maxl=0.d0
	  	minl=1.d4
	  	mint=1.d4
	  	maxt=0.d0
	  	nData=0
	  	ave_sm=0.d0
	  	DO
	  		READ(301,*,iostat=ierror) t, mt, sm, res0
	  		IF (ierror<0) EXIT
	  			ntE=(t-Event(ic)%t0)/Event(ic)%tE
	  			lData=dsqrt( Event(ic)%p2 + (t-Event(ic)%t0)**2/Event(ic)%ts2 )
	  		IF (res0<3.d0*sm.and.lData>lmin.and.ntE<tmax.and.sm<max_sm) THEN
	  			ntE=dabs(t-event(ic)%t0)/event(ic)%tE
	  			mint=min(ntE, mint)
	  			maxt=max(ntE,maxt)
	  			ave_sm=ave_sm+sm
	  			j=j+1
	  		ENDIF
	  		nData=nData+1	
	  	ENDDO
	  	!ave_sm=ave_sm/j
	  	np=j
	  	ALLOCATE(mtLc(1:np), lLc(1:np), sigma(1:np), FIl0(1:np))
	  	ALLOCATE( mI0(1:np), res(1:np) )
	  	REWIND(301)
	  	j=0
	  	DO i=1, nData
	  		READ(301,*) t, mt, sm, res0
	  		ntE=(t-Event(ic)%t0)/Event(ic)%tE
	  		lData=dsqrt( Event(ic)%p2 + (t-Event(ic)%t0)**2/Event(ic)%ts2 )
	  		IF (res0<3.d0*sm.and.lData>lmin.and.ntE<tmax.and.sm<max_sm) THEN
	  			mins=min(sm,mins)
	  			maxs=max(sm,maxs)
	  			j=j+1
	  			mtLc(j)=mt
	  			sigma(j)=sm
	  			res(j)=res0
	  			lLc(j)=lData
	  			maxl=MAX(lLC(j),maxl)
	  			minl=MIN(lLC(j),minl)
	  		ENDIF
	 	 ENDDO	
	  	CLOSE (301)
	  	PRINT*, ic,'Minimum and Maximum sigma_m=', mins, maxs
	  	!PRINT*,'Minimum and Maximum of |t-t0|/tE=', mint, maxt

		c1=Event(ic)%u1fit(nlc)
		c2=0.d0
	  	FIl0=mFLUX(rmin, llc, np)	!A for I=I(r)
		!first guess for m_b
	  	mI0=-2.5d0*dlog10(FIl0)+Event(ic)%mbl		
	  	!it is exact for blending=1, blending=(F0)/(F0+FB)
	  	
		!------- finding χ^2 for first geuss of m_b and blending=1 ------!
		DoF=np-2
		Chi2=0.d0	  	
	  	DO j=1, np	  	
	  		Chi2=Chi2+( (mtLc(j)-mI0(j) )/sigma(j))**2
	  	END DO	  		  	
	  	IF (np<2) THEN
	  		write(*,*) 'Just',np,' points with res<=3.sigma and t<',tmax,' found.'	  		
	  		STOP
	  	END IF
	  	
	  	!WRITE(*,*) 'No. of points used to find blending and m_bl, with res<3sigma:',np	
	  	minChi2=Chi2/(np-2)
	  	!WRITE(*,*) '!!!!	first geuss: b=1, m_bl=',Event(ic)%mbl,'	1-χ^2=',  &
	  				!mindChi2, '	u_1= ',Event(ic)%u1fit(nlc)
	  					  	
	  	minMbl=Event(ic)%mbl
		minb=1.d0
		  			  	 	  	 
	  	delm=.5d0 
	  	delb=2.d0
	 	dm=1.d-3
	 	db=1.d-3
        nmb=nint(delm/dm)
        nb=nint(delb/db)
        
      	DO k=1, nmb+1
       		mbl1=Event(ic)%mbl+(k-1)*dm-delm/2 	  			    
	   		DO j=1, nb+1
	   			b0=(j-1)*db
	  			mI0=Mbl1-2.5d0*dlog10(b0*FIl0+1.d0-b0)
	  			Chi2=0.d0
	  			DO i=1, np
	  			 	Chi2=Chi2+((mtLc(i)-mI0(i))/sigma(i))**2  		
	  			END DO	  		
	  			Chi2=Chi2/(np-2)
	  			!dChi2=dabs(1.d0-Chi2)
	  			!IF (mindChi2>dChi2) THEN
	  			IF (minChi2>Chi2) THEN
	  				!write(*,*) j,'Chi2, min Chi2',Chi2,'m_b, b', mbl1,b0
	  				!mindChi2=dChi2
	  				minb=b0
	  				minMbl=mbl1
	  				minChi2=Chi2
	  				!write(*,*) 'dChi2, b:',dChi2, minb
	  			END IF
	   		END DO	
	  	END DO
	  	
		!Now we want to derive 1 sigma of our fit:
      	!first of blending:
       		mbl1=minMbl
       		b0=minb	  			    
	   		DO j=1, nb+1
	   			b0=b0-db
	  			mI0=Mbl1-2.5d0*dlog10(b0*FIl0+1.d0-b0)
	  			Chi2=0.d0
	  			DO i=1, np
	  			 	Chi2=Chi2+((mtLc(i)-mI0(i))/sigma(i))**2  		
	  			END DO	  		
	  			Chi2=Chi2/(np-2)
	  			IF (Chi2>(1.d0+minChi2)) EXIT
	   		END DO
	   		db_n=minb-b0
	   		
       		mbl1=minMbl
       		b0=minb	  			    
	   		DO j=1, nb+1
	   			b0=b0+db
	  			mI0=Mbl1-2.5d0*dlog10(b0*FIl0+1.d0-b0)
	  			Chi2=0.d0
	  			DO i=1, np
	  			 	Chi2=Chi2+((mtLc(i)-mI0(i))/sigma(i))**2  		
	  			END DO	  		
	  			Chi2=Chi2/(np-2)
	  			IF (Chi2>(1.d0+minChi2)) EXIT
	   		END DO	   		
	   		db_p=b0-minb
	
	! Then of m_BL   		
       		b0=minb	
       		Mbl1=minMbl  			    
	   		DO k=1, nmb+1
	   			Mbl1=Mbl1-dm
	  			mI0=Mbl1-2.5d0*dlog10(b0*FIl0+1.d0-b0)
	  			Chi2=0.d0
	  			DO i=1, np
	  			 	Chi2=Chi2+((mtLc(i)-mI0(i))/sigma(i))**2  		
	  			END DO	  		
	  			Chi2=Chi2/(np-2)
	  			IF (Chi2>(1.d0+minChi2)) EXIT
	   		END DO
	   		dMbl_n=minMbl-Mbl1
	   		
       		b0=minb	
       		Mbl1=minMbl  			    
	   		DO k=1, nmb+1
	   			Mbl1=Mbl1+dm
	  			mI0=Mbl1-2.5d0*dlog10(b0*FIl0+1.d0-b0)
	  			Chi2=0.d0
	  			DO i=1, np
	  			 	Chi2=Chi2+((mtLc(i)-mI0(i))/sigma(i))**2  		
	  			END DO	  		
	  			Chi2=Chi2/(np-2)
	  			IF (Chi2>(1.d0+minChi2)) EXIT
	   		END DO	   		
	   		dMbl_p=Mbl1-minMbl
	   			   			
		!WRITE(*,*) ''//achar(27)//'[1m''Object No.',ic,' m_b by fitting:',minMbl,'±',dm,'χ^2=',minChi2,' No. of data points=',np	  
		!WRITE(*,*) ''//achar(27)//'[1m','of ',nData,' Object No.',ic,' blending:',minb,'±',db,''//achar(27)//'[0m'

		DEALLOCATE(llc)
	  	ALLOCATE(llc(1:4), Ap(1:4)) 
	  	llc(1)=Event(ic)%p
	  	llc(2)=1.d0
	  	llc(3)=minl
	  	llc(4)=maxl
	  	!c1=Event(ic)%u1fit(nlc)
	  	!c2=0.d0
	  	!Ap=mFLUX(rmin, llc, 4_I2B)

		!WRITE(*,*) ''//achar(27)//'[1m''	(l=p, A_p):',llc(1),Ap(1),'	(l=1, A):',llc(2),Ap(2),'	(l=min_l, A):'&
					!,llc(3),Ap(3),'(l=max_l, A):',llc(4), Ap(4), ''//achar(27)//'[0m'
			
	
		!IF (MINChi2>3) PRINT*, " Chi2 of b and m_BL fit > 3" 
		 
		Event(ic)%mbl=minMbl
		Event(ic)%blndng=minb
												
		DEALLOCATE(mtLc, lLc, sigma, FIl0, mI0)

		!END DO

	END SUBROUTINE MblFit	
!============================!
	
END MODULE lensFs2

