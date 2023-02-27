MODULE ELLIPTIC
USE nrtype
USE nrutil, ONLY : assert
IMPLICIT NONE

CONTAINS

FUNCTION ellf_d(phi,ak) 
	REAL(DP), INTENT(IN) :: phi,ak 
	REAL(DP) :: ellf_d
!⋆⋆⋆
!Legendre elliptic integral of the 1st kind F (φ, k), evaluated using Carlson's function RF .
!The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ ksinφ ≤ 1. 
	REAL(DP) :: s
	s=dsin(phi) 
	ellf_d=s*rf_d(dcos(phi)**2,(1.0_dp-s*ak)*(1.0_dp+s*ak),1.0_dp) 
END FUNCTION ellf_d
!----------------------------!
FUNCTION elle_d(phi,ak) 
	REAL(DP), INTENT(IN) :: phi,ak 
	REAL(DP) :: elle_d
!Special Functions
!￼Legendre elliptic integral of the 2nd kind E(φ,k), evaluated using Carlson's functions RD
!andRF. Theargumentrangesare0≤φ≤π/2,0≤ksinφ≤1. 
	REAL(DP) :: cc,q,s
	s=dsin(phi)
	cc=dcos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak) 
	elle_d=s*(rf_d(cc,q,1.0_dp)-((s*ak)**2)*rd_d(cc,q,1.0_dp)/3.0_dp) 
END FUNCTION elle_d
!----------------------------!
FUNCTION ellpi_d(phi,en,ak) 
	REAL(DP), INTENT(IN) :: phi,en,ak 
	REAL(DP) :: ellpi_d
!Seenotetoerf v,p.1094above.
!Legendre elliptic integral of the 3rd kind Π(φ, n, k), evaluated using Carlson's functions RJ and RF . 
!(Note that the sign convention on n is opposite that of Abramowitz and Stegun.) The ranges of φ and 
!k are 0 ≤ φ ≤ π/2, 0 ≤ ksinφ ≤ 1.
	REAL(DP) :: cc,enss,q,s
	s=dsin(phi)
	enss=en*s*s
	cc=dcos(phi)**2
	q=(1.0_dp-s*ak)*(1.0_dp+s*ak) 
	ellpi_d=s*(rf_d(cc,q,1.0_dp)-enss*rj_d(cc,q,1.0_dp,1.0_dp+enss)/3.0_dp) 
END FUNCTION ellpi_d
!----------------------------!
FUNCTION rc_d(x,y)
REAL(DP), INTENT(IN) :: x,y
REAL(DP) :: rc_d
REAL(DP), PARAMETER :: ERRTOL=0.04_dp,TINY=1.69e-38_dp,&
					SQRTNY=1.3e-19_dp,BIG=3.0e37_dp,TNBG=TINY*BIG,& 
					COMP1=2.236_dp/SQRTNY,COMP2=TNBG*TNBG/25.0_dp,& 
					THIRD=1.0_dp/3.0_dp,& 
					C1=0.3_dp,C2=1.0_dp/7.0_dp,C3=0.375_dp,C4=9.0_dp/22.0_dp
!Computes Carlson’s degenerate elliptic integral, RC(x,y). x must be nonnegative 
!and y must be nonzero. If y < 0, the Cauchy principal value is returned. 
!TINY must be at least 5 times the machine underflow limit, BIG at most one-fifth the machine maximum overflow limit.
	REAL(DP) :: alamb,ave,s,w,xt,yt
	
	call assert( (/x >= 0.0,y /= 0.0,x+abs(y) >= TINY,x+abs(y) <= BIG, &
			y >= -COMP1 .or. x <= 0.0 .or. x >= COMP2/),'rc_s') 
if (y > 0.0) then
    xt=x
    yt=y
    w=1.0
else 
	xt=x-y
	yt=-y
	w=dsqrt(x)/dsqrt(xt) 
end if
do
	alamb=2.0_dp*dsqrt(xt)*dsqrt(yt)+yt 
	xt=0.25_dp*(xt+alamb) 
	yt=0.25_dp*(yt+alamb) 
	ave=THIRD*(xt+yt+yt) 
	s=(yt-ave)/ave
	if (dabs(s) <= ERRTOL) exit
end do 
	rc_d=w*(1.0_dp+s*s*(C1+s*(C2+s*(C3+s*C4))))/dsqrt(ave) 
END FUNCTION rc_d
!----------------------------!
FUNCTION rf_d(x,y,z)

	REAL(DP), INTENT(IN) :: x,y,z
	REAL(DP) :: rf_d
	REAL(DP), PARAMETER ::  ERRTOL=0.08_dp,TINY=1.5e-38_dp,BIG=3.0e37_dp,&
							THIRD=1.0_dp/3.0_dp,&
							C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
!Computes Carlson's elliptic integral of the first kind, rf_d(x,y,z). x, y, and z must be nonnegative, 
!and at most one can be zero. TINY must be at least 5 times the machine underflow limit, BIG at most 
!one-fifth the machine overflow limit.
	REAL(DP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt 
	call assert(min(x,y,z) >= 0.0d0, min(x+y,x+z,y+z) >= TINY, &
			max(x,y,z) <= BIG, 'rf_d args') 
	xt=x
	yt=y 
	zt=z 
do
	sqrtx=dsqrt(xt)
	sqrty=dsqrt(yt)
	sqrtz=dsqrt(zt) 
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz 
	xt=0.25_dp*(xt+alamb) 
	yt=0.25_dp*(yt+alamb) 
	zt=0.25_dp*(zt+alamb) 
	ave=THIRD*(xt+yt+zt) 
	delx=(ave-xt)/ave
	dely=(ave-yt)/ave
	delz=(ave-zt)/ave
	if (max(dabs(delx),dabs(dely),dabs(delz)) <= ERRTOL) exit
end do
	e2=delx*dely-delz**2
	e3=delx*dely*delz 
	rf_d=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave) 
END FUNCTION rf_d
!----------------------------!

FUNCTION rd_d(x,y,z)

	REAL(DP), INTENT(IN) :: x,y,z
	REAL(DP) :: rd_d
	REAL(DP), PARAMETER :: ERRTOL=0.05_dp,TINY=1.0e-25_dp,BIG=4.5e21_dp,&
							C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
							C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
!Computes Carlson's elliptic integral of the second kind, rd_d(x,y,z). x and y must be nonnegative, 
!and at most one can be zero. z must be positive. TINY must be at least twice the negative 2/3 power 
!of the machine overflow limit. BIG must be at most 0.1 × ERRTOL times the negative 2/3 power of the machine underflow limit.
	REAL(DP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,& 
				ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
				
	call assert(min(x,y) >= 0.0d0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, & 
			'rd_d args')
	xt=x
	yt=y
	zt=z
	sum=0.0d0
	fac=1.0d0
do
	sqrtx=dsqrt(xt)
	sqrty=dsqrt(yt)
	sqrtz=dsqrt(zt) 
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz 
	sum=sum+fac/(sqrtz*(zt+alamb)) 
	fac=0.25_dp*fac
	xt=0.25_dp*(xt+alamb) 
	yt=0.25_dp*(yt+alamb) 
	zt=0.25_dp*(zt+alamb) 
	ave=0.2_dp*(xt+yt+3.0_dp*zt) 
	delx=(ave-xt)/ave 
	dely=(ave-yt)/ave 
	delz=(ave-zt)/ave
	if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit 
end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_dp*eb
	ee=ed+ec+ec 
	rd_d=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
			+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*dsqrt(ave)) 
END FUNCTION rd_d
!----------------------------!

FUNCTION rj_d(x,y,z,p)

	REAL(DP), INTENT(IN) :: x,y,z,p
	REAL(DP) :: rj_d
	REAL(DP), PARAMETER :: ERRTOL=0.05_dp,TINY=2.5e-17_dp,BIG=9.0e15_dp,&!,TINY=2.5e-13_dp,BIG=9.0e11_dp,&
					C1=3.0_dp/14.0_dp,C2=1.0_dp/3.0_dp,C3=3.0_dp/22.0_dp,& 
					C4=3.0_dp/26.0_dp,C5=0.75_dp*C3,C6=1.5_dp*C4,C7=0.5_dp*C2,& 
					C8=C3+C3
!Computes Carlson's elliptic integral of the third kind, RJ (x, y, z, p). x, y, 
!and z must be nonnegative, and at most one can be zero. p must be nonzero. If p < 0, the Cauchy
!principal value is returned. TINY must be at least twice the cube root of the machine
!underflow limit, BIG at most one-fifth the cube root of the machine overflow limit. 
	REAL(DP) :: a,alamb,alpha,ave,b,bet,delp,delx,&
				dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
				sm,tau,xt,yt,zt
	call assert(min(x,y,z) >= 0.0d0, min(x+y,x+z,y+z,abs(p)) >= TINY, &
						max(x,y,z,abs(p)) <= BIG, 'rj_d args')
	sm=0.0d0
	fac=1.0d0
if (p > 0.0d0) then
    xt=x
    yt=y
    zt=z
    pt=p
else 
	xt=min(x,y,z)
	zt=max(x,y,z) 
	yt=x+y+z-xt-zt 
	a=1.0_dp/(yt-p) 
	b=a*(zt-yt)*(yt-xt) 
	pt=yt+b 
	rho=xt*zt/yt 
	tau=p*pt/yt
end if 

do
	sqrtx=dsqrt(xt)
	sqrty=dsqrt(yt)
	sqrtz=dsqrt(zt) 
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz 
	alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2 
	bet=pt*(pt+alamb)**2
	sm=sm+fac*rc_d(alpha,bet) 
	fac=0.25_dp*fac 
	xt=0.25_dp*(xt+alamb) 
	yt=0.25_dp*(yt+alamb) 
	zt=0.25_dp*(zt+alamb) 
	pt=0.25_dp*(pt+alamb) 
	ave=0.2_dp*(xt+yt+zt+pt+pt) 
	delx=(ave-xt)/ave 
	dely=(ave-yt)/ave 
	delz=(ave-zt)/ave 
	delp=(ave-pt)/ave
	if (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL) exit 
end do
	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.0_dp*ec
	ee=eb+2.0_dp*delp*(ea-ec) 
	rj_d=3.0_dp*sm+fac*(1.0_dp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
					+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave)) 
	if (p <= 0.0d0) rj_d=a*(b*rj_d+3.0_dp*(rc_d(rho,tau)-rf_d(xt,yt,zt))) 
END FUNCTION rj_d
!----------------------------!
END MODULE ELLIPTIC