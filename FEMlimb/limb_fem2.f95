MODULE limb_fem2
USE Mesh
USE lensFs2
USE femFs2

IMPLICIT NONE

	INTEGER	:: lenfld, used_sets
!(ic, nLc):(2,1) (2,9) (3,2) (3,5) (4,1) (5,5) (6,4) (6,5) 
!           (6,6) (6,7) (8,6) (9,2)  (9,9)
!(2,1,Nran100)(2,9)(3,2)(3,5 not smooth) (5,5) (7,7) (8,6) (9,2) (9,9)

!good results:  (7,7):OB04254_SAAO:fig.1 (2,1):MB07233_OGLE:fig.2 (2,9)_MOA:fig.3 !it has some data points with res0>3sm
! 				 (4,1):MB10436_MOA:fig.4 (5,5):MB11093_Canopus:fig.5 
!				(8,6):OB110990_Pico:fig.6 (9,9):OB111101_LT:fig.7

	INTEGER, DIMENSION(1:2,1:8)	:: good_result=reshape( (/7,7,&
														  2,1,&
														  2,9,&
														  4,1,&
														  5,5,&
														  8,6,&
														  9,9,&
														  10,1/), (/2,8/)) !these are (ic, NLc) pairs.
													  														  
	REAL(DP), DIMENSION(1:8)	:: delta_trsh=(/0.4d0, 0.5d0, 0.5d0, 0.5d0, 0.4d0, 0.3d0, 0.3d0, 0.5d0/), &
										lmbd_choice=(/1.d-1, 1.d2, 5.d-1, 5.d2, 1.5d-2, 3.d-1, 1.d-1 ,3.d-2/), &
										 dsprn_bnd=(/3.d0,2.d0,1.01d0,2.d0,1.5d0,1.05d0,1.1d0,2.5d0/) ! ic=10, Ndata=20 -> dsprn_bnd=4.5
										!final: dsprn_bnd=(/1.1d0,2.d0,1.01d0,2.d0,1.5d0,1.05d0,1.1d0,2.5d0/) ! ic=10, Ndata=20 -> dsprn_bnd=4.5
										
						!lmbd_choice=(/8.d-2, 5.d0, 5.d-1, 5.d2, 1.d-3, 1.d-1, 1.d-1 ,1.d-2/)
						!lmbd_choice=(/0.1d0, 0.1d0, 3.16d-1, 3.1623d1, 3.2d-2, 1.d0, 3.d-3, 1.d-2/)	!SQrt optimization and run2, the first LC was analysed with r-histograms									
						!lmbd_choice=(/8.d-2,1.d0,5.d-1,1.d1,5.d-3,1.d-2,1.d-1, 1.d-2/)
								   
								   
								   !delta_trsh=(/0.4d0, 0.5d0, 0.5d0, 0.5d0, 0.4d0, 0.5d0, 0.3d0, 0.5d0/)
	REAL(DP), parameter	::	dlmbd=5.d-2								   
	INTEGER, PARAMETER	:: Nlambda=160, jlambda=5 !minimum λ =10^(-j_λ), 
												 !maximum λ =10^(-j_λ+N_λ*dlmbd)
	
	!REAL(DP), DIMENSION(1:Nlambda)	::	lmbd_bnd=(/ 1.d-6, 3.5d-6, 1.d-5, 3.5d-5, 1.d-4, 3.5d-4, 1.d-3/)!, 3.5d-3, 1.d-2/)
	!REAL(DP), DIMENSION(1:Nlambda)	::	lmbd_bnd=(/1.d-4, 3.5d-4, 1.d-3, 3.5d-3, 1.d-2, 3.5d-2, 1.d-1, 3.5d-1/)
	!REAL(DP), DIMENSION(1:Nlambda)	::	lmbd_bnd=(/1.d-3, 3.5d-3, 1.d-2, 3.5d-2, 1.d-1, 3.5d-1, 1.d0, 3.5d0, 1.d1, 3.5d1, 1.d2/)	
													
	INTEGER, DIMENSION(1:8)	::	Nbox_min=(/10,5,15,20,20,20,5,30/), Nbox_max=(/20,20,20,25,25,50,10,60/), &
								Nbins=(/100,100,100,100,100,100,9,20/) !Nbins=(/19,8,15,23,25,18,9,50/)!Nbins=(/50,6,18,22,25,50,8,50/)
	CHARACTER(len=50)	:: FMT1	
	REAL(DP)	::  dsprn_r
	LOGICAL	::	SIM, calmbl, outornot, Istep, NSqrt_mode, talk
CONTAINS

!=======================! 
SUBROUTINE lstsctrng(Xel, Q_ls)
	REAL(DP), DIMENSION(:), INTENT(IN)	::	Xel
	REAL(DP), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)	:: Q_ls
	INTEGER(I2B)	:: sizex, i,j,k,ii,jj, kk
	REAL(DP), DIMENSION(:), ALLOCATABLE	::  xratio, hi
	REAL(DP), DIMENSION(1:3,1:3) :: Q33
	REAL(DP)	::	invh1h2
	sizex=size(Xel)
	ALLOCATE( hi(2:sizex), xratio(2:sizex-1), Q_ls(1:sizex,1:sizex) )
	DO i=2, sizex
		hi(i)=Xel(i)-Xel(i-1)
	ENDDO
		
	DO i=2, sizex-1
		xratio(i)=(Xel(i)-Xel(i-1))/(Xel(i+1)-Xel(i-1))
		!write(*,*) 'lstsctrng sub. s_i:',sigmaI_cntn(i)
	END DO	
	
	IF (sizex<4) THEN
		WRITE(*,*) 'Not enough point to calculate least scattering MATRIX.'
		STOP
	END IF
	Q_ls=0.d0
	
	DO i=2, sizex-1
		invh1h2=1.d0/(hi(i)*hi(i+1))
		Q33=4.d0*invh1h2**2 * Qii(xratio(i))
		DO j=1, 3
			ii=i-2+j
			DO k=1, 3
				jj=i-2+k
				Q_ls(ii,jj)=Q_ls(ii,jj)+Q33(j,k)
			END DO
		END DO
	END DO		
	
END SUBROUTINE lstsctrng
!=======================! 
SUBROUTINE lstsctrng_3(Xel, Q_ls)
	REAL(DP), DIMENSION(:), INTENT(IN)	::	Xel
	REAL(DP), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)	:: Q_ls
	INTEGER(I2B)	:: sizex, i,j,k,ii,jj, kk
	REAL(DP), DIMENSION(:), ALLOCATABLE	::  xratio, hi
	REAL(DP), DIMENSION(1:3,1:3) :: Q33
	REAL(DP)	::	invh1h2!, x3dx32
	
	sizex=size(Xel)
	
	ALLOCATE( hi(2:sizex), xratio(2:sizex-1), Q_ls(1:sizex,1:sizex) )
	
	DO i=2, sizex
		hi(i)=Xel(i)-Xel(i-1)
	ENDDO
		
	DO i=2, sizex-1
		xratio(i)=(Xel(i)-Xel(i-1))/(Xel(i+1)-Xel(i-1))
		!write(*,*) 'lstsctrng sub. s_i:',sigmaI_cntn(i)
	END DO	
	
	IF (sizex<4) THEN
		WRITE(*,*) 'Not enough point to calculate least scattering MATRIX.'
		STOP
	END IF
	Q_ls=0.d0	
	DO i=2, sizex-1
		invh1h2=1.d0/(hi(i)*hi(i+1))
		Q33=4.d0*invh1h2**2 * Qii(xratio(i))
		DO j=1, 3
			ii=i-2+j
			DO k=1, 3
				jj=i-2+k
				Q_ls(ii,jj)=Q_ls(ii,jj)+Q33(j,k)
			END DO
		END DO
	END DO
			
	
	
	
!I_0=0	condition:
	invh1h2=1.d0/(hi(2)**4)
	Q_ls(1,1)=Q_ls(1,1)+invh1h2
	Q_ls(1,2)=Q_ls(1,2)-2.d0*invh1h2
	Q_ls(2,1)=Q_ls(2,1)-2.d0*invh1h2
	Q_ls(2,2)=Q_ls(2,2)+4.d0*invh1h2
!I_{N+1}=0	condition:
	!invh1h2=1.d0/(hi(sizex-1)**4)
	!Q_ls(sizex-1,sizex-1)=Q_ls(sizex-1,sizex-1)+4.d0*invh1h2
	!Q_ls(sizex-1,sizex)=Q_ls(sizex-1,sizex)-2.d0*invh1h2
	!Q_ls(sizex,sizex-1)=Q_ls(sizex,sizex-1)-2.d0*invh1h2
	!Q_ls(sizex,sizex)=Q_ls(sizex,sizex)+invh1h2
	
END SUBROUTINE lstsctrng_3
!=======================!
FUNCTION Qii(xratio)
	REAL(DP), DIMENSION(1:3, 1:3) :: Qii
	REAL(DP)	:: xratio
	INTEGER(I2B):: j,k
	
	Qii=0.d0
	
	Qii(1,1)=(xratio-1)**2
	Qii(2,1)=(xratio-1); Qii(2,2)=1
	Qii(3,1)=-xratio*(xratio-1); Qii(3,2)=-xratio; Qii(3,3)=xratio*xratio
	Qii(1,2)=Qii(2,1); Qii(1,3)=Qii(3,1); Qii(2,3)=Qii(3,2)	


END FUNCTION Qii
!=======================! 
SUBROUTINE LC_FEM(FNAME,r_cntn,I_cntn,l_all,A_all,dA_all, A_FEM,Chi2)
!this sub. calculates light curve of I_FEM 
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN)	:: FNAME
REAL(DP), DIMENSION(:), INTENT(IN)	:: r_cntn,I_cntn, l_all, A_all, dA_all
REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)	::	A_FEM
REAL(DP), INTENT(OUT)	::	Chi2

INTEGER		:: i, id, i2, i1, dim1, dim2

REAL(DP)	:: r,IFEM,li, t, A_mdl,A_exact, &
				Cnrml, A, sA, erI
REAL(DP), ALLOCATABLE, DIMENSION(:)	::	Iexact, Aeq
INTEGER(I1B)	:: nf
	
	OPEN(103, FILE=FNAME)
	WRITE(103,'(a)') '#     l      ,    A_FEM   ,  A_all    '
	dim1=size(A_all)
	dim2=size(I_cntn)
	print*, 'dim2, size_r',dim2, size(r_cntn)
	ALLOCATE(A_FEM(1:dim1), Inodes(1:dim2), rnodes(1:dim2))
    Inodes=I_cntn
    rnodes=r_cntn
    nf=12_I1B !->  I_FEM_modified(r) will be use in flux calculations
    Chi2=0.d0
    DO i=1, dim1  	
    	li=l_all(i)
		A_FEM(i)=FI0l(nf,rmin, li)
		Chi2=Chi2+((A_FEM(i)-A_all(i))/dA_all(i))**2
		WRITE(103,'(3(f12.5,x))') li, A_FEM(i), A_all(i)
		!print*,'i,l_i:',i,li,' A_FEM, A_DATA:',A_FEM(i), A_all(i)
		!PRINT*
	END DO
	Chi2=Chi2/dim1
	
	WRITE(103, *)'#	Chi^2=',Chi2,' DoF=',dim1
	CLOSE(103)
	print*,'LC_FEM SUB.: ic, nLc, Chi2: ',ic,nLc,Chi2
	DEALLOCATE( Inodes, rnodes )
END SUBROUTINE LC_FEM							          
!========================!
SUBROUTINE FEM_regularisation( lambda0, setNos, l_all, A_all, dA_all, I_input1,r_cntn, IFEM_1, lmbd_opt1)
 										 
 ! when sim is true, IFEM_1	is the solution for A_all with least error (lmbd_opt1).
 ! when .not.sim IFEM_1 is from A_all with lambda0.								       	         		    
	IMPLICIT NONE	
	REAL(DP), INTENT(IN)	::	lambda0		!we don't use it in case of SIM
	INTEGER , DIMENSION(:, :), INTENT(IN)	:: setNos
	REAL(DP), DIMENSION(:)   , INTENT(IN)   :: l_all, A_all, dA_all, I_input1
	REAL(DP), INTENT(OUT)	::	lmbd_opt1	!we don't need them in case of .not.SIM
	REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: r_cntn, IFEM_1
																											    												         									
	INTEGER , ALLOCATABLE, DIMENSION(:) ::  indx, index		
	REAL(DP), ALLOCATABLE, DIMENSION(:) ::   A_DATA1, dA1, IFEM1, li, uMITA1, mu, MImA, & 
			    		 					 MI, QI, I_cntn1, dsprn, var_dsprn, Chi2, err1 
			    		 					   		 	    						           			                                         			           			                                         																
	REAL(DP), ALLOCATABLE, DIMENSION(:, :)	:: uMI, Mstiff, Mlambda, uMIT, M2, dummy, &
								      	       Q_ls, I_lambda1, mean_I2	
								      	       								   	
	REAL(DP)	:: d, t, r, rb, V1, V2, wR2, lsq1, norm_dA, varI, varI_n, varI_p, maxvarI_n, &
			       AFEMmA, lambda, lambda1, pwr, min_dsprn, min_dchi2, max_dChi2, minvarI_p, &
			       dChi2, min_Chi2, all_varI, rDsprn				   
	    	    	
	INTEGER :: ierror, i, j, k, i1, idum1, iset, Ni, Nj, nsetsf, used0, Npoints, Nls, &
											jp, jeq, minjp, minjeq, i_opt, jn, maxjn
											
	INTEGER(I2B)	::	DoF		   						   
	INTEGER(I2B), PARAMETER	:: ng=2_I2B	
	CHARACTER(LEN=100)	::	Flambda

!-----------------------------------------!

	!ALLOCATE(ri(1:ng), wi(1:ng))
	!CALL gauleg(-1.d0, 1.d0, ri, wi, ng)
!IF (SIM) THEN	
!	WRITE(Flambda,'(a,I1,a,I1,a)')'lambda-dAFEM-dA-res-ic',ic,'nLc',nLc,'.dat'	
!	OPEN(201,file=Flambda)
!ENDIF
	Nsetsf=size( setNos(:,1:1) )
	Ni=size(setNos(1:1,:))
!print*,'SUBROUTINE FEM_regularisation_step:Nsetsf,ni=',Nsetsf,ni	
	lmbd_opt1=0.d0	
	IF (I1eqI2) THEN
		Nj=Ni+1
	ELSE
		Nj=Ni
	ENDIF
	
	Nls=size(l_all)
	Npoints=Nls+2
!print*,'SUBROUTINE FEM_regularisation_step: Nls=',Nls	
	ALLOCATE( dsprn(1:Nlambda+1), var_dsprn(1:Nlambda) )
	ALLOCATE( Chi2(1:Nlambda+1), err1(1:Nlambda+1) )
	ALLOCATE( r_cntn(1:Npoints), IFEM_1(1:Npoints) )
	ALLOCATE( I_cntn1(1:Npoints) )
	IFEM_1=0.d0	
	r_cntn=0.d0	
	DO i=2, Npoints-1
		r_cntn(i)=l_all(i-1)!(i-1)*deltar
	ENDDO
	r_cntn(Npoints)=1.d0
		
! print*,'Nsetsf, ni,nj',Nsetsf,Ni,Nj
 	Allocate( li(1:Ni), A_DATA1(1:Nj), IFEM1(1:Nj) )
	ALLOCATE( I_lambda1(1:Nlambda+1,1:Npoints))
	ALLOCATE( uMI(1:Nj, 1:Nj), uMIT(1:Nj,1:Nj),M2(1:Nj, 1:Nj) )		
	ALLOCATE( Mstiff(1:Ni, 1:Nj), Mlambda(1:Nj, 1:Nj)  , dummy(1:Nj,1:Nj))
	ALLOCATE( MI(1:Nj), indx(1:Nj), uMITA1(1:Nj),dA1(1:Nj) )
	
	I_cntn1=0.d0
	I_lambda1=0.d0
	lambda1=1.d-1**jlambda
	DO iset=1, Nsetsf	
		!print*,'iset:',iset
		li=0.d0
		!dA1=0.d0
		A_DATA1=0.d0
		
		DO i=1, Ni
			j=setNos(iset,i)
			!IF (j==Nls) print*,iset,i,'last data point was used:',j
			IF (j>Nls.or.j<=0) THEN
				write(*,'(2I3,a,I3,a,2I3)') i,Ni,' FEM_regularisation_step SUB. Alert!!! Bad number of data in setNos'&
						  ,j,'!  data seri No.:',iset,Nsetsf
				print*,' SIM:',SIM,'Istep',Istep
				!print*, setNos(1,:)
				print*,setNos(iset,:)
				STOP
			ENDIF
				
			li(i)=l_all(j)
			IF (I1eqI2) THEN
				!dA1(i+1)=dA_all(j)
				A_DATA1(i+1)=A_all(j)			
			ELSE
				!dA1(i)=dA_all(j)
				A_DATA1(i)=A_all(j)				
			ENDIF	
								
		ENDDO
		!IF ( I1eqI2 )	dA(1) = dsqrt( dA(2)**2 + dA(3)**2 )	!because the first row is for the I_2-I_1=0 constraint
					
		CALL mesh_l(rmin, li)
		ALLOCATE(mu(1:Nj))
		mu(1)=1.d0
		mu(2:Nj-1)=dsqrt(1.d0-rv(2:Nj-1))*dsqrt(1.d0+rv(2:Nj-1))
		mu(Nj)=0.d0

		Mstiff=0.d0				
		IF (Nd1==2.or.Nd1==3) THEN 
			CALL GlM_limbr2( 10_I2B, Mstiff)
		ELSE IF (Nd1==1.and.Nd2>1) THEN
			CALL GlM_limbr1( Mstiff)
		ELSE IF (Nd1==1.and.Nd2==1)	THEN
			CALL GlM_limbr0( Mstiff)	
		END IF
		uMI=0.d0
		IF (I1eqI2) THEN							
			uMI(2:Nj,1:Nj)=Mstiff		
			uMI(1,1)=1.d0
			uMI(1,2)=-1.d0
		ELSE
			uMI=Mstiff		
		ENDIF
		
		uMIT=0.d0
		uMITA1=0.d0			
		DO i=1, Nj
			DO j=1, Nj
				uMIT(j,i)=uMI(i,j)
				uMITA1(i)=uMITA1(i)+A_DATA1(j)*uMI(j,i)   !!!!!  MT.A = AT.M
				!uMITdA1(i)=uMITdA1(i)+dA1(j)*uMI(j,i)   !!!!!  MT.dA = dAT.M
			ENDDO
		ENDDO
		M2=matmul(uMIT, uMI)
		
		!norm_dA=dot_product(dA,dA)	
					
		CALL lstsctrng(mu,Q_ls)
		!CALL lstsctrng(rv,Q_ls)		
		!CALL lstsctrng_3(rv,Q_ls)
		!CALL lstsctrng_3(mu,Q_ls)

IF (SIM) THEN
	DO i1=1, Nlambda+1
		!lambda=nint(dsqrt(1.d1**i1))/(1.d1**jlambda)
		pwr=(i1-1)*dlmbd
		lambda=1.d1**pwr*lambda1
		!print*,i1,lambda,lambda1
		Mlambda=M2+lambda*Q_ls					! To solve the (MT.M+λQ).I=AT.M
		dummy=Mlambda
		CALL ludcmp(dummy,Nndr,indx,d)		
		IFEM1=uMITA1		
		CALL lubksb(dummy,Nndr,indx,IFEM1)
		
		!dummy=Mlambda
		!CALL ludcmp(dummy,Nndr,indx,d)		
		!dI=uMITdA1		
		!CALL lubksb(dummy,Nndr,indx,dI)

		DO k=2, Npoints-1			
			r=r_cntn(k)						
			j=1
			DO WHILE (r>rv(j)) 
				j=j+1
			END DO
			j=j-1
				
			IF (abs(r)<eps)	j=1
			rb=2.d0*(r-rv(j))/delr(j)-1.d0	
			IF (rb>1.d0.or.rb<-1.d0) THEN
				PRINT*, ' FEM_regularisation_step SUB. bad rb &
			 	& element:',j,' rb=',rb, ' r=',r
				STOP
			END IF
				
			V1=SHAPE1D1(rb,1_I1B)
			V2=SHAPE1D1(rb, 2_I1B)
			I_cntn1(k) = IFEM1(j)*V1+IFEM1(j+1)*V2
			!dI_FEM(k) =  (dI(j)*V1)**2 + (dI(j+1)*V2)**2 
		END DO !Npoints
		
		I_cntn1(1)=IFEM1(1)
		I_cntn1(Npoints)=IFEM1(Nj)
		!dI_FEM(1)=dabs(dI(1))
		!dI_FEM(Npoints)=dabs(dI(Nj))
		
		I_lambda1(i1,1:Npoints)=I_lambda1(i1,1:Npoints)+I_cntn1(1:Npoints)
		!dI_lambda1(i1,1:Npoints)=dI_lambda1(i1,1:Npoints)+dI_FEM(1:Npoints)
	ENDDO !Nlambda
ELSE
	lambda=lambda0
	
	Mlambda=M2+lambda*Q_ls					! To solve the (MT.M+λQ).I=AT.M
	dummy=Mlambda
	CALL ludcmp(dummy,Nndr,indx,d)		
	IFEM1=uMITA1		
	CALL lubksb(dummy,Nndr,indx,IFEM1)
		
	!dummy=Mlambda
	!CALL ludcmp(dummy,Nndr,indx,d)		
	!dI=uMITdA1		
	!CALL lubksb(dummy,Nndr,indx,dI)

	DO k=2, Npoints-1			
		r=r_cntn(k)						
		j=1
		DO WHILE (r>rv(j)) 
			j=j+1
		END DO
		j=j-1
				
		IF (abs(r)<eps)	j=1
		rb=2.d0*(r-rv(j))/delr(j)-1.d0	
		IF (rb>1.d0.or.rb<-1.d0) THEN
			PRINT*, ' FEM_regularisation_step SUB. bad rb &
		 	& element:',j,' rb=',rb, ' r=',r
			STOP
		END IF
				
		V1=SHAPE1D1(rb,1_I1B)
		V2=SHAPE1D1(rb, 2_I1B)
		I_cntn1(k) = IFEM1(j)*V1+IFEM1(j+1)*V2
		!dI_FEM(k) =  (dI(j)*V1)**2 + (dI(j+1)*V2)**2 
	END DO !Npoints
		
	I_cntn1(1)=IFEM1(1)
	I_cntn1(Npoints)=IFEM1(Nj)
	!dI_FEM(1)=dabs(dI(1))
	!dI_FEM(Npoints)=dabs(dI(Nj))
		
	IFEM_1=IFEM_1+I_cntn1(1:Npoints)
	!sigmaI=sigmaI+dI_FEM(1:Npoints)
	
ENDIF
		DEALLOCATE(dell, l, delr, rv, mu, nglobl, nglobr, Q_ls)	
	ENDDO	!iset
	
	
IF (SIM) THEN

	ALLOCATE(mu(1:Npoints),QI(1:Npoints))
	mu(1)=1.d0
	DO i=2, Npoints-1
		mu(i)=dsqrt(1.d0-r_cntn(i))*dsqrt(1.d0+r_cntn(i))
	ENDDO
	mu(Npoints)=0.d0
	CALL lstsctrng(mu,Q_ls)
	!CALL lstsctrng(r_cntn,Q_ls)
	
	I_lambda1=I_lambda1/dfloat(Nsetsf)
	!dI_lambda1=dsqrt(dI_lambda1)/dfloat(Nsetsf)
	min_dsprn=1.d14
	DO i=1, Nlambda+1
		QI=0.d0
		DO j=1, Npoints
			QI(j)=dot_product(Q_ls(j,:),I_lambda1(i,:))
		ENDDO
		dsprn(i)=dot_product(I_lambda1(i,:),QI)	
		IF (dsprn(i)<min_dsprn)	min_dsprn=dsprn(i)
		CALL Chi2_SIM_FEM(l_all, A_all, dA_all,r_cntn,I_lambda1(i,1:Npoints),DoF, Chi2(i), wR2)
	ENDDO
	
	!err1=0.d0
	!min_dChi2=1.d12
	!max_dChi2=0.d0
   ! DO i=1, Nlambda+1
		!DO k=2, Npoints-1
			!r=r_cntn(k)
			!IF (I_input1(k)>0) err1(i)=err1(i)+dabs(1.d0-I_lambda1(i,k)/I_input1(k))
		!ENDDO
		!IF (Npoints>2) THEN
			!err1(i)=err1(i)/(Npoints-2)
		!ENDIF	
   ! ENDDO		
		
!-------------- below is to find the first local minimum/maximum of dsprn -------------!
	
!	lmbd_opt1=0.d0
!	min_dsprn=1.d14
		
!	DO i=1, Nlambda
!		var_dsprn(i)=1.d0-dsprn(i)/dsprn(i+1)
!	ENDDO
	
!	DO i=1, Nlambda+1	
		!IF (var_dsprn(i)>0.d0.and.var_dsprn(i+1)<0.d0.and.dabs(var_dsprn(i))>1.d-8) THEN	!local maximum
		!IF (var_dsprn(i)<0.d0.and.var_dsprn(i+1)>0.d0.and.dabs(var_dsprn(i))>1.d-8) THEN	!local minimum
			!pwr=(i1-1)*dlmbd
			!lambda=1.d1**pwr*lambda1
			!lmbd_opt1=lambda
			!IFEM_1=I_lambda1(i,1:Npoints)
			!i_opt=i
			!EXIT
		!ELSE 
!		IF (dsprn(i)<min_dsprn.and.dabs(var_dsprn(i))<1.d-8.and.Chi2(i)<1.5) THEN	!	It could be local min or max but usually 
												!	the first local extermum is a minimum so we don't worry about this
!			IFEM_1=I_lambda1(i,1:Npoints)
!			i_opt=i
!			pwr=(i1-1)*dlmbd
!			lambda=1.d1**pwr*lambda1
!			lmbd_opt1=lambda
!			min_dsprn=dsprn(i)
					
!		ENDIF 
!	ENDDO
!	IF (i<Nlambda) print*,'λ No.',i,' out of ',Nlambda+1
!	IF (lmbd_opt1>eps) THEN
!		print*,"FEM_regularisation Sub. lambda of &
!		local extermum dispersion :",lmbd_opt1,' Chi^2:',Chi2(i_opt)
!		print*,'	Fractional error:',err1(i_opt)
!		print*
!	ENDIF
!------------------ below is to find optimised solution with dispersion<dsprn_r*min_dsprn --------------------!	
IF (lmbd_opt1<eps) THEN	
	min_dChi2=1.d14
	i_opt=1
	DO i=1, Nlambda+1	
		pwr=(i-1)*dlmbd
		lambda=1.d1**pwr*lambda1
		dchi2=dabs(1.d0-Chi2(i))
		!rChi2=nint(dChi2/min_dChi2)
		rDsprn=dsprn(i)/min_dsprn
		!if (talk) write(*,'(a,3(SE12.4),a,SE12.4)')'             ', lambda, Chi2(i),&
		 										!dsprn(i),' Dsprn/min_Dsprn: ',rDsprn
		!IF (dsprn(i)<min_dsprn.or.dChi2<min_dchi2) THEN
		IF ((dChi2<=min_dchi2).and.rDsprn<=dsprn_r) THEN
			IF (dChi2<min_dchi2) min_dChi2=dChi2
			!IF (dsprn(i)<min_dsprn) min_dsprn=dsprn(i)
			!if (talk) print*,lambda,dChi2/min_dChi2, dsprn(i)/min_dsprn,'dChi2*dsprn:',dsprn(i)*dChi2
			!if (talk) print*, lambda, rChi2, rDsprn,'dChi2*dsprn:',dsprn(i)*dChi2
			
			!IF (talk) THEN
				!print*
				!write(*,'(a,3(SE12.4),a,SE12.4)')' Best so far:', lambda, Chi2(i), dsprn(i),
				!												' Dsprn/min_Dsprn: ',rDsprn
				!print*
			!ENDIF	
			!print*
			IFEM_1=I_lambda1(i,1:Npoints)
			lsq1=err1(i)
			lmbd_opt1=lambda
			i_opt=i	
			if (talk) print*,'optimum so far:', i_opt, lambda,rDsprn,dchi2 		
		ENDIF	
	ENDDO
	
	IF (talk) THEN
		print*,"FEM_regularisation Sub. lambda of &
		& global minimum with min |Chi^2-1|:",lmbd_opt1,' Chi^2:',Chi2(i_opt)
		print*,' dispersion:',dsprn(i_opt),' min_dsprn:',min_dsprn
		print*
	ENDIF	
ENDIF	
	

	!close(201)
ELSE
	IFEM_1=IFEM_1/dfloat(Nsetsf)
	!sigmaI=dsqrt(sigmaI)/dfloat(Nsetsf)
	lmbd_opt1=lambda0
ENDIF	
	
END SUBROUTINE FEM_regularisation
!
!=======================!	
SUBROUTINE OUTPUT(FNAM, OUT_array, HEADER)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN)	:: FNAM, HEADER
	REAL(DP), DIMENSION(:,:), INTENT(IN)	:: OUT_array
	INTEGER	:: i, j, dim1, dim2
	CHARACTER(LEN=50)	::  FMTout 
	
	open(300, file=FNAM)

	dim1=size(OUT_array(:,1))
	dim2=size(OUT_array(1,:))
		
	WRITE(FMTout,'(a,I4,a)') '(',dim1,'(SE15.8,x))'
	
	WRITE(300,*) HEADER
	DO i=1, dim2
		WRITE(300,FMTout) (OUT_array(j,i), j=1,dim1)
	ENDDO
	
	CLOSE(300)						   
PRINT*										
PRINT*,'	OUTPUT SUB.	Results are written in file:'
PRINT*,'	           	',FNAM
PRINT*

	
END SUBROUTINE OUTPUT
!-----------------------!
SUBROUTINE Chi2_FEM(F_A,r,IFEM,lmax,DOF, Chi2,wR2)	!ic and nLc should be known
IMPLICIT NONE
	CHARACTER(LEN=*)	::	F_A
	REAL(DP), DIMENSION(:), INTENT(IN)	:: r, IFEM
	REAL(DP), INTENT(IN)	::	lmax
	INTEGER(I2B), INTENT(OUT)	::	DOF
	REAL(DP), INTENT(OUT)	::	Chi2,wR2
	REAL(DP), DIMENSION(:), ALLOCATABLE	:: t_all, l_all, mti, smi
	REAL(DP), DIMENSION(:), ALLOCATABLE :: A_FEM, A_all, dA_all
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: OUT_array
	CHARACTER (LEN=200)	:: FNAMC, HEADER
	REAL(DP):: dAi
	INTEGER(I2B)	::	NData
	INTEGER	::	i, n

	
	WRITE(FNAMC,'(a,a,a)')trim(DATA_PATH),trim(event(ic)%event_Fldr),trim(event &
															  (ic)%tel_file(nLC))	
	CALL t_to_l_DATA(FNAMC, lmax, t_all, l_all, mti, smi)
	NData=size(l_all)
	ALLOCATE(A_FEM(1:NData))
	
	n=size(r)
	ALLOCATE( INodes(1:n), rNodes(1:n) )
	INodes=IFEM
	rNodes=r
	A_FEM=FLUX_FEM(rmin, l_all, NData)
	!PRINT*,'SUBROUTINE Chi2_FEM: Data imported from file: '
	!PRINT*,FNAMC,' No. of data points:',NData	
	ALLOCATE( dA_all(1:NData), A_all(1:NData) )
	
	A_all=( 10**(0.4d0*( Event(ic)%mbl-mti )) - 1 )/Event(ic)%blndng+1	
	dA_all=smi*( A_all-1.d0 + 1/Event(ic)%blndng )*ctm!0.4d0*dlog(1.d1)
	
	DOF=NData!-n
	!PRINT*,'Chi2_FEM SUBROUTINE: No. Data (l<',lmax,') =',n,' DOF=',DOF	
	wR2=0.d0
	DO i=1, NData
		dAi=( A_all(i)-A_FEM(i) )/dA_all(i)
		wR2=wR2+dAi**2	
	ENDDO
	
	IF (DOF>1) THEN			
		Chi2=wR2/dfloat(DOF)
		!PRINT*, '	wR2:',wR2,'	DoF:',DoF,'	Chi2:',Chi2
	ELSE
		PRINT*,'There is not enough data points to calculate Chi^2, We suggest increasing lmax.'
		Chi2=-1.d0	   			   			
	ENDIF
	
	HEADER='#	  t	  	l	  	A	  	A_{FEM}		dA'	
	ALLOCATE(OUT_array(1:5,1:NData))
	OUT_array(1,:)=t_all
	OUT_array(2,:)=l_all
	OUT_array(3,:)=A_all
	OUT_array(4,:)=A_FEM
	OUT_array(5,:)=dA_all
	CALL OUTPUT(F_A, OUT_array, HEADER)
	DEALLOCATE(INodes, rNodes)
END SUBROUTINE Chi2_FEM
!-----------------------!
SUBROUTINE Chi2_Nsqrt(lmax,DOF, Chi2,wR2)	!ic and nLc should be known
IMPLICIT NONE
	REAL(DP), INTENT(IN)	::	lmax
	INTEGER, INTENT(OUT)	::	DOF
	REAL(DP), INTENT(OUT)	::	Chi2,wR2
	REAL(DP), DIMENSION(:), ALLOCATABLE	:: t_all, l_all, mti, smi
	REAL(DP), DIMENSION(:), ALLOCATABLE :: A_Nsqrt, A_all, dA_all
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: OUT_array
	CHARACTER (LEN=200)	:: FNAMC, HEADER
	REAL(DP):: dAi
	INTEGER(I2B)	::	NData
	INTEGER	::	i, n

	
	WRITE(FNAMC,'(a,a,a)')trim(DATA_PATH),trim(event(ic)%event_Fldr),trim(event &
															  (ic)%tel_file(nLC))	
	CALL t_to_l_DATA(FNAMC, lmax, t_all, l_all, mti, smi)
	NData=size(l_all)
	ALLOCATE(A_Nsqrt(1:NData),dA_all(1:NData), A_all(1:NData))
	
	A_Nsqrt=mFLUX(rmin, l_all, NData)
	
	A_all=( 10**(0.4d0*( Event(ic)%mbl-mti )) - 1 )/Event(ic)%blndng+1	
	dA_all=smi*( A_all-1.d0 + 1/Event(ic)%blndng )*ctm!0.4d0*dlog(1.d1)
	
	DOF=NData!-n
	!PRINT*,'Chi2_FEM SUBROUTINE: No. Data (l<',lmax,') =',n,' DOF=',DOF	
	wR2=0.d0
	DO i=1, NData
		dAi=( A_all(i)-A_Nsqrt(i) )/dA_all(i)
		wR2=wR2+dAi**2	
	ENDDO
	
	IF (DOF>1) THEN			
		Chi2=wR2/dfloat(DOF)
		!PRINT*, '	wR2:',wR2,'	DoF:',DoF,'	Chi2:',Chi2
	ELSE
		PRINT*,'There is not enough data points to calculate Chi^2, We suggest increasing lmax.'
		Chi2=-1.d0	   			   			
	ENDIF
	
END SUBROUTINE Chi2_Nsqrt

!-----------------------!
SUBROUTINE Chi2_SIM_FEM(l_SIM, A_SIM, dA_SIM,r,IFEM,DOF, Chi2,wR2)	
IMPLICIT NONE

	REAL(DP), DIMENSION(:), INTENT(IN)	:: l_SIM, A_SIM, dA_SIM, r, IFEM
	INTEGER(I2B), INTENT(OUT)	::	DOF
	REAL(DP), INTENT(OUT)	::	Chi2,wR2
	REAL(DP), DIMENSION(:), ALLOCATABLE	:: t_all, l_all, mti, smi
	REAL(DP), DIMENSION(:), ALLOCATABLE :: A_FEM, A_all, dA_all
	REAL(DP), DIMENSION(:,:), ALLOCATABLE :: OUT_array
	CHARACTER (LEN=200)	:: FNAMC, HEADER
	REAL(DP):: dAi
	INTEGER(I2B)	::	NData
	INTEGER	::	i, n

	n=size(r)
	NData=size(l_SIM)
	ALLOCATE( INodes(1:n), rNodes(1:n) )
	INodes=IFEM
	rNodes=r
	ALLOCATE(A_FEM(1:NData))
	A_FEM=FLUX_FEM(rmin, l_SIM, NData)
	
	DOF=NData!-n
	!PRINT*,'Chi2_FEM SUBROUTINE: No. Data (l<',lmax,') =',n,' DOF=',DOF	
	wR2=0.d0
	DO i=1, NData
		!print*,i,A_SIM(i),A_FEM(i) 
		dAi=( A_SIM(i)-A_FEM(i) )/dA_SIM(i)
		wR2=wR2+dAi**2	
	ENDDO
	
	IF (DOF>1) THEN			
		Chi2=wR2/dfloat(DOF)
		!PRINT*, '	wR2:',wR2,'	DoF:',DoF,'	Chi2:',Chi2
	ELSE
		PRINT*,'There is not enough data points to calculate Chi^2, We suggest increasing lmax.'
		Chi2=-1.d0	   			   			
	ENDIF
	DEALLOCATE(INodes, rNodes)
END SUBROUTINE Chi2_SIM_FEM
!-----------------------!
SUBROUTINE compareIs_and_sigmaI(F_Ir)
IMPLICIT NONE
	CHARACTER(LEN=200), INTENT(OUT)	::	F_Ir
	LOGICAL	::	FLAG
	REAL(DP), ALLOCATABLE, DIMENSION(:) ::  I_data, r_data, sigma_I, A_Data, dA_Data, mean_I2, &
											std_I, ri, wi, t_DATA, l_DATA, mti, smi, mean_I,   &										 
											I_SIM, A_SIM, dA_SIM, I_input1, I_step, A_step,    &
											dA_step, std2_I, r_optm, I_optm, I_input2,I_SIM2, I_step2
											
	REAL(DP)	::  threshold, CNrml0, CNrml1, CNrml2, Chi2_IFEM1, wR2_1, lambda, min_a12, &
					lmbd, lmbd2, mi, dmi, dlambda, MedR
	CHARACTER(len=20)	::	part			      														     
	CHARACTER(len=200)	::  FMT, F_A, FNAMC, F_Istep
	CHARACTER(LEN=600)	::	HEADER
	INTEGER			:: i, iran, Nran, j, j1, t1, t2, dt_run, Nbin1, Nbin2, Nsubs, used, &
						MinNbox, MaxNbox, number
					    
	INTEGER(I2B)	::	NData, d1, d2, ng, DOF	
	INTEGER(I1B), parameter	::  Nmin=1_I1B			
	INTEGER(I1B)	::  nI1, nI2, nI3
	REAL(DP), ALLOCATABLE, DIMENSION(:,:)	:: OUT_array
	INTEGER, ALLOCATABLE, DIMENSION(:,:)	:: setNos, setNos_opt
	INTEGER, ALLOCATABLE, DIMENSION(:)	::	index
	LOGICAL	::	success
	
	t1=time()
 
	Nran=10!20!50
	SIM=.FALSE.
	Verbal=.FALSE.
	
	ng=2_I2B	
	ALLOCATE(ri(1:ng), wi(1:ng))
	CALL gauleg(-1.d0, 1.d0, ri, wi, ng)
			
	i=0
	DO WHILE (i<9)
		i=i+1 
		!print*,good_result(1,i),ic,'.and.',good_result(2,i),nLc
		IF (good_result(1,i)==ic.and.good_result(2,i)==nLc) EXIT
	ENDDO
	IF (i>8) THEN
		print*, 'This ic and nLc are not listed in good_result. We take δ threshold to be 0.4'
	ENDIF
	number=i	
	PRINT*![91m red		[92m green		[93m yellow		[94m blue		[95m pink		[96m light blue	
	
	PRINT*,number,''//achar(27)//'[95m Optimize_and_sigmaI SUB. ic,nLC:',ic,nLC,&
			' minimum impact parameter:'//achar(27)//'[0m',Event(ic)%p 
	PRINT*,'	c1, c2:',c1,c2
	    
	lenfld=len(trim(event(ic)%event_Fldr))
	IF (i<9) THEN
		threshold=delta_trsh(i)	! We started from 0.5, then we reduced it in cases it didn't work
		Nbin1=Nbins(i)
		lambda=lmbd_choice(i)
		dsprn_r=dsprn_bnd(i)
	ELSE
		threshold=0.4d0
		Nbin1=100
		lambda=1.d-2
		dsprn_r=1.1d0
	ENDIF

	IF (ic/=10)	THEN	  	 	
		WRITE(FNAMC,'(a,a,a)')trim(DATA_PATH),trim(event(ic)%event_Fldr), &
											trim(event(ic)%tel_file(nLC))														  	   
		CALL t_to_l_DATA(FNAMC, 1.d0, t_DATA, l_DATA, mti, smi)	
		Ndata=size(l_DATA)
		IF (ic==5.or.ic==8.or.ic==7) THEN
		!IF (ic==7) THEN		
		 	used=0
		 	j=0
		 	DO	WHILE (j<10)	
				CALL data_groups_mu( Nbin1, threshold ,l_DATA, setNos, used, min_a12)
				IF (used==NData) EXIT
				DEALLOCATE(setNos)					
				j=j+1
			ENDDO
			IF (used<Ndata) THEN
				print*,' Not all data used, just ',used,' out of ',Ndata
				STOP
			ENDIF	
		ELSE
		 	used=0
		 	j=0
		 	DO	WHILE (j<10)	
				CALL data_groups( Nbin1, threshold ,l_DATA, setNos, used, min_a12)
				IF (used==NData) EXIT
				DEALLOCATE(setNos)
				j=j+1
			ENDDO
			IF (used<Ndata) THEN
				print*,' Not all data used, just ',used,' out of ',Ndata
				STOP
			ENDIF	
		ENDIF
	ELSE	!ic=10 is a simulated light curve with small p (5*10^-8) and data concentration at limb	
		Nbin1=101
		Ndata=100
		ALLOCATE(l_DATA(1:Ndata))
		dmi=(1.d0-Event(ic)%p2)/(Ndata+1)
		l_Data(1)=Event(ic)%p
		DO i=1, Ndata-1
			mi=i*dmi
			l_Data(Ndata-i+1)=dsqrt(1.d0-mi**2)
		ENDDO
			CALL data_groups_mu( Nbin1, threshold ,l_DATA, setNos, used, min_a12)		
	ENDIF
	
	Nbin2=size(setNos(1,:))
	Nsubs=size(setNos(:,1))
	d2=Ndata+2
print*,'No. of data, No. of profile sample points :', NData, d2
print*,'No. Histogram bins and full bins :',Nbin1, Nbin2
print*,'No. subsets :', Nsubs
IF (used<Ndata) THEN
	print*,'!!!!!!!!!!!!!!  Some of the data points were not used in data subsets. !!!!!!!!!'
	STOP
ENDIF	
	ALLOCATE( OUT_array(1:12,1:d2 ) , I_input1(1:d2 ), I_input2(1:d2 ), r_data(1:d2), std_I(1:d2) )
	ALLOCATE(mean_I(1:d2), mean_I2(1:d2), sigma_I(1:d2) )
	OUT_array=0.d0
		
	r_data(1)=0.d0
	r_data(2:Ndata+1)=l_data(1:Ndata)	
	r_data(d2)=1.d0
	
	!------------ I_input will be used for sigma_I estimation and calculating median(R) --------!
	!print*, 'r, I_tanh'		
	DO i=1, d2		
	 	I_input2(i)=tanhLD(r_data(i))
	 	I_input1(i)=NSQrtLD(r_data(i), c1, c2)
	 	!print*, r_data(i), I_input2(i)
    ENDDO	
	!-----------------------------------------------------------------!
	!!!!!!!!!!!!! so the FEM_regularisation_step SUB. will be run with the specified lambda and Nbin
	talk=.true.
	IF (ic/=10) THEN
		SIM=.TRUE.
		print*,'Simulation of square root profile...'	
		j1=-time() ! iran>1 is simulation with input profile of a square root law
		iran=nint(ran2(j1)*1000)+1													  														  										
		CALL m_to_A_DATA(iran, l_DATA, mti, smi, A_sim, dA_sim)
		!CALL FEM_regularisation_step( lambda, setNos, l_DATA, A_sim, dA_sim, A_step, dA_step, &
									  ! I_input1, I_input2, r_data, I_sim1, I_step1, lmbd1, lmbd2 )
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_sim, dA_sim, &
										I_input1, r_data, I_sim, lmbd )										
	  	
		!SIM=.FALSE.		
		!lambda=lmbd2	
		iran=-2	! iran=-2 is amplification from an tanhLD profile 
		!IF (Istep) iran=-1	! iran=-1 is amplification from an Intensity profile with a step function at limb 
		print*,'Simulation of tanh profile...'
		CALL m_to_A_DATA(iran, l_DATA, mti, smi, A_step, dA_step)
		!CALL FEM_regularisation_step( lambda, setNos, l_DATA, A_step, dA_step, &
										!I_input, r_data, I_step, std2_I, lmbd )
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_step, dA_step, &
										I_input2, r_data, I_step, lmbd )										
		iran=0 ! this is real data
		print*,'Data...'
		CALL m_to_A_DATA(iran, l_DATA, mti, smi, A_Data, dA_Data)
		!CALL FEM_regularisation_step( lambda, setNos, l_DATA, A_Data, dA_Data, A_step, dA_step, &
										!I_input1, I_input2, r_data, I_Data, I_step1, lmbd1, lmbd2 )
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_data, dA_data, &
										I_input1, r_data, I_data, lmbd2 )										
		print*,'A range:',A_data(Ndata),A_data(1)
		print*,'I range from data:',I_data(1),I_data(d2)
		print*,'l range:',r_data(2),r_data(d2-1)
	ELSE
		SIM=.TRUE.
		iran=-2	! iran=-2 is amplification from an tanhLD profile 
		CALL SIM_A_DATA(iran, l_DATA, A_Data, dA_Data)
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_data, dA_data, &
										I_input2, r_data, I_data, lmbd2 )
	  								
		print*,'A range from a tanh profile:',A_data(Ndata),A_data(1)
		print*,'I range:',I_data(1),I_data(d2)
		print*,'l range:',r_data(2),r_data(d2-1)
																				
		!SIM=.FALSE.
		!lambda=lmbd2
		j1=-time() ! iran>1 is simulation with input profile of a square root law
		iran=nint(ran2(j1)*1000)+1													  														  										
		CALL SIM_A_DATA(iran, l_DATA, A_sim, dA_sim)
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_sim, dA_sim, &
										I_input1, r_data, I_sim, lmbd )										
		
		iran=-1 ! this is step function
		CALL SIM_A_DATA(iran, l_DATA, A_step, dA_step)
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_step, dA_step, &
										I_input1, r_data, I_step, lmbd )										
										
	ENDIF									
	talk=.false.		
	OUT_array(1,1:d2)=r_data(1:d2)  
	OUT_array(2,1:d2)=I_data(1:d2)
	OUT_array(3,1:d2)=I_sim(1:d2) 
	OUT_array(4,1:d2)=I_step(1:d2)
    
	OUT_array(6,1:NData)=l_Data(1:NData)
	OUT_array(7,1:NData)=A_Data(1:NData)
	OUT_array(8,1:NData)=dA_Data(1:NData)		    	    	    		
	OUT_array(9,1:NData)=A_sim(1:NData)
	OUT_array(10,1:NData)=dA_sim(1:NData)
	OUT_array(11,1:NData)=A_step(1:NData)
	OUT_array(12,1:NData)=dA_step(1:NData)
	    	    	    
	IF (Nbin2>0) THEN
	    nI2=INT(log10(float(Nbin2)))+1
	    print*,nI2
	ELSE
	    nI2=1
	ENDIF
	    	
	IF (lmbd2>=1.d0) THEN
	    nI1=INT(log10(lmbd2))+1
	    write(fmt,'(a,I1,a)')'(a,I',nI1,',a)'
	    write(part,fmt)'_lambda',NINT(lmbd2),'_'	
	ELSE IF (lmbd2>=1.d-2) THEN
	    nI1=INT(log10(lmbd2*1.d3))+1
	    write(fmt,'(a,I1,a)')'(a,I',nI1,',a)'
	    write(part,fmt)'_lambda',NINT(lmbd2*1.d3),'d-3_'
	ELSE
	    nI1=INT(log10(lmbd2*1.d5))+1
	    !print*, lmbd2, nI1
	    write(fmt,'(a,I1,a)')'(a,I',nI1,',a)'
	    write(part,fmt)'_lambda',NINT(lmbd2*1.d5),'d-5_'
	ENDIF	

	nI1=INT(log10(float(Nbin1)))+1
	nI3=INT(log10(float(NData)))+1
	    
	WRITE(FMT,'(a,I1,a,I1,a,I1,a)')'(a,I',nI1,',a,I',nI2,',a,I',nI3,'a,a,a,a)'
 
!final:	WRITE(F_Ir, FMT) 'Intensity-defence-tanh-min_dsprn_dChi2-hist',Nbin1,'_size',Nbin2,'_No-data',NData,trim(part) &
!    					, trim( event(ic)%event_Fldr(1:lenfld-5) ),trim( event(ic)%tel_file(nLC) )
	WRITE(F_Ir, '(a,I2,a,a,a)') 'Intensity-10alpha',nint(10*dsprn_r),'event' &
    					, trim( event(ic)%event_Fldr(1:lenfld-5) ),trim( event(ic)%tel_file(nLC) )	    								  
	    								  
!final:	WRITE(F_A, FMT) 'LC-A_FEM-stdI-tanh-min_dsprn_dChi2-hist',Nbin1,'_size',Nbin2,'_No-data',NData,trim(part) &
!	    				, trim( event(ic)%event_Fldr(1:lenfld-5) ),trim( event(ic)%tel_file(nLC) )
	WRITE(F_A, '(a,I2,a,a,a)') 'LC-A_FEM-10alpha',nint(10*dsprn_r),'event' &
    					, trim( event(ic)%event_Fldr(1:lenfld-5) ),trim( event(ic)%tel_file(nLC) )	    								  

	WRITE(F_Istep,'(a,a,a)')'I_tanh-defence-', trim( event(ic)%event_Fldr(1:lenfld-5) ),trim( event(ic)%tel_file(nLC) ) 
	OPEN (202,file=F_Istep)	
	DO i=1, 1000
		write(202,'(3(SE12.5,x))')i*1.d-3,tanhLD(i*1.d-3),NSQrtLD(i*1.d-3, c1, c2)
	ENDDO 
	CLOSE(202)   				
print*,F_Ir

	CNrml1=Nrml_factor(ri, wi, r_data, I_data)														
	PRINT*,' Optimize_and_sigmaI SUB. Data: 2π∫I_data.rdr=',CNrml1
	CNrml2=Nrml_factor(ri, wi, r_data, I_step)
	CNrml0=Nrml_factor(ri, wi, r_data, I_input2)
	PRINT*,' Optimize_and_sigmaI SUB. Simulation: 2π∫I_(FEM,tanh).rdr=',CNrml2,' 2π∫I_tanh.rdr=',CNrml0
	!CNrml0=Nrml_factor(ri, wi, r_data, I_input1)
	!PRINT*,' Optimize_and_sigmaI SUB. Simulation:  2π∫I_SQrtLD.rdr=',CNrml0

IF (ic/=10) THEN
	CALL Chi2_FEM(F_A,r_data,I_data,1.d2,DOF,Chi2_IFEM1,wR2_1)	   
	PRINT*,' compareIs_and_sigmaI SUB. 	DOF1:',DOF,' Weighted Residuals^2=',wR2_1		
	IF (wR2_1>2*DoF) print*,'	Alert!!!! Chi^2 is too high, changing λ range may help. !!!!!!'
ELSE
	CALL Chi2_SIM_FEM(l_Data, A_SIM, dA_SIM,r_data,I_sim,DOF, Chi2_IFEM1,wR2_1)
	IF (wR2_1>2*DoF) print*,'	Alert!!!! Chi^2 is too high, changing λ range may help. !!!!!!'
ENDIF	   

	!dlambda=lambda/2.d1
IF (ic/=10) THEN													  														  
	SIM=.TRUE.	
	sigma_I=0.d0		
	mean_I=0.d0
	mean_I2=0.d0														  
	DO i=1, Nran
		
		j1=-time()
		iran=nint(ran2(j1)*100+i)
		CALL m_to_A_DATA(iran, l_DATA, mti, smi, A_sim, dA_sim)
		IF (ic==5.or.ic==8.or.ic==7) THEN
			CALL data_groups_mu( Nbin1, threshold ,l_DATA, setNos, used, min_a12)
		ELSE
			CALL data_groups( Nbin1, threshold ,l_DATA, setNos, used, min_a12)
		ENDIF
		CALL FEM_regularisation( lambda, setNos, l_DATA, A_sim, dA_sim, &
										I_input1, r_data, I_sim, lmbd )											    
	    mean_I=mean_I+I_sim
	    mean_I2=mean_I2+I_sim**2
		print*,i,iran,I_sim(1),' λ=',lmbd
	    
	    DEALLOCATE(A_sim, dA_sim, r_data, I_sim)
	ENDDO
	mean_I2=mean_I2/Nran
	mean_I=mean_I/Nran
	sigma_I=dsqrt(mean_I2-mean_I**2)
	ALLOCATE(index(1:d2))
	CALL indexx_dp(sigma_I,index)
	MedR=sigma_I(index(d2/2))
	PRINT*,'		median(sigma_I_FEM)=',MedR
ELSE
	SIM=.TRUE.	
	sigma_I=dabs(I_input2-I_Data)		
	ALLOCATE(index(1:d2))
	CALL indexx_dp(sigma_I,index)
	MedR=sigma_I(index(d2/2))
	PRINT*,'		median(sigma_I_FEM)=',MedR
	print*,'r,		IFEM_sqrt,			I_sqrt			Residual'
	DO i=1, d2
		print*,r_data(i),I_sim(i),I_input1(i),dabs(I_sim(i)-I_input1(i))
	ENDDO
ENDIF
	 
	OUT_array(5,:)=sigma_I
	
	!WRITE(HEADER,'(2(a,SE11.5),a,I4,2(a,SE14.5),a,2(f8.3,x))')'#	1:r		2:I_data		3:I_sqrt	4:I_step	&
		!&5:σ_I(from std(I_sqrt,I_input) ) 6:l	7:A_data		8:dA_data	9:A_sqrt		10:dA_sqrt 	9:A_sqrt	&
		!&	10:dA_sqrt		 ,    2π∫I_data.rdr=',CNrml1,' 2π∫I_tanh.rdr=',CNrml2,'DOF',DOF,' Weighted Residual^2:', &		
		!wR2_1,'   median(σ_I)=',MedR,' c, d=',c1,c2
	WRITE(HEADER,'(a)')'1:r		         2:I_data		3:I_sqrt	  4:I_tanh	        5:σ_I         &
	 	 &6:l	            7:A_data        8:dA_data	    9:A_sqrt		10:dA_sqrt 	   11:A_tanh       12:dA_tanh'		                                                                                                                                                                                                                                                                                                                
	
	CALL OUTPUT(F_Ir, OUT_array, HEADER)
						
	t2=time()
	dt_run=t2-t1
	print*,''//achar(27)//'[95m compareIs_and_sigmaI SUB. Optimisation and estimating sigma_I for this event by ',&
			Nran,' LCs, took ', dt_run,'sec'//achar(27)//'[0m'

END SUBROUTINE compareIs_and_sigmaI
!-----------------------!
FUNCTION Nrml_factor(ri, wi, rFEM, IFEM)	!ri and wi are from gauleg subroutine
IMPLICIT NONE

	REAL(DP), DIMENSION(:), INTENT(IN)	::	rFEM, IFEM
	REAL(DP), DIMENSION(:), INTENT(IN)	::	ri, wi
	REAL(DP)	::	IFEM_ri, Nrml_factor
	INTEGER		::	i, j, nIFEM, iel
	INTEGER(I2B)::	ng
	REAL(DP)	::	sumi, sumj, r_i, r_i1, r_i2, dr_i, rb, &
					V1, V2, IFEM_r
	
	ng=size(wi)!2_I2B
	nIFEM=size(IFEM)
	IF (nIFEM/=size(rFEM)) THEN
		print*,'Error in using Nrml_factor function.'
		STOP
	ENDIF
	sumi=0
	DO i=1, nIFEM-1
		sumj=0
		r_i1=rFEM(i)
		r_i2=rFEM(i+1)
		dr_i=r_i2-r_i1
		DO j=1, ng
			!rb=2.d0*(ri(j)-r_i1)/dr_i-1.d0
			rb=ri(j)
			r_i=dr_i*(rb+1.d0)/2.d0+r_i1
			V1=SHAPE1D1(rb,1_I1B)
			V2=SHAPE1D1(rb, 2_I1B)
			IFEM_r = IFEM(i)*V1+IFEM(i+1)*V2
		 	!sumj=sumj+wi(j)*ri(j)*IFEM_r
		 	sumj=sumj+wi(j)*r_i*IFEM_r*dr_i/2.d0
		 	!print*, rb, r_i, sumj
		ENDDO
		sumi=sumi+sumj
	ENDDO
	
	Nrml_factor=sumi*2.d0*PI_D

END FUNCTION Nrml_factor
!-----------------------!
SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,sig)
IMPLICIT NONE

	REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(DP), INTENT(OUT) :: a,b,siga,sigb,chi2
	REAL(DP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
	
	!Given a set of data points in same-size arrays x and y, fit them to a straight line y = a+bx
	!by minimizing χ2. sig is an optional array of the same length containing the individual
	!standard deviations. If it is present, then a,b are returned with their respective probable
	!uncertainties siga and sigb, the chi-square chi2. If sig is not present the normalization of 
	!chi2 is to unit standard deviation on all points.
	
	INTEGER(I4B) :: ndata
	REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
	REAL(SP), DIMENSION(size(x)), TARGET :: t
	REAL(SP), DIMENSION(:), POINTER :: wt
	
	if (present(sig)) then
		ndata=assert_eq(size(x),size(y),size(sig),'fit')
		wt=>t !Use temporary variable t to store weights.
		wt(:)=1.0_sp/(sig(:)**2)
		ss=sum(wt(:)) !Accumulate sums with weights.
		sx=dot_product(wt,x)
		sy=dot_product(wt,y)
	else
		ndata=assert_eq(size(x),size(y),'fit')
		ss=real(size(x),sp) !Accumulate sums without weights.
		sx=sum(x)
		sy=sum(y)
	end if
	sxoss=sx/ss
	t(:)=x(:)-sxoss
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		b=dot_product(t/sig,y)
	else
		b=dot_product(t,y)
	end if
	st2=dot_product(t,t)
	b=b/st2 !Solve for a, b, σa, and σb.
	a=(sy-sx*b)/ss
	siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
	sigb=sqrt(1.0_sp/st2)
	t(:)=y(:)-a-b*x(:)

	if (present(sig)) then
		t(:)=t(:)/sig(:)
		chi2=dot_product(t,t) !Calculate χ2.
	else
		chi2=dot_product(t,t)
		sigdat=sqrt(chi2/(size(x)-2)) !For unweighted data evaluate typical
									  !sig using chi2, and adjust the
									  !standard deviations.
		siga=siga*sigdat
		sigb=sigb*sigdat
	end if
END SUBROUTINE fit
!-----------------------!
SUBROUTINE FIT_SQrtLD(x, y, a, b, c, siga, sigb, sigc, Chi2, sig)
IMPLICIT NONE
! This subroutine find the best fit of y=a+b.x+c.√x for a, b, c, the ranges are [0,1]
	REAL(DP), DIMENSION(:), INTENT(IN):: x, y, sig
	REAL(DP), INTENT(OUT)	::	a, b, c, siga, sigb, sigc, Chi2
	
	INTEGER	::	DoF, na, nb, nc	
	REAL(DP)	::	ai, bi, ci, da, db, dc, WR2, minWR2	
	REAL(DP), DIMENSION(:), ALLOCATABLE	::	y_par, WR
	CHARACTER(LEN=30)	::	FChi2			
	INTEGER	::	i, n, ia, ib, ic
	
	IF (outornot) THEN
		WRITE(FChi2,'(a,I1,a,I1,a)')'Chi2-SQFIT-Nsqrt-ic',ic,'-Nlc',Nlc,'.txt'
		OPEN(104, file=FChi2)
		WRITE(104,*) '#	  a   ,   b   ,     c,      Chi2      '
	ENDIF

	n=size(x)
	DoF=n-3
	a=0.d0
	b=0.d0
	c=0.d0
	IF (DoF>0) THEN
		ALLOCATE(y_par(1:n), WR(1:n))
		da=2.d-3
		db=2.d-3
		dc=2.d-3
		na=501
		nb=501
		nc=501
		Chi2=1.d12
		DO ia=1, na
			ai=(ia-1)*da
			DO ib=1, nb
				bi=(ib-1)*db
				DO ic=1,nc
					ci=(ic-1)*dc
					y_par(:)=ai+bi*x+ci*dsqrt(x)
					WR(:)=(y_par-y)/sig
					WR2=dot_product(WR,WR)/DoF
					IF (WR2<Chi2) THEN
						a=ai
						b=bi
						c=ci
						Chi2=WR2
					ENDIF
				ENDDO
			ENDDO
		ENDDO		
!---------------- siga derivation ------------------!	
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.ai>=-1)
			ai=ai-da
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=a-ai
		
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.ai<=2)
			ai=ai+da
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=MAX(siga,ai-a)
		
!---------------- sigb derivation ------------------!	
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.bi>=-1)
			bi=bi-db
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigb=b-bi
		
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.bi<=2)
			bi=bi+db
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigb=MAX(siga,bi-b)

!---------------- sigc derivation ------------------!	
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.ci>=-1)
			ci=ci-dc
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigc=c-ci
		
		WR2=0
		ai=a
		bi=b
		ci=c		
		DO WHILE (WR2<1+Chi2.and.ci<=2)
			ci=ci+dc
			y_par(:)=ai+bi*x+ci*dsqrt(x)
			WR(:)=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigc=MAX(sigc,ci-c)
		
	ELSE
		PRINT*,'Not enough data points to fit y=a+b.x+c.√x'	
		
	ENDIF	


END SUBROUTINE FIT_SQrtLD
!-----------------------!
SUBROUTINE FIT_NSQrtLD(x, y, a, b, siga, sigb, Chi2, sig)
IMPLICIT NONE
! This subroutine find the best fit of y=a+b.x+c.√x for a, b, c, the ranges are [0,1]
	REAL(DP), DIMENSION(:), INTENT(IN):: x, y, sig
	REAL(DP), INTENT(OUT)	::	a, b, siga, sigb, Chi2
	
	INTEGER	::	DoF, na, nb
	REAL(DP)	::	a0, ai, b0, bi, da, db, WR2, minWR2	
	REAL(DP), DIMENSION(:), ALLOCATABLE	::	y_par, WR
	CHARACTER(LEN=30)	::	FChi2			
	INTEGER	::	i, n, ia, ib, ic
	
	IF (outornot) THEN
		WRITE(FChi2,'(a,I1,a,I1,a)')'Chi2-NSQFIT-Nsqrt-ic',ic,'-Nlc',Nlc,'.txt'
		OPEN(104, file=FChi2)
		WRITE(104,*) '#	  a   ,   b   ,     Chi2      '
	ENDIF

	n=size(x)
	DoF=n-2
	a=0.d0
	b=0.d0
	IF (DoF>0) THEN
		ALLOCATE(y_par(1:n), WR(1:n))
		a0=-0.55d0	!0.3d0
		da=1.d-3
		b0=-0.25d0
		db=1.d-3	!0.0d0
		na=1650!51
		nb=1550!81
		Chi2=1.d12
		DO ia=1, na
			ai=a0+(ia-1)*da
			DO ib=1, nb
				bi=b0+(ib-1)*db
				DO i=1, n
					y_par(i)=NSQrtmu(x(i), ai, bi)
				ENDDO	
				WR=(y_par-y)/sig
				WR2=dot_product(WR,WR)/DoF
				IF (WR2<Chi2) THEN
					a=ai
					b=bi
					Chi2=WR2
				ENDIF
			ENDDO
		ENDDO		
!---------------- siga derivation ------------------!	
		WR2=0
		ai=a
		bi=b		
		DO WHILE (WR2<1+Chi2.and.ai>=-1.d0)
			ai=ai-da
			DO i=1, n
				y_par(i)=NSQrtmu(x(i), ai, bi)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=a-ai
		
		WR2=0
		ai=a
		bi=b		
		DO WHILE (WR2<1+Chi2.and.ai<=2.d0)
			ai=ai+da
			DO i=1, n
				y_par(i)=NSQrtmu(x(i), ai, bi)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=MAX(siga,ai-a)
		
!---------------- sigb derivation ------------------!	
		WR2=0
		ai=a
		bi=b		
		DO WHILE (WR2<1+Chi2.and.bi>=-1.d0)
			bi=bi-db
			DO i=1, n
				y_par(i)=NSQrtmu(x(i), ai, bi)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigb=b-bi
		
		WR2=0
		ai=a
		bi=b		
		DO WHILE (WR2<1+Chi2.and.bi<=2.d0)
			bi=bi+db
			DO i=1, n
				y_par(i)=NSQrtmu(x(i), ai, bi)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		sigb=MAX(sigb,bi-b)
		
	ELSE
		PRINT*,'Not enough data points to fit y=a+b.x+c.√x'	
		
	ENDIF	


END SUBROUTINE FIT_NSQrtLD
!-----------------------!
SUBROUTINE FIT_NLLD(x, y, a, siga, Chi2, sig)
IMPLICIT NONE
! This subroutine find the best fit of y=a+b.x+c.√x for a, b, c, the ranges are [0,1]
	REAL(DP), DIMENSION(:), INTENT(IN):: x, y, sig
	REAL(DP), INTENT(OUT)	::	a, siga, Chi2
	
	INTEGER	::	DoF, na
	REAL(DP)	::	a0, ai, da, WR2, minWR2	
	REAL(DP), DIMENSION(:), ALLOCATABLE	::	y_par, WR
	CHARACTER(LEN=30)	::	FChi2			
	INTEGER	::	i, n, ia
	
	IF (outornot) THEN
		WRITE(FChi2,'(a,I1,a,I1,a)')'Chi2-NSQFIT-Nsqrt-ic',ic,'-Nlc',Nlc,'.txt'
		OPEN(104, file=FChi2)
		WRITE(104,*) '#	  a   ,   b   ,     Chi2      '
	ENDIF

	n=size(x)
	DoF=n-1
	a=0.d0
	IF (DoF>0) THEN
		ALLOCATE(y_par(1:n), WR(1:n))
		da=1.d-3	
		na=801
		a0=0.1!0.3
		Chi2=1.d12
		DO ia=1, na
			ai=a0+(ia-1)*da
			DO i=1, n
				y_par(i)=NLLD(x(i), ai, 0.d0)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
			IF (WR2<Chi2) THEN
				a=ai
				Chi2=WR2
			ENDIF
		ENDDO		
!---------------- siga derivation ------------------!	
		WR2=0
		ai=a		
		DO WHILE (WR2<(1.d0+Chi2).and.ai>=-1.d0)
			ai=ai-da
			DO i=1, n
				y_par(i)=NLLD(x(i), ai, 0.d0)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=a-ai
		
		WR2=0
		ai=a		
		DO WHILE (WR2<(1.d0+Chi2).and.ai<=2)
			ai=ai+da
			DO i=1, n
				y_par(i)=NLLD(x(i), ai, 0.d0)
			ENDDO	
			WR=(y_par-y)/sig
			WR2=dot_product(WR,WR)/DoF
		ENDDO
		siga=MAX(siga,ai-a)
		
		
	ELSE
		PRINT*,'Not enough data points to fit y=a+b.x+c.√x'	
		
	ENDIF	


END SUBROUTINE FIT_NLLD

!-----------------------!
FUNCTION Ordr_Mgntd(realNo)
IMPLICIT NONE
	REAL(DP), INTENT(IN)	::	RealNo	
	INTEGER	:: Ordr_Mgntd
	INTEGER	::	intlog
	
	intlog=NINT(dlog10( dabs(realNo) ))
	
END FUNCTION Ordr_Mgntd
!-----------------------!
						
END MODULE Limb_fem2
!***********************!

PROGRAM Limb
USE Limb_fem2
IMPLICIT NONE

	INTEGER	:: t1, t2, dt_run, len_file				
	REAL(DP), ALLOCATABLE, DIMENSION(:) :: rFEM, IFEM1, IFEM2, IFEM3, sigma_I, mu, I_input
	CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE	::	F_Irs
	CHARACTER(LEN=200)::F_Ir, FileName, F_test
	INTEGER, DIMENSION(:), ALLOCATABLE	:: ici, nLci
	INTEGER	::  i, imax, j, n, lenF, ierror, nfiles, DOF									     
	REAL(DP)	:: lmin, tmax, a, b, siga, sigb, u_LLD, du, u_LLD_sim, du_sim, &
					u1, u2, du1, du2,u1_sim, u2_sim, du1_sim, du2_sim, c, sigc, Cnrml, Flux, std,&
					minMbl, dMbl_n, dMbl_p, minb, db_n, db_p, minChi2, Chi2_sim, max_sm, wR2
	REAL(DP)	:: Chi2
	CHARACTER	:: user_input
	LOGICAL	::	BestFit
	
	I1eqI2=.TRUE.
	Nd1=2
	Nd2=2
	t1=time()
	outornot=.FALSE.	
	calmbl=.FALSE.
	BestFit=.FALSE.
	Istep=.False.	
	CALL Events_info
	PRINT*, 'Events information imported from Events.f95'
!-----------------------------------------------------------------------------------------	
	!print *, char(7)
	!PRINT*,' Do you want to calculate blending parameters and base line magnitudes? (Y/N)'
	!READ(*,*) user_input
	!user_input='n'
	!IF (user_input=='Y'.or.user_Input=='y') THEN
		!calmbl=.TRUE.
	!ELSE IF (user_input=='N'.or.user_Input=='n') THEN
		!calmbl=.FALSE.
	!ELSE
		!print*,'Invalid character, we take it as a no.'
		!calmbl=.FALSE.
	!ENDIF
	
	!IF (calmbl) THEN
		!tmax=3 !in terms of tE
		!max_sm=0.1d0
		!write(Filename,'(a,I1,a)')'mBL-b-tmax',NINT(tmax),'tE-lminp-max_sm5d-2.txt'	
		!OPEN (101, file=filename)
		!WRITE(101,*) '# ic, nLc,  Mbl  , dMbl_n , dMbl_p ,    b   ,  db_n  ,   db_p , min χ^2,  DoF   '
		!DO i=1,9
			!ic=i!good_result(1,i)
			!NLc=Event(ic)%nref
			!nLc=Event(ic)%nref_figs
			!IF (Event(ic)%p<eps) THEN
				!lmin=Event(ic)%p!1.5d0
			!ELSE
				!lmin=Event(ic)%p
			!ENDIF		
			
			!CALL MblFit(lmin, tmax, max_sm, minMbl, dMbl_n, dMbl_p, minb, db_n, db_p, minChi2, DoF)
			!WRITE(101,'(2(x,I3),7(x,f8.3),x,I6)') ic, NLc, minMbl, dMbl_n, dMbl_p, minb, db_n, db_p, minChi2, DoF
		!ENDDO
		!CLOSE(101)
	!ENDIF
!-----------------------------------------------------------------------------------------	
	!print *, char(7)
		
	PRINT*,' Do you want to run optimisation on LCs (from Events.f95) ?(Y/N)'
	READ(*,*) user_input
	!user_input='y'
	IF (user_input=='Y'.or.user_Input=='y') THEN
		
		Nfiles=1
		ALLOCATE(F_Irs(1:Nfiles),ici(1:Nfiles),nLci(1:Nfiles))
		BestFit=.FALSE.
		DO i=1, Nfiles
			ic=good_result(1,i)
			nLc=good_result(2,i)
			c2=0.0d0!Event(ic)%u1fit(nLc)
			c1=0.6d0!0.9d0-c2!0.d0!c1/2.d0				
			!c1=Event(ic)%u1fit(nLc)
			!c2=0.d0			
			CALL compareIs_and_sigmaI(F_Ir)
			ici(i)=ic
			nLci(i)=nLc
			F_Irs(i)=F_Ir
			print*,F_Irs(i)
			print*,' Dispersion <',dsprn_bnd(i),"*min_dispersion.  δ threshold=",delta_trsh(i)
		ENDDO
		PRINT*,' Output files:'
		DO i=2, Nfiles
			print*,F_Irs(i)
		ENDDO	
		
	ELSE IF (user_input=='N'.or.user_Input=='n') THEN
		print *, char(7)
		print*,'Do you have the optimised recovered profiles in separate files and want to fit them to square root profiles?(Y/N)'
		READ(*,*) user_input	
		IF (user_input=='Y'.or.user_Input=='y') THEN
		
			Print*,' Please make a file containing file names and ic and nLc (according to Events.f95) in each row '
			print*," And enter the file name:"
			READ(*,*) FILENAME
			
			OPEN(1001, file=FILENAME, action='read')
			n=0
			DO 
				read(1001,*,iostat=ierror) F_test, ic, nLc
				IF (ierror<0) EXIT
				n=n+1
			ENDDO
			Nfiles=n	
			REWIND(1001)
			ALLOCATE(F_Irs(1:Nfiles),ici(1:Nfiles),nLci(1:Nfiles))
			DO i=1,Nfiles
				read(1001,*,iostat=ierror)F_Irs(i), ici(i), nLci(i)
				print*,F_Irs(i), ici(i), nLci(i)
			ENDDO	
			BESTfit=.TRUE.
			print*,'No. of profiles:',Nfiles  		
		ELSE IF (user_input=='N'.or.user_Input=='n') THEN
			print*,' We will update our code please come back later for more. See you!'
		ELSE
			print*,'Invalid character, See you I guess!'
		ENDIF
		
	ELSE
		print*,'Invalid character, End of program.'
	ENDIF
!-----------------------------------------------------------------------------------------	
	 		
	IF (BestFit) THEN
		outornot=.TRUE.
		PRINT*, 'Do you want to fit with a Nsqrt profile or a LLD? (s/l)'
		read(*,*) user_input
		IF (user_input=='S'.or.user_Input=='s') THEN
			Nsqrt_mode=.TRUE.
		ELSE IF (user_input=='L'.or.user_Input=='l') THEN
			Nsqrt_mode=.FALSE.
		ELSE
			print*,'Invalid character! We will fit LLD'
			NSqrt_mode=.FALSE.
		ENDIF

		IF (Nsqrt_mode) THEN
			WRITE(F_test,'(a,I1,a,I1,a)')'Data-fit-sqrt-2.txt'
			OPEN(102, file=F_test)
			WRITE(102,'(a,f7.4,a,f7.4)')'# The input profile parameters in case of simulatopn are: c,d'
			WRITE(102,*)' #ic,nlc, u1   ,    u2     ,   du1   ,   du2   ,  minChi2,&
					&    u1_sim   , u2_sim ,   du1_sim ,  du2_sim , minChi2_sim,DOF,    c    ,    d '
		ELSE
			WRITE(F_test,'(a,I1,a,I1,a)')'Data-fit-linear-2.txt'
			OPEN(102, file=F_test)
			WRITE(102,*)'#ic,nLc, u_LLD   ,    σ_u   ,  minChi2  , u_LLD_sim ,σ_u_sim ,minChi2_sim,DOF,    c    ,    d '
		ENDIF
		
		
		DO	j=1, Nfiles
			
			ic=ici(j)
			nLc=nLci(j)
			c2=Event(ic)%u1fit(nLc)
			c1=0.9d0-c2!0.d0!c1/2.d0				
	
			lenF=len( trim( event(ic)%event_fldr ) )
			OPEN(200,file=F_Irs(j),ACTION='READ')	
!1:r		2:I_FEM		3:I_sim		4:σ_I 		 5:l	6:A_data		7:dA_data	8:A_SIM		9:dA_SIM
			n=0
			READ(200,*)
			DO
	  			READ(200,*,iostat=ierror) a,a,a, a,a,a, a,a,a, a,a,a
	  			IF (ierror<0) EXIT
	  			n=n+1	
			ENDDO
					
			REWIND(200)
			ALLOCATE( rFEM(1:n),IFEM1(1:n),Sigma_I(1:n),IFEM2(1:n),mu(1:n), I_input(1:n), IFEM3(1:n) )
			
			READ(200,*)
			std=0.d0
			imax=0
			DO i=1, n
				READ(200,*,iostat=ierror) rFEM(i),IFEM1(i),IFEM2(i), IFEM3(i),sigma_I(i),a, a,a,a, a,a,a
				IF (imax==0.and.rFEM(i)>0.99) imax=i
				mu(i)=dsqrt(1.d0-rFEM(i)**2)
	  			IF (ierror<0) EXIT	
			ENDDO
			imax=imax-1
			!print*,'imax:',imax,' r_max:',rFEM(imax),'	c,d:',c1,c2
			!sigma_I=sigma_I*std
			!sigma_I=1.d0
			!mu(1)=1.d0
			!mu(n)=0.d0
			!mu(2:n-1)=dsqrt(1.d0-rFEM(2:n-1)**2)
			IF (NSqrt_mode)	THEN
				DoF=imax-2
				IF (DoF>0) THEN
					
					CALL FIT_NSQrtLD(mu(1:imax), IFEM1(1:imax), u1, u2, du1, du2, Chi2, sigma_I(1:imax))

					CALL FIT_NSQrtLD(mu(1:imax), IFEM2(1:imax), u1_sim, u2_sim, du1_sim, du2_sim, Chi2_sim, sigma_I(1:imax))
					!I0=1.d0/(PI_D*(1.d0-u1/3.d0-u2*0.2d0))
					write(102,'(2I3,10(f10.5,x), I4,2(f10.5,x))')ic,nLc, u1, u2, du1, du2, Chi2, u1_sim, u2_sim, du1_sim, &
																du2_sim, Chi2_sim, DOF, c1, c2
					!write(*,'(2I3,10(f10.5,x), I4,2(f10.5,x))')ic,nLc, u1, u2, du1, du2, Chi2, u1_sim, u2_sim, du1_sim, &
																!du2_sim, Chi2_sim, DOF, c1, c2
				c1=u1
				c2=u2
				ENDIF
			ELSE				
				DoF=imax-1
				IF (DoF>0) THEN

				 	CALL FIT_NLLD(rFEM(1:imax), IFEM1(1:imax), u_LLD, du, Chi2, sigma_I(1:imax))
				 	CALL FIT_NLLD(rFEM(1:imax), IFEM2(1:imax), u_LLD_sim, du_sim, Chi2_sim, sigma_I(1:imax))			
					write(102,'(2I3,x,6(f10.5,x), I4,2(f10.5,x))') ic,nlc,u_LLD, du, Chi2,u_LLD_sim, du_sim, Chi2_sim, DoF, c1,c2
					c1=u_LLD
					c2=0.d0
				ENDIF
					
			ENDIF
			CALL Chi2_Nsqrt(1.d2,DOF, Chi2,wR2)
			print*
			print*,'ic,nlc:',ic,nlc
			print*,'	best-fitted parameter:',c1,c2
			Print*,'	Weighted Residual of best-fitted parameter:',wR2,' No. of data points with l<100:',DoF
			print*
			!WRITE(*,'(a,f7.4,a,f7.4)')'The input profile parameters are: c=',c1,' d=',c2
			!PRINT*
			
			!PRINT*,'	Best fit of least square:'
			!CALL FIT_lsq(rFEM, IFEM, u_lsq)
			!PRINT*
			!PRINT*,'	And for simulated LC:'
			!CALL FIT_lsq(rFEM, IFEM_SIM, u_lsq)	
			DEALLOCATE( rFEM, IFEM1, Sigma_I, IFEM2, mu, I_input, IFEM3 )
		ENDDO
		CLOSE(102)
		PRINT*,'Fitting results exported to ',F_test
		
	ENDIF
	t2=time()	
	dt_run=t2-t1
	print*,'		Program was running for ',dt_run,'s'
	CALL execute_command_line( 'say All done!')
	print *, char(7),char(7),char(7)
STOP		

END PROGRAM Limb
!*************************!
