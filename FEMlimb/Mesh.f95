MODULE Mesh
USE lensFs2
USE Events
IMPLICIT NONE

LOGICAL	:: I1eqI2, Verbal
REAL(DP), PARAMETER	:: ctm=0.4d0*dlog(1.d1)
REAL(SP), PARAMETER	:: Nsigma=3.0
	
CONTAINS 


FUNCTION sigma_l( a_u, a_rho, a_tE,li)
	REAL(DP), INTENT(IN)	:: a_u, a_rho, a_tE, li
	REAL(DP)	:: sigma_l
		! note that :a_u=s_u0*p/rho_st		a_rho=s_rho/rho_st
		!			a_tE=s_tE/tE

	sigma_l=sqrt( (a_u/li)**2	+  (a_rho*li)**2  + a_tE**2*( 1 - Event(ic)%p/li) &
															  *(1+Event(ic)%p/li) )
		
END FUNCTION sigma_l
!=======================!
SUBROUTINE  t_to_l_DATA(INfile, lmax, ti, li, mti, smi)
	
	CHARACTER(len=*), INTENT(IN)	:: INfile
	REAL(DP), INTENT(IN)	::	lmax
	REAL(DP), DIMENSION(:), Allocatable, INTENT(OUT)	:: ti, li, mti, smi	
			
	REAL(DP), DIMENSION(:), Allocatable	:: ti0, mt0, sm0, l0, resmi, resm0	
	INTEGER(I4B), ALLOCATABLE, DIMENSION(:)	:: index 
	
	REAL(DP)	:: dt, dt0, dtmin, dtmax, t_trans, tb, tf, &
					t, mt, sm, res0, delt2, ave_sm, ave_res,&
					l_i, dAl, a_u, a_rho, a_tE, ave_dAl, ave_off
	INTEGER	:: i,  j,  Nls, ierror
	
	t_trans= dsqrt(lmax-Event(ic)%p)*dsqrt(lmax+Event(ic)%p)*Event(ic)%ts
	tb= Event(ic)%t0 -  t_trans			
	tf= tb+2*t_trans
			
	OPEN(121, file=INfile)
	ierror=0
	j=0
	ave_sm=0.d0
	ave_res=0.d0
	DO
	  	READ(121,*,iostat=ierror) t, mt, sm, res0
	  	IF (ierror<0) EXIT	
			
	  	IF (t>=tb.and.t<=tf) THEN
	  		j=j+1
	  		ave_sm=ave_sm+sm 
	  		ave_res=ave_res+abs(res0)
	  	ENDIF	 
	ENDDO
	REWIND(121)
	Nls=j
	ave_sm=ave_sm/Nls
	ave_res=ave_res/Nls
	
	IF (Verbal) THEN
		PRINT*
		PRINT*, 'SUBROUTINE  t_to_l_DATA: data file:',INfile,' data counted',Nls
		PRINT*, 'SUBROUTINE  t_to_l_DATA: average of σ_m =',ave_sm
		PRINT*, 'SUBROUTINE  t_to_l_DATA: average of residual =',ave_res
		PRINT*
	ENDIF
			
	ALLOCATE(l0(1:Nls), index(1:Nls), ti0(1:Nls), mt0(1:Nls), sm0(1:Nls))
	ALLOCATE(li(1:Nls), ti(1:Nls), mti(1:Nls), smi(1:Nls), resmi(1:Nls),resm0(1:Nls) )
			
	ierror=0
	j=0
	i=0
	ave_off=0.d0
	DO
	  	READ(121,*,iostat=ierror) t, mt, sm, res0
	  	IF (ierror<0) EXIT	  			
	  	IF (t>=tb.and.t<=tf) THEN
	  		j=j+1
	  		ti0(j)=t
	  		delt2=(Event(ic)%t0-t)**2
	  		l0(j)=dsqrt(Event(ic)%p2+delt2/Event(ic)%ts2)
	  		mt0(j)=mt
	  		sm0(j)=sm
	  		resm0(j)=res0
	  		IF (res0>sm) THEN
	  			i=i+1
	  			ave_off=ave_off+res0/sm
	  		ENDIF
	  	END IF	
	 ENDDO
	 						
	 IF (i>0) ave_off=ave_off/i
	 
	 IF (Verbal) THEN
	 	print*
	 	print*, ' t_to_l_DATA SUB. average of Data residual/sigma_m > 1:',ave_off,&
	 			' number of data with residual>sigma',i
	 ENDIF	
			
	CLOSE(121)
			
	CALL indexx_dp(l0,index)			
			
	DO i=1, Nls
		j=index(i)
		li(i)=l0(j)
		mti(i)=mt0(j)
		smi(i)=sm0(j)
		resmi(i)=resm0(j)
		ti(i)=ti0(j)
	ENDDO
	
	IF (Verbal) THEN
		PRINT*, ' t_to_l_DATA SUB. max l',li(Nls)
		PRINT*
	ENDIF
				
END SUBROUTINE  t_to_l_DATA	
!=======================!
SUBROUTINE  m_to_A_DATA(iran, l_DATA, mti, smi, A_Data, dA_Data)!, A_model)
	
	REAL(DP), DIMENSION(:), INTENT(IN)	:: l_DATA, mti, smi	
	INTEGER, INTENT(IN)	::	iran
	REAL(DP), DIMENSION(:), Allocatable, INTENT(OUT)	:: A_Data, dA_Data
	REAL(DP), DIMENSION(:), Allocatable	:: A_model
	INTEGER	:: idum1
	INTEGER(I2B):: n_DATA
	
	n_DATA=size(l_DATA)
	ALLOCATE(A_DATA(1:n_DATA), dA_DATA(1:n_DATA), A_model(1:n_DATA))
	!iran=0 -> Data		
	A_Data=( 10**(0.4d0*( Event(ic)%mbl-mti )) - 1.d0 )/Event(ic)%blndng+1
	dA_Data=smi*( A_Data-1.d0 + 1/Event(ic)%blndng )*ctm	
	
	! iran/=0 -> Simulation
	IF (iran>0) THEN		
		idum1=-time()-iran !change it in each run
		IF (iran==10) idum1=-iran
		A_model=mFLUX(rmin, l_DATA, n_DATA)	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1)
    ELSE IF (iran==-1) THEN	
    	idum1=-time()-1 !change it in each run
		A_model=FLUX_step(rmin, l_DATA, n_DATA)	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1) 
    ELSE IF (iran==-2) THEN	
    	idum1=-time()-2 !change it in each run
		A_model=FLUX_tanh(rmin, l_DATA, n_DATA)	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1)    	   	   	
    ENDIF	 
				
END SUBROUTINE  m_to_A_DATA	
!=======================!
SUBROUTINE  SIM_A_DATA(iran, l_DATA, A_Data, dA_Data)
	
	REAL(DP), DIMENSION(:), INTENT(IN)	:: l_DATA
	INTEGER, INTENT(IN)	::	iran
	REAL(DP), DIMENSION(:), Allocatable, INTENT(OUT)	:: A_Data, dA_Data
	REAL(DP), DIMENSION(:), Allocatable	:: A_model
	INTEGER	:: i,idum1
	INTEGER(I2B):: n_DATA
	
	n_DATA=size(l_DATA)
	ALLOCATE(A_DATA(1:n_DATA), dA_DATA(1:n_DATA), A_model(1:n_DATA))

	! iran/=0 -> Simulation
	IF (iran>0) THEN		
		idum1=-100!-time()-iran
		IF (iran==10) idum1=-iran
		A_model=mFLUX(rmin, l_DATA, n_DATA)
		dA_Data=5.d-3*A_model	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1)

    ELSE IF (iran==-1) THEN	
    	idum1=-time()-1 !change it in each run
		A_model=FLUX_step(rmin, l_DATA, n_DATA)	
		dA_Data=5.d-3*A_model	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1)

    ELSE IF (iran==-2) THEN	
    	idum1=-200!-time()-1 !change it in each run
		A_model=FLUX_tanh(rmin, l_DATA, n_DATA)
		dA_Data=5.d-3*A_model	
    	A_Data=FRAND0(A_model, n_DATA, dA_Data, nsigma, idum1)	    	   	   	
    ENDIF
    	 
END SUBROUTINE  SIM_A_DATA	
	
!=======================!
SUBROUTINE data_groups(Nbox0, threshold ,l_all, setNos, used, min_a12)
IMPLICIT NONE	           
	REAL(DP), INTENT(IN)	:: threshold
	REAL(DP), DIMENSION(:), INTENT(IN)	:: l_all
	INTEGER, INTENT(IN)		:: Nbox0
	
	INTEGER, INTENT(OUT)		:: used
	INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT)	:: setNos
	REAL(DP), INTENT(OUT)	:: min_a12
	
	INTEGER, ALLOCATABLE, DIMENSION(:,:)	:: setNos0
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: useNo	
	REAL(DP), ALLOCATABLE, DIMENSION(:,:)	:: setls	
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: setNob, setNof, histogram
	
	REAL(DP)	::	dl0, deltal, med_a12
	CHARACTER*20	:: FMT1	
	LOGICAL		:: Quality, CHOSE
	INTEGER, PARAMETER	:: amax=10**3
	INTEGER		:: size0 ,i,j,k, ntot, idum1, idum2, minbar, &
				    idata, itot, imax, maxhist, Nbox, Nsetsf
! in this subroutine we create 	a histogram of data in terms of projected !
! impact parameter:l
! we want to know which data point belongs to which histogram box, so we  !
! save "begin and end numbers" of data points in each box of histogram in !
! arrays:setNob, setNof , so the histogram would be:setNof-setNob		  !
! we return the maximum height of the histogram i.e. Nbox too.			  !
! min_a12 is the minimum of all subset's med_δ=median(min(Δl_i,Δl_i+1)/max(Δl_i,Δl_i+1))
size0=size(l_all)
IF (Verbal) print*,size0,'SUBROUTINE data_groups, lmax:',l_all(size0)
min_a12=1.d2	
IF (Nbox0>1) THEN
	
	!Nbox0=int((1-Event(ic)%p)/dl0)+1
	deltal=1.d0-Event(ic)%p
	dl0=deltal/dfloat(Nbox0)

	ALLOCATE(useNo(1:size0))
	IF (Verbal) Write(*,*) 'Width of boxes of histogram=',dl0
	Write(FMT1,'(a,I4,a)') "(",Nbox0+1,"(I3,x))"

	ALLOCATE( setNob(1:Nbox0), setNof(1:Nbox0), histogram(1:Nbox0) )	
    setNob=0;	setNof=0;	histogram=0
    k=1
    DO i=1, size0-1
    	!j=nint( (l_all(i)-Event(ic)%p)/dl0)
    	j=int( (l_all(i)-Event(ic)%p)/dl0)+1
		!if (j==0) write(*,*) 'j=0 for ',Event(ic)%p, l_all(i)
    	histogram(j)=histogram(j)+1
	ENDDO
	j=int( (l_all(size0)-Event(ic)%p)/dl0)+1
	IF (j==Nbox0+1) THEN
		histogram(Nbox0)=histogram(Nbox0)+1
	ELSE IF (j<Nbox0+1) THEN
		histogram(j)=histogram(j)+1	
	ENDIF
!print*,histogram		
    imax=imaxloc(histogram)
    maxhist=histogram(imax)

	j=0
	Nbox=0
	DO i=1, Nbox0
	
		IF (histogram(i)/=0) then
		
			setNob(i)=j+1
			setNof(i)=setNob(i)+histogram(i)-1
			j=setNof(i)			
			!ntot=ntot*histogram(i)
			Nbox=Nbox+1
			
		ENDIF	
		
	ENDDO
	
	!Write(FMT1,'(a,I4,a)') "(",Nbox,"(I4,x))"

	!print*,'data_groups SUB, data histogram:',histogram
	ALLOCATE(setls(1:maxhist,1:Nbox),setNos0(1:maxhist,1:Nbox))
		
	setls=0.d0
	setNos0=0
	useNo=0
	i=0
	itot=0
	used=0
	
	DO WHILE(itot<amax*maxhist.and.i<maxhist.and.used<size0)
	
		idum1=-(itot+1)-time()  
		k=0
		i=i+1
		itot=itot+1
		
		DO j=1, Nbox0
		
			IF (histogram(j)/=0) THEN
			    CHOSE=.FALSE.
				DO WHILE (.NOT.CHOSE)
					idata=nint( ran2(idum1)*( histogram(j)-1 ) )+setNob(j)
					!idata=int( ran2(idum1)*histogram(j)+setNob(j) )
					IF (idata>size0.or.idata<1.or.idata<setNob(j).or.idata>setNof(j)) THEN
						print*,' Bad No. of data in module Mesh, SUB: data_groups',idata, &
							   ' accepted range is:',setNob(j),setNof(j)
						STOP
					ENDIF

					minbar=minval( useNo( setNob(j):setNof(j) ) )
				
					IF (useNo(idata)<=minbar) THEN
						CHOSE=.TRUE.
				    	K=K+1
				 		useNo(idata)=useNo(idata)+1
						setNos0(i,k)=idata				
						setls(i,k)=l_all(idata)
						IF ( k>1.and.setls(i,k)<setls(i,k-1) ) THEN
							print*, i,'Sub. data_groups: faulty point order',k,k-1,' ls:',setls(i,k-1), setls(i,k)
							print*, setls(i,:)
							STOP
						ENDIF
						!PRINT*, i,idata
					ENDIF	!useNo<=minbar
				
				ENDDO	!while chose
				
			ENDIF	!histogram>0
						
		END DO	!j

		IF (k<Nbox) THEN
			write(*,*) '!!!!! data_groups SUB. Only ', k,' No. of data assigned to the',i,&
					   'th group. it should be ',Nbox
			STOP
		ENDIF
		
		CALL mesh_quality(threshold, setls(i,:), med_a12, quality)
		min_a12=MIN(min_a12,med_a12)
		IF ( .NOT.quality) THEN
			DO j=1, Nbox
				useNo(setNos0(i,j))=useNo(setNos0(i,j))-1
			ENDDO
			setls(i,:)=0.d0
			setNos0(i,:)=0
			i=i-1
			!WRITE(*,*) itot,'th data group rejected.'			
		ENDIF
		
		used=0
		DO j=1, size0
			IF (useNo(j)>0) used=used+1
		ENDDO
		
	ENDDO	
	Nsetsf=i
		
	IF (Verbal.and.Nsetsf==0) THEN
		PRINT*
		PRINT*, ' data_groups SUB. : No subsets could be found after',itot,' attempts.'
		PRINT*, "                    minimum of all subset's med_δ=median(min(Δl_i,Δl_i+1)/max(Δl_i,Δl_i+1)) was:",min_a12
		PRINT*
	ELSE IF (Verbal) THEN
		PRINT*
		PRINT*,' data_groups SUB.: Max height of histogram=',maxhist
		PRINT*,'                   No. of data series: ',Nsetsf,' after ',itot,' attempts.'
		PRINT*,'                   No. of used data:',used,' out of',size0
		PRINT*
	ENDIF

	IF (Nsetsf>0) THEN					
		ALLOCATE(setNos(1:Nsetsf,1:Nbox))
		DO i=1, Nsetsf
			DO j=1, Nbox
				IF (setNos0(i,j)>0) THEN
				setNos(i,j)=setNos0(i,j)
				ELSE
					print*,' data_groups SUB. zero value in subset data numbers!!!!'
					STOP
				ENDIF		
			ENDDO
		ENDDO		
	ELSE
		ALLOCATE(setNos(1:1,1:1))
		setNos=setNos0(1:1,1:1)		
	ENDIF	

ELSE !Nbox0=1
	
	Nsetsf=1
	Nbox=1
	ALLOCATE(setNos(1:1,1:1))
	setNos(1,i)=1	
	
ENDIF
	
END SUBROUTINE data_groups
!=======================!
SUBROUTINE data_groups_mu(Nbox0, threshold ,l_all, setNos, used, min_a12)
IMPLICIT NONE	           
	REAL(DP), INTENT(IN)	:: threshold
	REAL(DP), DIMENSION(:), INTENT(IN)	:: l_all
	INTEGER, INTENT(IN)		:: Nbox0
	
	INTEGER, INTENT(OUT)		:: used
	INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT)	:: setNos
	REAL(DP), INTENT(OUT)	:: min_a12
	
	INTEGER, ALLOCATABLE, DIMENSION(:,:)	:: setNos0
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: useNo	
	REAL(DP), ALLOCATABLE, DIMENSION(:,:)	:: setls	
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: setNob, setNof, histogram
	REAL(DP), ALLOCATABLE, DIMENSION(:)	:: l_bins
	
	REAL(DP)	::	dl0, deltal, med_a12, dmi, mi
	CHARACTER*20	:: FMT1	
	LOGICAL		:: Quality, CHOSE
	INTEGER, PARAMETER	:: amax=10**3
	INTEGER		:: size0 ,i,j,k, ntot, idum1, idum2, minbar, &
				    idata, itot, imax, maxhist, Nbox, Nsetsf
! in this subroutine we create 	a histogram of data in terms of projected !
! impact parameter:l
! we want to know which data point belongs to which histogram box, so we  !
! save "begin and end numbers" of data points in each box of histogram in !
! arrays:setNob, setNof , so the histogram would be:setNof-setNob		  !
! we return the maximum height of the histogram i.e. Nbox too.			  !
! min_a12 is the minimum of all subset's med_δ=median(min(Δl_i,Δl_i+1)/max(Δl_i,Δl_i+1))
size0=size(l_all)
IF (Verbal) print*,size0,'SUBROUTINE data_groups_mu, lmax:',l_all(size0)
min_a12=1.d2	
IF (Nbox0>2) THEN
		ALLOCATE(l_bins(1:Nbox0))
		dmi=(1.d0-Event(ic)%p2)/(Nbox0-1)
		l_bins(1)=Event(ic)%p
		DO i=1, Nbox0-2
			mi=i*dmi
			l_bins(Nbox0-i)=dsqrt(1.d0-mi**2)
		ENDDO
		l_bins(Nbox0)=1.d0

	ALLOCATE(useNo(1:size0))
	IF (Verbal) Write(*,*) 'Width of boxes of histogram=',dl0
	Write(FMT1,'(a,I4,a)') "(",Nbox0+1,"(I3,x))"

	ALLOCATE( setNob(1:Nbox0-1), setNof(1:Nbox0-1), histogram(1:Nbox0-1) )	
    setNob=0;	setNof=0;	histogram=0
    k=1
    DO i=1, size0
    	DO j=2, Nbox0
    		IF (l_all(i)<l_bins(j)) EXIT
    	ENDDO
    	histogram(j-1)=histogram(j-1)+1
	ENDDO
	
	!DO i=1, Nbox0-1
	!	print*,i, l_bins(i), histogram(i)
	!ENDDO		
    imax=imaxloc(histogram)
    maxhist=histogram(imax)

	j=0
	Nbox=0
	DO i=1, Nbox0-1
	
		IF (histogram(i)/=0) then
		
			setNob(i)=j+1
			setNof(i)=setNob(i)+histogram(i)-1
			j=setNof(i)			
			!ntot=ntot*histogram(i)
			Nbox=Nbox+1
			
		ENDIF	
		
	ENDDO
	
	!Write(FMT1,'(a,I4,a)') "(",Nbox,"(I4,x))"
	!print*,'data_groups_mu SUB, data histogram_mu:',histogram
	ALLOCATE(setls(1:maxhist,1:Nbox),setNos0(1:maxhist,1:Nbox))
		
	setls=0.d0
	setNos0=0
	useNo=0
	i=0
	itot=0
	used=0
	
	DO WHILE(itot<amax*maxhist.and.i<maxhist.and.used<size0)
	
		idum1=-(itot+1)-time()  
		k=0
		i=i+1
		itot=itot+1
		
		DO j=1, Nbox0-1
		
			IF (histogram(j)/=0) THEN
			    CHOSE=.FALSE.
				DO WHILE (.NOT.CHOSE)
					idata=nint( ran2(idum1)*( histogram(j)-1 ) )+setNob(j)
					!idata=int( ran2(idum1)*histogram(j)+setNob(j) )
					IF (idata>size0.or.idata<1.or.idata<setNob(j).or.idata>setNof(j)) THEN
						print*,' Bad No. of data in module Mesh, SUB: data_groups',idata, &
							   'No. of bin:',j,' accepted range is:',setNob(j),setNof(j)
						STOP
					ENDIF

					minbar=minval( useNo( setNob(j):setNof(j) ) )
				
					IF (useNo(idata)<=minbar) THEN
						CHOSE=.TRUE.
				    	K=K+1
				 		useNo(idata)=useNo(idata)+1
						setNos0(i,k)=idata				
						setls(i,k)=l_all(idata)
						IF ( k>1.and.setls(i,k)<setls(i,k-1) ) THEN
							print*, i,'Sub. data_groups: faulty point order',k-1,k,' ls:',setls(i,k-1), setls(i,k)
							print*, setls(i,:)
							STOP
						ENDIF
						!PRINT*, i,idata
					ENDIF	!useNo<=minbar
				
				ENDDO	!while chose
				
			ENDIF	!histogram>0
						
		END DO	!j

		IF (k<Nbox) THEN
			write(*,*) '!!!!! data_groups SUB. Only ', k,' No. of data assigned to the',i,&
					   'th group. it should be ',Nbox
			STOP
		ENDIF
		
		CALL mesh_quality(threshold, setls(i,:), med_a12, quality)
		min_a12=MIN(min_a12,med_a12)
		IF ( .NOT.quality) THEN
			DO j=1, Nbox
				useNo(setNos0(i,j))=useNo(setNos0(i,j))-1
			ENDDO
			setls(i,:)=0.d0
			setNos0(i,:)=0
			i=i-1
			!WRITE(*,*) itot,'th data group rejected.'			
		ENDIF
		
		used=0
		DO j=1, size0
			IF (useNo(j)>0) used=used+1
		ENDDO
	ENDDO	
	Nsetsf=i
		
	!IF (Verbal.and.Nsetsf==0) THEN
		!PRINT*
		!PRINT*, ' data_groups_mu SUB. : No subsets could be found after',itot,' attempts.'
		!PRINT*, "                    minimum of all subset's med_δ=median(min(Δl_i,Δl_i+1)/max(Δl_i,Δl_i+1)) was:",min_a12
		!PRINT*
	!ELSE 
	IF (Verbal) THEN
		PRINT*
		PRINT*,' data_groups SUB.: Max height of histogram=',maxhist
		PRINT*,'                   No. of data series: ',Nsetsf,' after ',itot,' attempts.'
		PRINT*,'                   No. of used data:',used,' out of',size0
		PRINT*
	ENDIF

	IF (Nsetsf>0) THEN					
		ALLOCATE(setNos(1:Nsetsf,1:Nbox))
		DO i=1, Nsetsf
			DO j=1, Nbox
				IF (setNos0(i,j)>0) THEN
				setNos(i,j)=setNos0(i,j)
				ELSE
					print*,' data_groups SUB. zero value in subset data numbers!!!!'
					STOP
				ENDIF		
			ENDDO
		ENDDO		
	ELSE
		ALLOCATE(setNos(1:1,1:1))
		setNos=setNos0(1:1,1:1)		
	ENDIF	


ELSE !Nbox0=1
	
	Nsetsf=1
	Nbox=1
	ALLOCATE(setNos(1:1,1:1))
	setNos(1,i)=1	
	
ENDIF
	
END SUBROUTINE data_groups_mu

!=======================!
SUBROUTINE mesh_quality(threshold,li, med_a12, quality)
	
	REAL(DP), INTENT(IN)	:: threshold
	REAL(DP), DIMENSION(:), INTENT(IN)	:: li
	LOGICAL, INTENT(OUT)	:: quality
	REAL(DP), INTENT(OUT)	::	med_a12
	
	INTEGER		:: Nls, i, imed
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: index 
	REAL(DP)	:: invthresh, delta1, delta2, &
				   mina12!, ave_dl, min_dl
	REAL(DP), ALLOCATABLE, DIMENSION(:)	::	a12
	
	quality=.True.
	Nls=size(li)
	
	invthresh=1.d0/threshold
	mina12=1.d0
	med_a12=0.d0

	IF (Nls>2) THEN

		ALLOCATE(a12(1:Nls-2))
			
		DO i=2, Nls-1		
			delta1=li(i+1)-li(i)
			delta2=li(i)-li(i-1)
			
			a12(i-1)=min(delta1,delta2)/max(delta2,delta1)
		ENDDO
		
		CALL indexx_dp(a12,index)
		
		imed=INT(Nls/2.0)
		med_a12=a12( index(imed) )
		IF (med_a12<threshold) quality=.FALSE.

	ELSE
	
		quality=.FALSE.
		
	ENDIF				
		
		 	
		
END SUBROUTINE mesh_quality
!=======================!	
SUBROUTINE  mesh_l(rmin, li)
	
	REAL(DP), INTENT(IN)	:: rmin	
	REAL(DP), DIMENSION(:), INTENT(IN)	:: li	
	REAL(DP)	:: l_max						
	INTEGER	:: i, id, i1, i2, j, Nsim, ierror
	
	Nsim=size(li)							
	Nndt=Nsim
			
	IF (Nd1>1) THEN
		Nt=(Nndt-1)/(Nd1-1)
	ELSE
		Nt=Nndt
	END IF														
	ALLOCATE( dell(1:Nt), l(1:Nndt))
	ALLOCATE( nglobl(1:Nt,1:Nd1) )
			
	!------ Global Numbering ------!		
	DO i=1, Nt
		IF (Nd1>1) THEN	
			DO id=1, Nd1
				Nglobl(i,id)=(i-1)*(Nd1-1)+id
			ENDDO
		ELSE
			Nglobl(i,1)=i
		END IF		
	ENDDO	
	!------------------------------!
		
	l=li
	l_max=l(Nndt)
	!WRITE(*,*) 'mseh_l sub: data coverage in terms of l',l(1),l(Nndt)
		 							
	DO i=1, Nt
		i2=Nglobl(i,Nd1)
		i1=Nglobl(i,1)
		dell(i)=l(i2)-l(i1)
	ENDDO
	!------ Global Numbering of r-domain ------!
	IF (I1eqI2) THEN
		Nndr=Nndt+1 ! one node is added because of the I1=I2 constraint
	ELSE
		Nndr=Nndt
	END IF
										
	IF (Nd2>1) Nr=(Nndr-1)/(Nd2-1)		!Change it for Nd2=1
			
	!WRITE(*,*)'Mesh_l SUB.:No. of simulated data points and annuli: ',Nndt, Nndr
	ALLOCATE( delr(1:Nr), rv(1:Nndr) )
	ALLOCATE( Nglobr(1:Nr,1:Nd2) )
			
	DO i=1, Nr	
		IF (Nd2>1) THEN
			DO id=1, Nd2
				Nglobr(i,id)=(i-1)*(Nd2-1)+id
			ENDDO
		ELSE
				Nglobr(i,1)=i	
		ENDIF	
	ENDDO

	!------------------------------------------!  
		
	rv(1)=rmin	
	
	DO i=2, Nndt	
		rv(i)=l(i-1)
	ENDDO	
	rv(Nndr)=1 !because l_max may be less than one, but we need r=1
	
	IF (l_max<1-dell(Nt)) rv(Nndr-1)=( l_max+rv(Nndr-1) )/2
	IF (Nd2==3) THEN
		print*,'!!!!!!Note that this meshing works for solely equaly spaced data when using Nd2=3!!!!!!!!'
		print*,'!!!!!!because in second degree elements the 2nd Node of each element should be in midpoint of the element!!!!! '
		rv(2)=(rv(3)+rv(1))/2 
	ENDIF	  

	IF (Nd2>1) THEN
		DO i=1, Nr
			i2=Nglobr(i,Nd2)	
			i1=Nglobr(i,1)	
			delr(i)=rv(i2)-rv(i1)
		ENDDO
	ELSE
		DO i=1, Nndr-1	
			delr(i)=rv(i+1)-rv(i)
		ENDDO
	END IF
					
END SUBROUTINE  mesh_l
!=======================!
FUNCTION SHAPEm(xb, id)
	REAL(DP), INTENT(IN)     :: xb
	INTEGER(I1B), INTENT(IN) :: id
	REAL(DP)	:: SHAPEm

	SELECT CASE(id)
	CASE(1)
		SHAPEm=(1.d0-xb)/2.d0
	CASE(2)
		SHAPEm=(1.d0+xb)/2.d0
	END SELECT
		
END FUNCTION SHAPEm
!=======================!
SUBROUTINE indexx_dp(arr,index)
	!! Indexes an array arr(:), i.e., outputs the array indx(:) such that arr(indx(j)) is in ascending order for j = 1, 2, . . . , N . The input quantity arr is not changed.
	USE nrutil, ONLY : arth,assert_eq,nrerror,swap 
	IMPLICIT NONE
	
	REAL(DP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: index 
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(DP) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,ll,r 
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	 
	n=size(arr)!assert_eq(size(index),size(arr),'indexx_dp')
	ALLOCATE(index(1:n)) 
	index=arth(1,1,n)
	jstack=0 
	ll=1
	r=n
	do
		if (r-ll < NN) then 
			do j=ll+1,r
				indext=index(j) 
				a=arr(indext) 
				do i=j-1,ll,-1
					if (arr(index(i)) <= a) exit
					index(i+1)=index(i) 
				ENDDO
				index(i+1)=indext 
			ENDDO
			if (jstack == 0) RETURN 
			r=istack(jstack) 
			ll=istack(jstack-1) 
			jstack=jstack-2
		else 
			k=(ll+r)/2
			call swap(index(k),index(ll+1))
			call icomp_xchg(index(ll),index(r)) 
			call icomp_xchg(index(ll+1),index(r)) 
			call icomp_xchg(index(ll),index(ll+1)) 
			i=ll+1
			j=r
			indext=index(ll+1)
			a=arr(indext)
			do
				do
					i=i+1
					if (arr(index(i)) >= a) exit 
				ENDDO
				do
					j=j-1
					if (arr(index(j)) <= a) exit 
				ENDDO
				if (j < i) exit
				call swap(index(i),index(j)) 
			ENDDO
			index(ll+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small') 
			if (r-i+1 >= j-ll) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1 
			else
				istack(jstack)=j-1 
				istack(jstack-1)=ll 
				ll=i
			end if 
		end if
		!write(*,*) i, ll
	ENDDO 

CONTAINS

SUBROUTINE icomp_xchg(i,j) 
	INTEGER(I4B), INTENT(INOUT) :: i,j 
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
    	swp=i
    	i=j
    	j=swp
	end if
END SUBROUTINE icomp_xchg 
	
END SUBROUTINE indexx_dp 
!========================!	
	
END MODULE Mesh