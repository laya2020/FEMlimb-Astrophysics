MODULE FEMFs2
	USE nrutil
	USE lensFs2
	USE Mesh
	IMPLICIT NONE
	
CONTAINS

!=========================!
	SUBROUTINE GlM_limbr2(ng, MInty)
		INTEGER(I2B), INTENT(IN)	:: ng
		!INTEGER(I2B), PARAMETER	::	 ng2=100	
		REAL(DP), DIMENSION(1:Nndt,1:Nndr), INTENT(OUT)	:: MInty
		REAL(DP), DIMENSION(1:Nd1, 1:Nd1)	:: Gii, dumGii, invGii
		REAL(DP), DIMENSION(1:Nd1, 1:Nd2)	:: An1, In1, Anm, Inm
		REAL(DP), DIMENSION(1:ng)	:: li, wl
		REAL(DP), DIMENSION(1:ng)	:: r1, w1	
		REAL(DP), DIMENSION(1:1)	:: ystart
		INTEGER(I2B)	:: i, j, k, ig, jg, i1, i2, j1, j2, ib 
		INTEGER(I1B)	:: id, jd
		REAL(DP)	:: BF1, BF2, Jr, Jt,  rb, rc, dr, lb, lc, accuracy, h1,&
					 hmin, delta, MAG, sum1, ll, r,  INTEG, wi, wj, wk, sumA,&
					rb1, rb2, dw, func, cte,  tol, rbp, test									

		CALL gauleg(-1.d0, 1.d0, li, wl, ng)	
		CALL gauleg(-1.d0, 1.d0, r1, w1, ng)
	!------ dw is the width of rA(l,r)'s divergency in r=li  ------!

	!li=min(1.5d0*p,1.d0)
	!ll=p
	!func=(ll+1.d-8)*IntegA(ll,ll+1.d-8)
	!tol=1.d-8
	!delta=delr(1)/4
	!cte=min(2*IntegA(ll,1.d0)/func,0.5d0)!0.2d0
	!dw=dw_brent(ll,cte,func,1.d-8,1.d0,tol)
	!dw=min(dw, delta)
	!write(*,*)'cte:',cte, 'delr(1)/4, dw:	',delta, dw
	!write(*,*)'rA(l,r), r=l+1.d-8:',func, 'rA(l,r), r=l+dw:	',(ll+dw)*IntegA(ll,ll+dw)
	!IF (dw<0)	STOP'dw<0'
	!ib=nint(dw/dell(1))+1
	!write(*,*) 'dw:',dw, '	ib:',ib
	ib=nr !ib=nr => all integrations would be runge kutta
	!-------------------------------------------!					
		Gii=0.d0;	invGii=0.d0
		MInty=0.d0; Anm=0.d0; Inm=0.d0
	!write(*,*)	
	!write(*,*)' SUB. Glm_limbr2: Global Matrix is being calculated ... '	
IF (Nd1==3) THEN
	invGii(1,1)=4.5d0;			invGii(1,2)=-0.75d0; 	   invGii(1,Nd2)=1.5d0
	invGii(2,1)=invGii(1,2);	invGii(2,2)=1.125d0; 	   invGii(2,Nd2)=-0.75d0
	invGii(Nd1,1)=invGii(1,Nd2); 	invGii(Nd1,2)=invGii(2,Nd2); invGii(Nd1,Nd2)=4.5d0
	write(*,*) ' invGii done for Nd=3'
ELSE IF (Nd1==2) THEN
	invGii(1,1)=2.d0;	invGii(1,2)=-1.d0
	invGii(2,1)=-1.d0;	invGii(2,2)=2.d0
END IF 	
	!-----------------------------!	
	DO i=1, Nt
		!IF (MOD(i,10)==0) write(*,*) ' ',i,' rows calculated.'
		delta =dell(i)/2
		i2=Nglobl(i,Nd1);	i1=Nglobl(i,1)
		lc=( l(i2) + l(i1) )/2
		!write(*,*) i,dell(i)
				
		DO j=1, Nr
			dr =delr(j)/2
			j2=Nglobr(j,Nd2);	j1=Nglobr(j,1)
			rc =( rv(j2) + rv(j1) )/2
			
			DO id=1, Nd1
			DO jd=1, Nd2
				sum1=0.d0
		!	lb integration						
				DO ig=1, ng
					lb =li(ig)
					wi =wl(ig)
					ll =delta*lb+lc
					BF1 =SHAPE1D1(lb, id)
	   !	rb integration	
					ystart=0
					CALL odeint(jd,ll,dr,rc , ystart, -1.d0, 1.d0)!, accuracy, h1, hmin)					
					sum1=sum1+wi*BF1*ystart(1)
				END DO	! l integration				
	!-----------------------------!
	
		Anm(id, jd) =sum1
			END DO	!jd
			END DO	!id
	    !write(*,*) j,'Runge Kutta: ', sum1		
		Anm =dr*Anm
	!------- Gii^-1 Anm -----------!
		Inm=0.d0
		    DO id=1, Nd1
		    DO jd=1, Nd2
		    	sum1 =0.d0
		    	DO k=1, Nd1
		    		sum1 =sum1+invGii(id, k)*Anm(k,jd)
		    	END DO
		    	Inm(id, jd) = sum1	
		    END DO
		    END DO	
		        	
	!------- assembling -----------!
		!WRITE(*,*) 'Global numbering of elements: ', i, j
			DO id=1, Nd1
			DO jd=1, Nd2
				MInty(Nglobl(i,id), Nglobr(j,jd))=MInty(Nglobl(i,id), Nglobr(j,jd))+Inm(id, jd)			
			END DO
			END DO
	!------------------------------!					
		END DO		!j:	r
	END DO		!i: l

	!---	continuity condition of light curve	-----!
	DO i=1, Nt-1
		MInty((Nd1-1)*i+1,:)=MInty((Nd1-1)*i+1,:)/2
	END DO
	
	END SUBROUTINE GlM_limbr2	

!=========================!
	SUBROUTINE GlM_limbr1(MInty)
	IMPLICIT NONE
			
		REAL(DP), DIMENSION(1:Nndt,1:Nndr), INTENT(OUT)	:: MInty
		
		REAL(DP), DIMENSION(1:Nd1, 1:Nd2)	:: An1, In1, Anm, Inm

		REAL(DP), DIMENSION(1:1)	:: ystart
		INTEGER(I2B)	:: i, j, k, ig, jg, i1, i2, j1, j2
		INTEGER(I1B)	:: id, jd
		REAL(DP)	:: BF1, BF2, Jr, Jt,  rb, rc, dr, lb, lc, accuracy, h1,&
					 hmin, delta, MAG, sum1, ll, r,  INTEG, wi, wj, wk, sumA,&
					rb1, rb2, dw, func, cte,  tol, rbp, test									


		MInty=0.d0; Anm=0.d0; Inm=0.d0

	!-----------------------------!	
	DO i=1, Nndt
		!IF (MOD(i,10)==0) write(*,*) ' ',i,' rows calculated.'
		lc=l(i)
		!write(*,*) i,dell(i)
				
		DO j=1, Nr
			dr =delr(j)/2
			j2=Nglobr(j,Nd2);	j1=Nglobr(j,1)
			rc =( rv(j2) + rv(j1) )/2
			!write(*,*) '		',j,j1,j2

			DO jd=1, Nd2
	   !	rb integration  !	
					ystart=0
					CALL odeint(jd,lc,dr,rc , ystart, -1.d0, 1.d0)
					!CALL odeint(jd,lc,dr,rc , ystart, -1.d0, 1.d0)!, accuracy, h1, hmin)							
	   !    ------------    !
	   !write(*,*)  j, sum1
					Anm(1, jd) =ystart(1)
			END DO	!jd

	    !write(*,*) j,'Runge Kutta: ', sum1							
		Anm =dr*Anm

		        	
	!------- assembling -----------!
		!WRITE(*,*) 'Global numbering of r-elements: ', j,jd,' to ', Nglobr(j,jd)
			
			DO jd=1, Nd2
				MInty(i, Nglobr(j,jd))=MInty(i, Nglobr(j,jd))+Anm(1, jd)			
			END DO
			
	!------------------------------!					
		END DO		!j:	r
	END DO		!i: l

	
	END SUBROUTINE GlM_limbr1
!=========================!
	SUBROUTINE GlM_limbr0(MInty)
	IMPLICIT NONE
		
		REAL(DP), DIMENSION( 1:Nndt-1, 1:Nr), INTENT(OUT)	:: MInty
		
		REAL(DP)	:: An1, In1, Anm, Inm

		REAL(DP), DIMENSION(1:1)	:: ystart
		INTEGER(I2B)	:: i, j, k, ig, jg, i1, i2, j1, j2
		INTEGER(I1B)	:: id, jd
		REAL(DP)	:: BF1, BF2, Jr, Jt,  rb, rc, dr, lb, lc, accuracy, h1,&
					 hmin, delta, MAG, sum1, ll, r,  INTEG, wi, wj, wk, sumA,&
					rb1, rb2, dw, func, cte,  tol, rbp, test									


		MInty=0.d0; Anm=0.d0; Inm=0.d0

	!-----------------------------!	
	DO i=1, Nndt-1
		!IF (MOD(i,10)==0) write(*,*) ' ',i,' rows calculated.'
		lc=l(i)
		!write(*,*) i,dell(i)
				
		DO j=1, Nr
			dr =delr(j)/2
			j2=Nglobr(j,Nd2);	j1=Nglobr(j,1)
			rc =( rv(j2) + rv(j1) )/2
			!write(*,*) '		',j,j1,j2

			!DO jd=1, Nd2
	   !	rb integration  !	
					ystart=0
					jd=0
					CALL odeint(jd,lc,dr,rc , ystart, -1.d0, 1.d0)!, accuracy, h1, hmin)							
	   !    ------------    !
	   !write(*,*)  j, sum1
					Anm =ystart(1)
			!END DO	!jd

	    !write(*,*) j,'Runge Kutta: ', sum1							
		Anm =dr*Anm

		        	
	!------- assembling -----------!
		!WRITE(*,*) 'Global numbering of r-elements: ', j,jd,' to ', Nglobr(j,jd)
			
			
				MInty(i, j)=MInty(i, j)+Anm			
			
			
	!------------------------------!					
		END DO		!j:	r
	END DO		!i: l

	
	END SUBROUTINE GlM_limbr0	
!===========================!
    SUBROUTINE ainvert(a,n,ainv,flag)
    implicit none
      integer(I2B), intent(in) :: n
      integer, dimension(1:n) :: indx
      integer :: j,i
      REAL(DP), dimension(1:n,1:n), intent(inout) :: a
      REAL(DP), dimension(1:n,1:n), intent(out) ::ainv
      logical, intent(out)	:: flag
      REAL(DP) d

      DO i=1,n
        DO j=1,n
          ainv(i,j)=0.d0
        ENDDO
        ainv(i,i)=1.d0
      ENDDO

      call ludcmp(a,n,indx,d)
      
	  flag=.TRUE.
      DO j=1,n
        IF (abs(a(j,j))<1.e-18) THEN
        	WRITE(*,*) 'femFs2.f95, ainvert: zero determinant'
			flag=.FALSE.
        ENDIF	
        call lubksb(a,n,indx,ainv(:,j)) !in book : ainv(1,j)!!!!!!!!!!
      ENDDO
	
      END SUBROUTINE ainvert
!===========================!
 	SUBROUTINE lubksb(a,n,indx,b)
    IMPLICIT NONE
     	 INTEGER(I2B), INTENT(IN) :: N
         REAL(DP), DIMENSION(1:N,1:N), INTENT(IN) :: a
         INTEGER, DIMENSION(1:N), INTENT(IN) :: indx
         REAL(DP), DIMENSION(1:N), INTENT(INOUT) :: b
         INTEGER :: i,ii,ll
         REAL(DP) :: summ
         ii=0
         DO i=1,n
                 ll=indx(i)
                 summ=b(ll)
                 b(ll)=b(i)
                 IF (ii /= 0) THEN
                         summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
                 else IF (summ /= 0.0) THEN
                         ii=i
                 end if
                 b(i)=summ
         ENDDO
         
         DO i=n,1,-1
                 IF (a(i,i)/=0) b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)               
         ENDDO
    
	END SUBROUTINE lubksb
!===========================!
	SUBROUTINE ludcmp(a,n,indx,d)
	IMPLICIT NONE
    	INTEGER(I2B), INTENT(IN) :: n
    	REAL(DP), INTENT(OUT) :: d
        REAL(DP), DIMENSION(1:N,1:N), INTENT(INOUT) :: a
    	INTEGER, DIMENSION(1:n), INTENT(OUT) :: indx
    	REAL(DP), DIMENSION(1:n) :: vv
    	REAL(DP), PARAMETER :: EPS=1.D-18    
        INTEGER :: i, j, k, IMAX
    	REAL(DP) :: dum, AAMAX


            d = 1.0
            vv = 1.0 / MAXVAL(ABS(a),DIM=2)

            DO  j = 1 , n

               DO i = 2 , j-1
                  a(i,j) = a(i,j) - DOT_PRODUCT(a(i,1:i-1),a(1:i-1,j))

               ENDDO

          aamax = 0.0
               DO i = j , n

                  IF (j > 1) &
                  a(i,j) = a(i,j) - DOT_PRODUCT(a(i,1:j-1),a(1:j-1,j))
                  dum = vv(i) * ABS(a(i,j))
             IF (dum .GE. aamax) THEN
                     imax = i
                     aamax = dum
                  ENDIF
               ENDDO

          IF (j .NE. imax) THEN
                  DO  k = 1 , n

                     dum = a(imax,k)
                     a(imax,k) = a(j,k)
                     a(j,k) = dum
                  ENDDO
                  d = -d
                  vv(imax) = vv(j)

               ENDIF
               indx(j) = imax 

          IF (a(j,j) == 0.0) a(j,j) =EPS !TINY(a(j,j))
               IF (j < n) a(j+1:n,j) = a(j+1:n,j) / a(j,j)
         ENDDO

	END SUBROUTINE ludcmp
		
!==================================!

END MODULE femFs2