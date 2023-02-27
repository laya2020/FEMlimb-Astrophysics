MODULE events
	USE Nrtype
	IMPLICIT NONE
	
 	CHARACTER(LEN=75)	::		DATA_PATH="./DataBase/"
 	TYPE microlensing
		CHARACTER(LEN=55), ALLOCATABLE, DIMENSION(:)	::	tel_file
		REAL(DP), ALLOCATABLE, DIMENSION(:)	::	u1fit, du1
		CHARACTER(LEN=150)	::  event_fldr, fit_file
		INTEGER				::	nTel, Nref, Nref_figs
		REAL(DP)			::  t0, u0, s_u0, tE, s_tE, rho_st, s_rho, El, El2, mbl2, blndng2 &
								, ts, ts2, p, p2, blndng, mbl, sigma_l, blndng_Choi, mbl_Choi								
	END TYPE microlensing
 										  
 	REAL(DP)	::	FBL, FBL0 =PI_d 	!FBL=PI_d*(1-climb1/3-climb2/5)	
 	INTEGER		::	ic, nLc
	
	TYPE(microlensing), DIMENSION(1:10), SAVE	:: event !the 10'th one is simulayion
	
CONTAINS
		
	SUBROUTINE EVENTS_info
	IMPLICIT NONE
		INTEGER	::	i
		
		nLc=size(event)
		
		event%ntel=(/10,11,10,4,8,7,8,7,12,1/)
		event%nref=(/1, 1, 1, 1, 2, 2, 1, 1, 1, 1/)	!due to the article saying: 
													!OGLE or MOA if both data presented OGLE
		event%nref_figs=(/1,1,2,1,1,6,1,6,1,1/)	!due to the Choi's email saying: the first 											
		!observatory listed in each figure is the reference observatory primary parameters		
		event%t0=(/4245.056d0, 4289.269d0, 4963.816d0, 5395.791d0, 5678.555d0,  &
				   				 5742.005d0, 3166.823d0, 5758.691d0, 5823.574d0,0.d0/)
				   				
		event%u0=(/0.0363d0, 0.0060d0, 0.0005d0, 0.0002d0, 0.0292d0, 0.0029d0,  &
				    					           0.0111d0, 0.0151d0, 0.0485d0,1.d-4/)
				    					
		event%s_u0=(/0.0005d0, 0.0002d0, 0.0001d0, 0.0002d0, 0.0002d0, .0001d0, &
				     							   0.0004d0, 0.0004d0, 0.0005d0, 0.d0/)
		    		   
		event%tE=(/8.13d0, 15.90d0, 64.99d0, 12.78d0, 14.97d0, 2.65d0,12.84d0,  &
																6.70d0, 29.06d0,7.29d0/)
																
		event%s_tE=(/0.07d0, 0.05d0, 0.61d0, 1.08d0, 0.05d0, 0.06d0, 0.09d0,    &
																 0.07d0, 0.11d0,0.d0/)
															
		event%rho_st=(/0.0590d0, 0.0364d0, 0.0020d0, 0.0041d0,0.0538d0,0.0129d0 &
												 , 0.0418d0, 0.0199d0, 0.0979d0, 0.02d0/)
												
		event%s_rho=(/0.0006d0, 0.0001d0, 0.0001d0, 0.0003d0,0.0002d0,0.0003d0, &
												   0.0004d0, 0.0003d0, 0.0006d0,0.d0/)
		event%El	= 1/event%rho_st
		event%El2	= event%El**2
		event%ts	= event%rho_st*event%tE
		event%ts2	= event%ts**2
		event%p		= event%u0/event%rho_st
		event%p2	= event%p**2
										   
!fitted values of m_bl and blending for light curves of ref_figs, t<3tE, sigma_m<0.1, res_m<3*sigma_m
!												  
		event%mbl2=(/17.462d0, 16.307d0, 19.267d0, 16.948d0, 16.691d0, 15.866d0, 16.33d0,&
															   18.484d0, 15.178d0,15.d0/)
		
		event%blndng2=(/0.54d0, 1.018d0, 0.810d0, 0.026d0, 1.781d0, 0.056d0, 1.040d0 &
															 &, 0.983d0, 1.039d0,1.d0/)
!fitted values of m_bl and blending for light curves of ref_figs by minimising 1-Chi2 and res<3sigma
!											  
		event%mbl=(/17.462d0, 16.307d0, 19.506d0, 16.948d0, 16.68d0, 15.874d0, 16.326d0,&
															   18.734d0, 15.178d0,15.d0/)
		
		event%blndng=(/0.54d0, 1.018d0, 1.0d0, 0.026d0, 1.759d0, 0.057d0, 1.023d0 &
															 &, 1.23d0, 1.039d0,1.d0/)
															 
!fitted values of m_bl and blending by lsq fitting of m_Choi_best_fit and A_choi_best_fit

		event%mbl_Choi=(/17.458d0, 16.307d0, 19.264d0, 16.948d0, 16.69d0, 15.874d0, 16.329d0,&
															   18.477d0, 15.178d0,15.d0/)
		
		event%blndng_Choi=(/0.538d0, 1.018d0, .806d0, 0.026d0, 1.783d0, 0.057d0, 1.04d0 &
															 &, 0.968d0, 1.040d0,1.d0/)
															 
															 
		event%sigma_l=(/0.01d0, 0.004d0, 0.042d0, 0.092d0, 0.005d0, 0.02d0, .01d0 &
																,0.023d0, 0.06d0,0.d0/)															   														
				
		event(1)%event_fldr="MB07176_data/"		
		event(2)%event_fldr="MB07233+OB07302_data/"		
		event(3)%event_fldr="MB09174_data/"
		event(4)%event_fldr="MB10436_data/"
		event(5)%event_fldr="MB11093_data/"
		event(6)%event_fldr="MB11274_data/"
		event(7)%event_fldr="OB04254_data/"		
		event(8)%event_fldr="OB110990+MB11300_data/"
		event(9)%event_fldr="OB111101+MB11325_data/"
		event(10)%event_fldr="_data/"
		
		event(1)%fit_file="plots/MB07176_LC.res"
		event(2)%fit_file="plots/MB07233+OB07302_LC.res"
		event(3)%fit_file="plots/MB09174_LC.res"
		event(4)%fit_file="plots/MB10436_LC.res"
		event(5)%fit_file="plots/MB11093_LC.res"
		event(6)%fit_file="plots/MB11274_LC.res"
		event(7)%fit_file="plots/OB04254_LC.res"		
		event(8)%fit_file="plots/OB110990+MB11300_LC.res"
		event(9)%fit_file="plots/OB111101+MB11325_LC.res"

		DO i=1, nLc
			ALLOCATE(event(i)%tel_file(1:event(i)%ntel), event(i)%u1fit(1:event &
									   (i)%ntel), event(i)%du1(1:event(i)%ntel) )
		ENDDO
		
		event(1)%u1fit=(/0.53d0,0.5d0,0.6d0, 0.51d0,0.6d0, 0.6d0, 0.6d0, 0.44d0, 0.6d0, 0.6d0/) !it had 11 elements, probably this was the source of segmentation fault
		event(2)%u1fit=(/0.53d0,0.56d0,0.6d0,.49d0,0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.56d0, 0.6d0, 0.6d0/)
		event(3)%u1fit=(/0.d0, 0.33d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0/)
		event(4)%u1fit=(/0.52d0, 0.6d0, 0.6d0, 0.6d0/)
		event(5)%u1fit=(/0.55d0, 0.51d0, 0.58d0, 0.69d0, 0.51d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0/)
		event(6)%u1fit=(/0.6d0, 0.6d0, 0.6d0, 0.48d0, 0.6d0, 0.51d0, 0.6d0/)
		event(7)%u1fit=(/0.7d0, 0.56d0, 0.d0, 0.69d0, 0.78d0, 0.d0, 0.55d0, 0.d0/)
		event(8)%u1fit=(/0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.56d0, 0.6d0/)
		event(9)%u1fit=(/0.74d0, 0.77d0, 0.81d0, 0.89d0, 0.6d0, 0.77d0, 0.78d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0/)
		event(10)%u1fit=(/0.6d0/)
		!If a du1 is 0.d0 it means that u1 was not derived for that LC, so we take it as u1=0.6±0 for simulation perposes
		event(1)%du1=(/0.04d0, 0.05d0, 0.d0, 0.05d0, 0.d0, 0.d0, 0.d0, 0.06d0, 0.d0, 0.d0/)
		event(2)%du1=(/0.04d0, 0.02d0, 0.d0, 0.02d0, 0.d0, 0.d0, 0.0d0, 0.d0, 0.04d0, 0.d0, 0.d0/)
		event(3)%du1=(/0.d0, 0.02d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
		event(4)%du1=(/0.1d0, 0.d0, 0.d0, 0.d0/)
		event(5)%du1=(/0.04d0, 0.1d0, 0.04d0, 0.05d0, 0.03d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
		event(6)%du1=(/0.d0, 0.d0, 0.d0, 0.02d0, 0.d0, 0.01d0, 0.d0/)
		event(7)%du1=(/0.07d0, 0.1d0, 0.d0, 0.1d0, 0.09d0, 0.d0, 0.06d0, 0.d0/)
		event(8)%du1=(/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.04d0, 0.d0/)
		event(9)%du1=(/0.07d0, 0.08d0, 0.07d0, 0.14d0, 0.d0, 0.06d0, 0.05d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
		event(10)%du1=(/0.d0/)
		
		event(1)%tel_file(1)="O_MOA_1.dat"
		event(1)%tel_file(2)="O_CTIO_2.dat"
		event(1)%tel_file(3)="O_CTIO_3.dat"
		event(1)%tel_file(4)="O_Auckland_4.dat"
		event(1)%tel_file(5)="O_Canopus_5.dat"
		event(1)%tel_file(6)="O_Danish_6.dat"
		event(1)%tel_file(7)="O_FCO_7.dat"
		event(1)%tel_file(8)="O_Lemmon_8.dat"
		event(1)%tel_file(9)="O_Steward_9.dat"
		event(1)%tel_file(10)="O_VLO_A.dat"
		
		event(2)%tel_file(1)="O_OGLE_1.dat"
		event(2)%tel_file(2)="O_CTIO_2.dat"
		event(2)%tel_file(3)="O_CTIO_3.dat"
		event(2)%tel_file(4)="O_Canopus_4.dat"
		event(2)%tel_file(5)="O_Canopus_5.dat"
		event(2)%tel_file(6)="O_Danish_6.dat"
		event(2)%tel_file(7)="O_FCO_7.dat"
		event(2)%tel_file(8)="O_Lemmon_8.dat"
		event(2)%tel_file(9)="O_MOA_9.dat"
		event(2)%tel_file(10)="O_Perth_A.dat"
		event(2)%tel_file(11)="O_SAAO_B.dat"
		
		event(3)%tel_file(1)="O_MOA_1.dat"
		event(3)%tel_file(2)="O_CTIO_2.dat"
		event(3)%tel_file(3)="O_CTIO_3.dat"
		event(3)%tel_file(4)="O_Bronberg_4.dat"
		event(3)%tel_file(5)="O_CAO_5.dat"
		event(3)%tel_file(6)="O_Canopus_6.dat"
		event(3)%tel_file(7)="O_Craigie_7.dat"
		event(3)%tel_file(8)="O_Kumeu_8.dat"
		event(3)%tel_file(9)="O_LT_9.dat"
		event(3)%tel_file(10)="O_Possum_A.dat"
		
		event(4)%tel_file(1)="O_MOA_1.dat"
		event(4)%tel_file(2)="O_FTS_2.dat"
		event(4)%tel_file(3)="O_SAAO_3.dat"
		event(4)%tel_file(4)="O_SAAO_4.dat"
		
		event(5)%tel_file(1)="O_MOA_1.dat"
		event(5)%tel_file(2)="O_OGLE_2.dat"
		event(5)%tel_file(3)="O_CTIO_3.dat"
		event(5)%tel_file(4)="O_CTIO_4.dat"
		event(5)%tel_file(5)="O_Canopus_5.dat"
		event(5)%tel_file(6)="O_Perth_6.dat"
		event(5)%tel_file(7)="O_FTN_7.dat"
		event(5)%tel_file(8)="O_FTS_8.dat"
		
		event(6)%tel_file(1)="O_MOA_1.dat"
		event(6)%tel_file(2)="O_OGLE_2.dat"
		event(6)%tel_file(3)="O_CTIO_3.dat"
		event(6)%tel_file(4)="O_Kumeu_4.dat"
		event(6)%tel_file(5)="O_Perth_5.dat"
		event(6)%tel_file(6)="O_Auckland_6.dat"
		event(6)%tel_file(7)="O_FCO_7.dat"
		
		event(7)%tel_file(1)="O_OGLE_1.dat"
		event(7)%tel_file(2)="O_CTIO_2.dat"
		event(7)%tel_file(3)="O_CTIO_3.dat"
		event(7)%tel_file(4)="O_Boyden_4.dat"
		event(7)%tel_file(5)="O_Canopus_5.dat"
		event(7)%tel_file(6)="O_FCO_6.dat"
		event(7)%tel_file(7)="O_SAAO_7.dat"
		event(7)%tel_file(8)="O_Danish_8.dat"
		
		event(8)%tel_file(1)="O_OGLE_1.dat"
		event(8)%tel_file(2)="O_MOA_2.dat"
		event(8)%tel_file(3)="O_SAAO_3.dat"
		event(8)%tel_file(4)="O_SAAO_4.dat"
		event(8)%tel_file(5)="O_Canopus_5.dat"
		event(8)%tel_file(6)="O_Pico_6.dat"
		event(8)%tel_file(7)="O_Possum_7.dat"
		
		event(9)%tel_file(1)="O_OGLE_1.dat"
		event(9)%tel_file(2)="O_MOA_2.dat"
		event(9)%tel_file(3)="O_CTIO_3.dat"
		event(9)%tel_file(4)="O_CTIO_4.dat"
		event(9)%tel_file(5)="O_Auckland_5.dat"
		event(9)%tel_file(6)="O_Canopus_6.dat"
		event(9)%tel_file(7)="O_FTS_7.dat"
		event(9)%tel_file(8)="O_FTN_8.dat"
		event(9)%tel_file(9)="O_LT_9.dat"
		event(9)%tel_file(10)="O_Possum_A.dat"
		event(9)%tel_file(11)="O_SSO_B.dat"
		event(9)%tel_file(12)="O_VLO_C.dat"
		
		event(10)%tel_file(1)="sim.dat"
		
		write(*,*)''//achar(27)//'[95m Parameters of ',nLc,'light-curves importe&
				  &d'//achar(27)//'[0m'
		
								  
	END SUBROUTINE EVENTS_info
	

END MODULE events
