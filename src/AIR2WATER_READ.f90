
    SUBROUTINE read_calibration

USE commondata
USE ifport              ! necessario per il comando makedirqq

IMPLICIT NONE
INTEGER:: i, j, status
CHARACTER(LEN=1) :: string
CHARACTER(LEN=4) :: string1, string2, string3
LOGICAL result          ! necessario per il comando makedirqq

! read input information
OPEN(unit=1,file='input.txt',status='old',action='read')
READ(1,*)		! header
READ(1,*) name
READ(1,*) air_station
READ(1,*) water_station
READ(1,*) series
READ(1,*) time_res
READ(1,*) version
READ(1,*) Tice_cover
READ(1,*) fun_obj
READ(1,*) beta
READ(1,*) flag_icecal
READ(1,*) flag_snowfall 
READ(1,*) depth
READ(1,*) num_mod
READ(1,*) run
READ(1,*) prc
READ(1,*) n_run
READ(1,*) mineff_index
READ(1,*) log_flag
READ(1,*) CFL
CLOSE(1)

station=TRIM(air_station)//'_'//TRIM(water_station)

WRITE(string1,'(F4.2)') beta
WRITE(string2,'(I1)') flag_icecal
WRITE(string3,'(I1)') flag_snowfall
folder = TRIM(name)//'/output_'//TRIM(version)//'_beta'//TRIM(string1)//'_icecal'//TRIM(string2)//'_snowfall'//TRIM(string3)//'/'
result=makedirqq(folder)
    
WRITE(*,*) 'Objective function ',fun_obj

IF (run .eq. 'FORWARD') THEN
    OPEN(unit=1,file=TRIM(name)//'/parameters_forward.txt',status='old',action='read')    
    READ(1,*) (par(i), i=1,n_par)  
ELSE IF (run .eq. 'PSO') THEN    
    ! read PSO parameters
    OPEN(unit=1,file='PSO.txt',status='old',action='read')
    READ(1,*)		    ! header
    READ(1,*) n_particles
    READ(1,*) c1,c2
    READ(1,*) wmax,wmin
    CLOSE(1)
END IF

IF (run .eq. 'PSO' .or. run .eq. 'LATHYP') THEN    
    ! read model parameters
    OPEN(unit=1,file=TRIM(name)//'/parameters.txt',status='old',action='read')

    READ(1,*) (parmin(i),i=1,n_par);	
    READ(1,*) (parmax(i),i=1,n_par);	

    ! parameters that are not used are zeroed
    flag_par=.true.
    IF (version=='ICES') THEN                    ! air2water with ice
        parmin(7)=0;	parmax(7)=0;    flag_par(7)=.false.;
        parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.;
    ELSE
        parmin(9)=0;    parmax(9)=0;    flag_par(9)=.false.;
        parmin(10)=0;    parmax(10)=0;    flag_par(10)=.false.;
        parmin(11)=0;    parmax(11)=0;    flag_par(11)=.false.;
        parmin(12)=0;	parmax(12)=0;    flag_par(12)=.false.;
        IF (version=='6') THEN                  ! air2water with 6 parameters
        	parmin(7)=0;    parmax(7)=0;    flag_par(7)=.false.;
	        parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.;
	    ELSEIF (version=='4') THEN              ! air2water 4 parameters
		    parmin(5)=0;    parmax(5)=0;    flag_par(5)=.false.;
		    parmin(6)=0;	parmax(6)=0;    flag_par(6)=.false.;
        	parmin(7)=0;    parmax(7)=0;    flag_par(7)=.false.;
	        parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.;
        END IF
    END IF
        
    n_parcal=0    
    DO i=1,n_par
        IF (flag_par(i)==.true.) THEN
            n_parcal=n_parcal+1
        END IF
    END DO
    norm_min=SQRT(n_parcal*0.01)   ! 0.01 --> 1%
    
    CLOSE(1)
    	
    ! write parameters
    OPEN(unit=2,file=TRIM(folder)//'/parameters.txt',status='unknown',action='write')
    WRITE(2,'(<n_par>(F10.5,1x))') (parmin(i),i=1,n_par)
    WRITE(2,'(<n_par>(F10.5,1x))') (parmax(i),i=1,n_par)
    CLOSE(2)
    
    IF (log_flag==1) THEN
        parmin(2)=DLOG(parmin(2)); parmax(2)=DLOG(parmax(2));
        parmin(3)=DLOG(parmin(3)); parmax(3)=DLOG(parmax(3));
!        IF (version .ne. '4') THEN
!            parmin(5)=DLOG(parmin(5)); parmax(5)=DLOG(parmax(5));
!        END IF
    END IF
    
END IF

! Limits for numerical stability
IF (num_mod .eq. 'RK4') THEN
    LIM=2.785d0*CFL
ELSEIF (num_mod .eq. 'RK2' .or.  num_mod .eq. 'EUL' ) THEN
    LIM=2.0d0*CFL
END IF


! read T series (calibration)
CALL read_Tseries('c')

RETURN
END

!-------------------------------------------------------------------------------
!				LETTURA PERIODO VALIDAZIONE
!-------------------------------------------------------------------------------
SUBROUTINE read_validation

USE commondata

IMPLICIT NONE

!DEALLOCATE(date, tt, Tair, Twat_obs, Twat_obs_agg, Twat_mod, Twat_mod_agg, delta, hice, hice_obs, hice_mod, hice_obs_agg, hice_mod_agg)
DEALLOCATE(date, tt, Tair, Twat_obs, Twat_obs_agg, Twat_mod, Twat_mod_agg, delta, hice, hice_obs, hice_mod, hice_obs_agg, hice_mod_agg, hbice_obs_agg, hbice_mod_agg, hsice, hbice, hsice_obs_agg, hsice_mod_agg, hsnow, hsnow_obs, hsnow_mod, hsnow_obs_agg, hsnow_mod_agg, prec, hslush, hsice_obs, hsice_mod, hbice_obs, hbice_mod, Sf, Rf, snowfall) 
DEALLOCATE(I_pos, I_inf, I_posI, I_infI,I_posIb, I_infIb,I_posIs, I_infIs,I_posS, I_infS)
!DEALLOCATE(I_pos, I_inf, I_posI, I_infI, I_posS, I_infS)
! read T series (validation)
CALL read_Tseries('v')

RETURN
END

!-------------------------------------------------------------------------------
!				LETTURA FILE TEMPERATURA
!-------------------------------------------------------------------------------
SUBROUTINE read_Tseries(p)

USE commondata

IMPLICIT NONE

INTEGER :: i, j, k, status
INTEGER :: leap, year_ini
CHARACTER(LEN=1),INTENT(IN) :: p
CHARACTER(LEN=10) :: period

n_tot=0;

IF (p=='c') THEN
    period='calibration'
ELSE
    period='validation'
END IF

OPEN(unit=3,file=TRIM(name)//'/'//TRIM(station)//'_'//series//p//'.txt',status='unknown',action='read', iostat=status)
openif3: IF (status==0) THEN
	readloop3: DO
		READ(3,*,iostat=status)
		IF (status/=0) EXIT
		n_tot=n_tot+1
	END DO readloop3
	readif3: IF(status>0) THEN
	END IF readif3	
END IF openif3
REWIND(3)

! allocation + replication of the 1st year
WRITE(*,1001)  n_tot/365.25,TRIM(period)
1001 FORMAT('There are ',f4.1,' years for ', a12)

IF (p=='v' .and. n_tot .lt. 365) THEN
    WRITE(*,*) 'Validation period < 1 year --> validation is skipped'
    GO TO 100
END IF

n_year=CEILING(n_tot/365.25)
n_tot=n_tot+365             ! the 1st year is replicated. The 1st year is always considered 365 days long
ALLOCATE(date(n_tot,3),stat=status)
ALLOCATE(Tair(n_tot),stat=status)
ALLOCATE(Twat_obs(n_tot),stat=status) 
ALLOCATE(Twat_obs_agg(n_tot),stat=status) 
ALLOCATE(Twat_mod(n_tot),stat=status) 
ALLOCATE(Twat_mod_agg(n_tot),stat=status)
ALLOCATE(hsice_obs(n_tot),stat=status) 
ALLOCATE(hsice_mod(n_tot),stat=status) 
ALLOCATE(hsice_obs_agg(n_tot),stat=status) 
ALLOCATE(hsice_mod_agg(n_tot),stat=status) 
ALLOCATE(hbice_obs(n_tot),stat=status) 
ALLOCATE(hbice_mod(n_tot),stat=status) 
ALLOCATE(hbice_obs_agg(n_tot),stat=status) 
ALLOCATE(hbice_mod_agg(n_tot),stat=status) 
ALLOCATE(hice_obs(n_tot),stat=status) 
ALLOCATE(hice_mod(n_tot),stat=status)
ALLOCATE(hice_mod_agg(n_tot),stat=status)
ALLOCATE(hice_obs_agg(n_tot),stat=status)
!ALLOCATE(hbice_obs(n_tot),stat=status) 
!ALLOCATE(hbice_mod(n_tot),stat=status)
!ALLOCATE(hsice_obs(n_tot),stat=status) 
!ALLOCATE(hsice_mod(n_tot),stat=status)
ALLOCATE(hsnow_obs(n_tot),stat=status) 
ALLOCATE(hsnow_mod(n_tot),stat=status)
ALLOCATE(hsnow_obs_agg(n_tot),stat=status) 
ALLOCATE(hsnow_mod_agg(n_tot),stat=status)
ALLOCATE(hbice(n_tot),stat=status) 
ALLOCATE(hsice(n_tot),stat=status) 
ALLOCATE(hslush(n_tot),stat=status) 
ALLOCATE(tt(n_tot),stat=status)
ALLOCATE(delta(n_tot),stat=status) 
ALLOCATE(hice(n_tot),stat=status) 
ALLOCATE(hbice(n_tot),stat=status) 
ALLOCATE(hsnow(n_tot),stat=status) 
ALLOCATE(prec(n_tot),stat=status) 
ALLOCATE(Rf(n_tot),stat=status) 
ALLOCATE(Sf(n_tot),stat=status) 
ALLOCATE(snowfall(n_tot),stat=status) 

DO i=366,n_tot
	!READ(3,*) (date(i,j),j=1,3),Tair(i),Twat_obs(i),hice_obs(i)
    READ(3,*) (date(i,j),j=1,3),Tair(i),Twat_obs(i),hice_obs(i),hbice_obs(i),hsice_obs(i),hsnow_obs(i),prec(i), snowfall(i) 
END DO
date(1:365,:)=-999
Tair(1:365)=Tair(366:730)
Twat_obs(1:365)=Twat_obs(366:730)
hice_obs(1:365)=hice_obs(366:730)
hbice_obs(1:365)=hbice_obs(366:730)
hsice_obs(1:365)=hsice_obs(366:730)
hsnow_obs(1:365)=hsnow_obs(366:730)
prec(1:365)=prec(366:730)
snowfall(1:365)=snowfall(366:730) 

CLOSE(3)

year_ini=date(366,1)

! check leap years + define tt
k=0
DO j=1,365
    tt(k+j)=REAL(j)/365.0d0
END DO
k=365
DO i=1,n_year
    CALL leap_year(year_ini+i-1,leap)
    IF(leap==0) THEN
        DO j=1,365
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/365.0d0
        END DO
        k=k+365
    ELSE
        DO j=1,366
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/366.0d0
        END DO
        k=k+366
    END IF
END DO


100 RETURN 
END