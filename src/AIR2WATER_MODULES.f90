MODULE commondata
IMPLICIT NONE
SAVE 

INTEGER, PARAMETER :: n_par = 12
REAL(KIND=8), PARAMETER :: pi = ACOS(0.d0)*2.d0
REAL(KIND=8), PARAMETER :: DoY = 1.0d0/365.0d0

REAL(KIND=8), PARAMETER :: rhow=1000.0d0
REAL(KIND=8), PARAMETER :: rhoi=917.0d0
REAL(KIND=8), PARAMETER :: rhos=300.0d0
REAL(KIND=8), PARAMETER :: rhosl=rhow+rhos*(1-rhow/rhoi)
REAL(KIND=8), PARAMETER :: cp=4186.0d0
REAL(KIND=8), PARAMETER :: Lf=334000.0d0
REAL(KIND=8), PARAMETER :: ni=1-rhos/rhoi
REAL(KIND=8), PARAMETER :: E=0.001d0
REAL(KIND=8), PARAMETER :: ki=2.0d0
REAL(KIND=8), PARAMETER :: ksnow=0.30d0

INTEGER :: ii
INTEGER :: n_tot, n_dat, n_datI ,n_datIb, n_datIs, n_datS!, n_totS
INTEGER :: qty,log_flag
INTEGER :: n_year
INTEGER :: n_parcal
INTEGER :: flag_icecal, flag_snowfall 
INTEGER, ALLOCATABLE, DIMENSION(:) :: I_pos, I_posI, I_posIb, I_posIs, I_posS
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_inf, I_infI, I_infIb, I_infIs, I_infS
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: date
REAL(KIND=8) :: Tice_cover, prc, beta, depth
REAL(KIND=8) :: mean_obs, TSS_obs, std_obs
REAL(KIND=8) :: meanI_obs, TSSI_obs, stdI_obs
REAL(KIND=8) :: meanIb_obs, TSSIb_obs, meanIs_obs, TSSIs_obs 
REAL(KIND=8) :: meanS_obs, TSSS_obs !, stdS_obs
REAL(KIND=8) :: Tmin
REAL(KIND=8) :: mineff_index,finalfit
REAL(KIND=8) :: CFL, LIM
REAL(KIND=8) :: norm_min
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: tt, Tair, Twat_obs_agg, Twat_obs
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Twat_mod, Twat_mod_agg, delta, hice,hbice,hsice
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: hice_mod, hice_obs, hice_obs_agg, hice_mod_agg
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: hbice_obs,hsice_obs,hbice_mod,hsice_mod,hbice_obs_agg,hsice_obs_agg,hbice_mod_agg,hsice_mod_agg
REAL(KIND=8) :: hsnowst, csnow, cice, cslush, Dslushf, Dslushm, DslushRf, Dhicetot
REAL(KIND=8) :: eps, fsin, fsinn, zeta, denbi, Bbi, Cbi, gammas, hsnowmelt, Dhice, Rnn, sur, gammai, denbinn
REAL(KIND=8) :: phi, psi, theta
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: hsnow, hsnow_mod, hsnow_obs, hsnow_mod_agg, hsnow_obs_agg
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: hslush, Sf, Rf, prec, snowfall 

REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmin
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmax
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: par, par_best
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: m, q, r2
LOGICAL,ALLOCATABLE,DIMENSION(:) :: flag_par

CHARACTER(len=100) :: folder
CHARACTER(LEN=60) :: name, air_station, water_station, station, run
CHARACTER(LEN=1) :: series, unit
CHARACTER(LEN=3) :: time_res
CHARACTER(LEN=3) :: fun_obj, num_mod
CHARACTER(LEN=4) :: version

INTEGER :: n_run,n_particles      ! PSO
REAL(KIND=8) ::c1,c2,wmin,wmax    ! PSO

END MODULE commondata