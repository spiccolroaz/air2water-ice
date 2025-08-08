    !-------------------------------------------------------------------------------

    !				sub_1
    !-------------------------------------------------------------------------------
    SUBROUTINE sub_1(ei,eiT,eiI,eiIb,eiIs,eiS)
    USE commondata

    ! Dichiarazione variabili
    implicit none
    INTEGER:: i, j, k
    REAL(KIND=8) :: ei
    REAL(KIND=8) :: eiT,eiI,eiIb,eiIs,eiS
    REAL(KIND=8) :: ind, TSS, mean_mod_ymo,mean_obs_ymo
    REAL(KIND=8) :: TSSI

    ! Richiamo modello conversione
    CALL model

    ! Calcolo dell'indice di efficienza
    CALL funcobj(ei,eiT,eiI,eiIb,eiIs,eiS)

    RETURN
    END

    !-------------------------------------------------------------------------------
    !				MODELLO
    !-------------------------------------------------------------------------------
    SUBROUTINE model
    USE commondata

    ! Dichiarazione variabili
    IMPLICIT NONE

    INTEGER :: j, k
    INTEGER :: nsub
    REAL(KIND=8) :: DD, Tw_neg
    REAL(KIND=8) :: K1, K2, K3, K4
    REAL(KIND=8) :: K1ice, K2ice, K3ice, K4ice
    REAL(KIND=8) :: pp, lambda, dt
    REAL(KIND=8) :: dTair , ttt
    REAL(KIND=8) :: Tairk, Tairk1, Twatk, ttk, hicek
    REAL(KIND=8) :: alpha, d1, d2, A, B, C, gamma, fn, fnn
    REAL(KIND=8) :: ratio, Dhsice, fluxes, Tairnewt, Tairnewtt, z, zz

    Tmin=MAX(0.0d0,MINVAL(Twat_obs))
    Tmin=MAX(4.d0,Tmin)	! Computation of Tmin

    IF (Twat_obs(1)==-999) THEN	! Initial condition
        Twat_mod(1)=Tmin
    ELSE
        Twat_mod(1)=Twat_obs(1)
    END IF

    hice=0.0d0 ! Initial condition
    hbice=0.0d0 ! Initial condition      
    hslush=0.0d0 ! Initial condition    
    hsice=0.0d0 ! Initial condition     
    hsnow=0.0d0
    
    DO j=1,n_tot-1

        IF (num_mod .eq. 'CRN') THEN
            dt=1.0d0
            CALL a2w(Tair(j), Twat_mod(j), tt(j), K1, DD, pp)
            delta(j)=DD
            Twatk=( Twat_mod(j)*2.0d0*DD +                     &
                dt*( pp + par(1) + par(2)*Tair(j+1) + par(5)*cos(2.d0*pi*(tt(j)+ttt-par(6))) ) ) &
                /(2.0d0*DD+dt*par(3)) ;
            !Tw_neg=MIN(0.0d0,Twatk)
            Twatk=MAX(Twatk,Tice_cover)
        ELSE
            ! Substepping procedure
            CALL a2w(Tair(j), Twat_mod(j), tt(j), K1, DD, pp)
            DD=MAX(DD,0.01d0)
            lambda=(pp/par(4) - par(3))/DD
            pp=-LIM/lambda
            IF (lambda .le. 0.0d0 .and. pp .lt. 1.0d0) THEN
                nsub=CEILING(1.0d0/pp)
                nsub=MIN(100,nsub)
                dt=1.0d0/nsub
            ELSE
                dt=1.0d0
                nsub=1
            END IF
            dTair=(Tair(j+1)-Tair(j))/DBLE(nsub)
            ttt = DoY/DBLE(nsub)

            Twatk=Twat_mod(j)
            DO k=1,nsub ! Substepping cycle
                Tairk = Tair(j) + dTair*DBLE(k-1)
                Tairk1= Tairk + dTair
                ttk=tt(j) + ttt*DBLE(k-1)

                IF (num_mod .eq. 'RK4') THEN
                    CALL a2w(Tairk, Twatk, ttk, K1, DD, pp)
                    delta(j)=DD
                    CALL a2w(0.5d0*(Tairk + Tairk1), Twatk + 0.5d0*K1, ttk + 0.5*ttt, K2, DD, pp)
                    CALL a2w(0.5d0*(Tairk + Tairk1), Twatk + 0.5d0*K2, ttk + 0.5*ttt, K3, DD, pp)
                    CALL a2w(Tairk1, Twatk + K3, ttk + ttt, K4, DD, pp)

                    Twatk=Twatk + 1.0d0/6.0d0*(K1 + 2.0d0*K2 + 2.0d0*K3 + K4 )*dt
                ELSEIF (num_mod .eq. 'RK2') THEN
                    CALL a2w(Tairk, Twatk, ttk, K1, DD, pp)
                    delta(j)=DD
                    CALL a2w(Tairk1, Twatk + K1, ttk + ttt, K2, DD, pp)

                    Twatk=Twatk + 0.5d0*(K1 + K2)*dt
                ELSEIF (num_mod .eq. 'EUL') THEN
                    CALL a2w(0.5d0*(Tairk + Tairk1), Twatk, ttk, K1, DD, pp)
                    delta(j)=DD

                    Twatk=Twatk + K1*dt
                ELSE
                    WRITE(*,*) 'Error in the choice of the numerical model'
                    STOP
                END IF
                !Tw_neg=MIN(Twatk,Tw_neg,0.0d0)
                Twatk=MAX(Twatk,Tice_cover)

            END DO
        END IF

        Twat_mod(j+1)=Twatk
        Tmin=MAX(4.d0,MIN(Tmin,Twat_mod(j+1)))

        IF (version=='ICES' .and. (Twat_mod(j+1) .le. Tice_cover .or. hice(j) .gt. 0.0d0)) THEN
            Twat_mod(j+1)=Tice_cover
            hbice(j)=MAX(hbice(j),1D-4)
            
            IF (flag_snowfall ==0) THEN !SNOWFALL
                !IF (Tair(j) .le. (par(9) + par(12))) THEN
                IF (Tair(j) .le. (par(9)+par(12))) THEN
                    !IF (Tair(j) .le. par(12)) THEN
                    Sf(j)=prec(j) ! m H20
                    Rf(j)=0.0d0
                ELSE
                    Sf(j)=0.0d0
                    Rf(j)=prec(j)
                END IF

                IF (Tair(j+1) .le. (par(9)+par(12))) THEN
                    Sf(j+1)=prec(j+1) ! m H20
                    Rf(j+1)=0.0d0
                ELSE
                    Sf(j+1)=0.0d0
                    Rf(j+1)=prec(j+1)
                END IF
                
            ELSE IF (flag_snowfall ==1) THEN
                Sf(j)=snowfall(j) 
                Rf(j)=prec(j)-snowfall(j)
				IF (Rf(j) .lt. 0.0d0) THEN
					Rf(j) = 0.0d0
				END IF
                Sf(j+1)=snowfall(j+1) 
                Rf(j+1)=prec(j+1)-snowfall(j+1) 
				IF (Rf(j+1) .lt. 0.0d0) THEN
					Rf(j+1) = 0.0d0
				END IF
            ELSE
                Sf(j)=-999.0d0 
                Rf(j)=-999.0d0 
                Sf(j+1)=-999.0d0 
                Rf(j+1)=-999.0d0 
            END IF 
            
            
            hsnowst=hsnow(j)+dt*rhow/rhos*(0.5d0*(Sf(j)+Sf(j+1))-E)
            hsnowst=MAX(hsnowst,0.0d0)
                
            IF (0.5d0*(Tair(j)+Tair(j+1)) .lt. par(9)) THEN !GROWTH EQUATION (only Crank-Nicholson)
                dt=1.0d0
                ratio=0.0d0

                IF (rhos*hsnowst .gt. ((rhow-rhoi)*(hbice(j)+hsice(j))+(rhow-rhosl)*hslush(j))) THEN
                    csnow=rhos*hsnowst
                    cice=(rhow-rhoi)*(hsice(j)+hbice(j))
                    cslush=(rhow-rhosl)*hslush(j)

                    Dslushf=(csnow-(cice+cslush))/(rhos+rhow-rhosl)
                    Dslushf=MAX(Dslushf,0.0d0) ! probably this is not needed! Just a double-check
                ELSE
                    Dslushf=0.0d0
                END IF

                DslushRf=0.5d0*(Rf(j)+Rf(j+1))/ni 
                hslush(j+1)=hslush(j)+Dslushf + DslushRf 
                hslush(j+1)=MIN(hslush(j)+hsnowst,hslush(j+1))  !This means no H20 layer above snow
                hsnow(j+1)=hsnowst-Dslushf-DslushRf  
                hsnow(j+1)=MAX(hsnow(j+1),0.0d0)   ! probably this is not needed! Just a double-check

                IF (hslush(j+1) .gt. 0.0d0) THEN !GROWTH SNOW ICE
                    eps=1.0d0/((rhoi-rhos)*Lf) !rho1-rhos=rhoi-rhoi(i-ni) -> (only fraction of water in slush needs to be frozen)
                    fsin=(par(9)-Tair(j))/((hsnow(j)/ksnow)+1.0d0/par(10))             !par(10)=ca
                    fsinn=(par(9)-Tair(j+1))/((hsnow(j+1)/ksnow)+1.0d0/par(10))
                    IF ((fsin+fsinn).lt. 0.0d0) THEN !Shouldn't be activated
                        Dhsice = 0.0d0
                        !STOP
                    ELSE
                        Dhsice=0.50d0*dt*86400.0d0*eps*(fsin+fsinn) 
                    END IF
                    !Dhsice=MAX(Dhsice,0.0d0)  ! Accounts for poor-man solution: (fsin+fsinn).lt. 0.0d0 --> no fictitious melting

                    IF (Dhsice .eq. 0.0d0) THEN
                        ratio=hslush(j+1)/(hslush(j+1)*1D-6)  !num
                    ELSE
                        ratio=hslush(j+1)/(Dhsice)
                    END IF

                    Dhsice=MIN(Dhsice,hslush(j+1))
                    hsice(j+1)=hsice(j)+Dhsice
                ELSE
                    hsice(j+1)=hsice(j)
                END IF

                IF (ratio .gt. 1.0d0) THEN
                    hslush(j+1)=hslush(j+1)-Dhsice
                    hbice(j+1)=hbice(j)
                ELSE
                    hslush(j+1)=0.0d0
                END IF

                IF (hslush(j+1) .eq. 0.0d0) THEN !GROWTH BLACK ICE 
                    zeta=0.5d0*(1-ratio)*dt*86400.0d0/(rhoi*Lf)
                    theta=(par(9)-Tair(j))/(hbice(j)/ki+hsice(j)/ki+hsnow(j)/ksnow+1/par(10))
                    phi=par(9)-Tair(j+1)
                    psi=hsice(j+1)/ki+hsnow(j+1)/ksnow+1/par(10)
                    Bbi=psi*ki-hbice(j)-zeta*theta
                    Cbi=-hbice(j)*psi*ki-zeta*phi*ki-zeta*theta*psi*ki
                    hbice(j+1)=0.5d0*(-Bbi+SQRT(Bbi**2-4.0d0*Cbi))
                END IF

            ELSE
                fn=par(1)+par(2)*Tair(j)+par(5)*COS(2.d0*pi*(tt(j)-par(6)))
                fnn=par(1)+par(2)*Tair(j+1)+par(5)*COS(2.d0*pi*(tt(j)+ttt-par(6)))
				
				!Buoyancy
				IF (rhos*hsnowst .gt. ((rhow-rhoi)*(hbice(j)+hsice(j))+(rhow-rhosl)*hslush(j))) THEN
					csnow=rhos*hsnowst
					cice=(rhow-rhoi)*(hsice(j)+hbice(j))
					cslush=(rhow-rhosl)*hslush(j)
					Dslushf=(csnow-(cice+cslush))/(rhos+rhow-rhosl)
					Dslushf=MAX(Dslushf,0.0d0) !!
				ELSE
					Dslushf=0.0d0
                END IF
				
				DslushRf=0.5d0*(Rf(j)+Rf(j+1))/ni 
                hslush(j+1)=hslush(j)+Dslushf + DslushRf 
                hslush(j+1)=MIN(hslush(j)+hsnowst,hslush(j+1))  !This means no H2O layer above snow
                hsnow(j+1)=hsnowst-Dslushf - DslushRf  
                hsnow(j+1)=MAX(hsnow(j+1),0.0d0) 

                IF ((fn+fnn) .lt. 0.0d0) THEN
                    ! poor man solution
                    dt=1.0d0
                    hbice(j+1)=hbice(j)
                    hsice(j+1)=hsice(j)

                ELSEIF (hice(j) .gt. 0.0d0) THEN !MELTING EQUATION (only Crank-Nicholson)
                    dt=1.0d0
                    IF (hsnow(j+1) .gt. 0.0d0) THEN 
                        gammas=(par(11)*depth*rhow*cp)/(rhos*Lf)
                        !gammas=(depth*rhow*cp)/(rhos*Lf)
                        Rnn=0.5d0*dt*gammas*(fn+fnn)

                        IF (Rnn .lt. 0.0d0) THEN
                            WRITE(*,*) 'Non può essere, non deve mai entrare qui' !!!ModifiedOct2023 rimuovere la condizione se non entra mai!
                            STOP
                            Rnn=0.0d0
                        END IF
						
						IF (Rnn .gt. hsnow(j+1)) THEN 
                            sur=Rnn-hsnow(j+1) 
                            Rnn=hsnow(j+1)  
							hsnow(j+1)=0.d0 
							IF (hslush(j+1) .gt. 0.0d0 .and. sur .lt. hslush(j+1)) THEN
								hslush(j+1)=hslush(j+1)-sur 
								sur=0.0d0
								Dhicetot=0.0d0
							ELSEIF (hslush(j+1) .gt. 0.0d0 .and. sur .gt. hslush(j+1)) THEN
								sur=sur-hslush(j+1)
								Dhicetot=sur*rhos/rhoi
								hslush(j+1)=0.0d0
							END IF 
						ELSE
							sur=0.0d0
							Dhicetot=0.0d0
							Dslushm=(Rnn*rhos/rhow)/ni				
							hsnow(j+1)=hsnow(j+1)-Rnn 
							hsnow(j+1)=MAX(hsnow(j+1),0.0d0)
							Dslushm=MIN(Dslushm,hsnow(j+1))
							hslush(j+1)=hslush(j+1)+Dslushm 
						END IF
                        !!!
                        

                    ELSEIF (hsnow(j+1) == 0.0d0) THEN 
                        gammai=(par(11)*depth*rhow*cp)/(rhoi*Lf)
                        hicek=hice(j)-0.5d0*dt*gammai*(fn+fnn)

                        hice(j+1)=MAX(hicek,0.0d0)
                        Dhicetot=hice(j)-hice(j+1)
                    END IF

                    IF (hsice(j) .ge. Dhicetot) THEN
                        hsice(j+1)=hsice(j)-Dhicetot
                        hbice(j+1)=hbice(j)
                    ELSE
                        Dhicetot=Dhicetot-hsice(j)
                        hsice(j+1)=0.0d0
                        hbice(j+1)=hbice(j)-Dhicetot
                        IF (hbice(j+1)<epsilon(1.0d0) .and. hbice(j+1) .ne. 0.0d0) THEN
                            hbice(j+1)=0.0d0
                        END IF
                    END IF

                ELSE
                    hbice(j+1)=0.0d0
                    hsice(j+1)=0.0d0
                    hslush(j+1)=0.0d0
                    hsnow(j+1)=0.0d0
                END IF
            END IF
            hice(j+1)=hbice(j+1)+hsice(j+1)

        ELSE
            hice(j+1)=0.0d0
            hbice(j+1)=0.0d0
            hsice(j+1)=0.0d0
            hsnow(j+1)=0.0d0
            hslush(j+1)=0.0d0
        END IF

    END DO
    delta(n_tot) = DD

    RETURN
    END
    !-------------------------------------------------------------------------------
    !				INTEGRAZIONE NUMERICA
    !-------------------------------------------------------------------------------
    SUBROUTINE a2w(Ta, Tw, time, K, DD, pp)

    USE commondata

    IMPLICIT NONE

    REAL(KIND=8) :: lambda
    REAL(KIND=8), INTENT(OUT) :: K, DD, pp
    REAL(KIND=8), INTENT(IN) :: Ta, Tw, time


    IF (Tw>=Tmin) THEN
        DD=DEXP( -(Tw-Tmin)/par(4) );
    ELSE
        IF (version=='8') THEN
            DD=DEXP( (Tw-Tmin)/par(7) ) + DEXP( -Tw/par(8) );
        ELSE
            DD=1.0d0
        END IF
    END IF

    pp = par(1) + par(2)*Ta -  par(3)*Tw + par(5)*COS(2.d0*pi*(time-par(6))) ! Note that if version == '4' par(5)==0

    !!! lower bound for numerical stability
    !!lambda=pp/par(4) - par(3)
    !!IF (lambda .le. 0.0d0) THEN
    !!    IF (num_mod .eq. 'RK4' .and. DD .lt. -lambda/LIM) THEN
    !!        DD=-lambda/LIM
    !!    ELSEIF ((num_mod .eq. 'RK2' .or. num_mod .eq. 'EUL') .and. DD .lt. -lambda/LIM ) THEN
    !!        DD=-lambda/LIM
    !!    END IF
    !!END IF

    K = pp/DD

    RETURN
    END SUBROUTINE


    !-------------------------------------------------------------------------------
    !				ICE MODEL
    !-------------------------------------------------------------------------------

    SUBROUTINE ice(Ta, h, K)

    USE commondata

    IMPLICIT NONE

    REAL(KIND=8), INTENT(OUT) :: K
    REAL(KIND=8), INTENT(IN) :: Ta, h

    IF (version=='ICE0') THEN
        K = par(9) - Ta
    ELSEIF (version=='ICE' .or. version=='ICEM' .or. version=='ICES') THEN
        K = par(11)*(par(9)-Ta)/(h + par(10))
    END IF

    RETURN
    END SUBROUTINE

    !-------------------------------------------------------------------------------
    !				FUNZIONE OBIETTIVO
    !-------------------------------------------------------------------------------
    !SUBROUTINE funcobj(ind)
    SUBROUTINE funcobj(ind,indT,indI,indIb,indIs,indS)
    ! Subroutine per il calcolo della funzione obiettivo della simulazione
    ! Dichiarazione delle variabili
    USE commondata
    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8),INTENT(OUT) :: ind,indT,indI,indIb,indIs,indS
    REAL(KIND=8) :: TSS, mean_mod, TSS_mod, std_mod, covar_mod
    REAL(KIND=8) :: TSSI, meanI_mod, TSSI_mod, stdI_mod, covarI_mod
    REAL(KIND=8) :: TSSIb, TSSIs, TSSS
    REAL(KIND=8) :: tmp, tmpI, tmpIb, tmpIs, tmpS, max_, max_err

    max_err=-9999
    Twat_mod_agg=-999
    hice_mod_agg=-999
    hbice_mod_agg=-999
    hsice_mod_agg=-999
    hsnow_mod_agg=-999
    indT=-999
    indI=-999
    indIb=-999
    indIs=-999
    indS=-999
    DO i=1,n_dat            ! 1st year = warm-up period
        tmp=0.0d0
        DO j=I_inf(i,1),I_inf(i,2)
            tmp=tmp+Twat_mod(I_pos(j))
        END DO
        Twat_mod_agg(I_inf(i,3))=tmp/REAL(I_inf(i,2)-I_inf(i,1)+1)
    END DO
    DO i=1,n_datI           ! 1st year = warm-up period
        tmpI=0.0d0
        DO j=I_infI(i,1),I_infI(i,2)
            tmpI=tmpI+hice(I_posI(j))
        END DO
        hice_mod_agg(I_infI(i,3))=tmpI/REAL(I_infI(i,2)-I_infI(i,1)+1)
    END DO
    DO i=1,n_datIb           ! 1st year = warm-up period
        tmpIb=0.0d0
        DO j=I_infIb(i,1),I_infIb(i,2)
            tmpIb=tmpIb+hbice(I_posIb(j))
        END DO
        hbice_mod_agg(I_infIb(i,3))=tmpIb/REAL(I_infIb(i,2)-I_infIb(i,1)+1)
    END DO
    DO i=1,n_datIs           ! 1st year = warm-up period
        tmpIs=0.0d0
        DO j=I_infIs(i,1),I_infIs(i,2)
            tmpIs=tmpIs+hsice(I_posIs(j))
        END DO
        hsice_mod_agg(I_infIs(i,3))=tmpIs/REAL(I_infIs(i,2)-I_infIs(i,1)+1)
    END DO 
    DO i=1,n_datS           ! 1st year = warm-up period
        tmpS=0.0d0
        DO j=I_infS(i,1),I_infS(i,2)
            tmpS=tmpS+hsnow(I_posS(j))
        END DO
        hsnow_mod_agg(I_infS(i,3))=tmpS/REAL(I_infS(i,2)-I_infS(i,1)+1)
    END DO
    
    IF (fun_obj=='NSE') THEN
        TSS=0.d0
        TSSI=0.0d0
        TSSIb=0.0d0
        TSSIs=0.0d0
        TSSS=0.0d0
        DO i=1,n_dat
            TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
        END DO
        DO i=1,n_datI
            TSSI=TSSI+(hice_mod_agg(I_infI(i,3))-hice_obs_agg(I_infI(i,3)))**2
        END DO
        DO i=1,n_datIb
            TSSIb=TSSIb+(hbice_mod_agg(I_infIb(i,3))-hbice_obs_agg(I_infIb(i,3)))**2
        END DO
        DO i=1,n_datIs
            TSSIs=TSSIs+(hsice_mod_agg(I_infIs(i,3))-hsice_obs_agg(I_infIs(i,3)))**2
        END DO
        DO i=1,n_datS
            TSSS=TSSS+(hsnow_mod_agg(I_infS(i,3))-hsnow_obs_agg(I_infS(i,3)))**2
        END DO        
        indT=(1.d0-TSS/TSS_obs)
        indI=(1.d0-TSSI/TSSI_obs)
        indIb=(1.d0-TSSIb/TSSIb_obs)
        indIs=(1.d0-TSSIs/TSSIs_obs)
        indS=(1.d0-TSSS/TSSS_obs)
        IF (version=='ICES') THEN 
            IF (flag_icecal ==0) THEN
                ind=indT*beta+(1.0d0-beta)*indI
            ELSE IF (flag_icecal ==1) THEN
                ind=indT*beta+(1.0d0-beta)*(indIb+indIs)*1.0d0/2.0d0
            END IF
            !ind=indT*REAL(n_dat)/REAL(n_dat + n_datIb + n_datIs + n_datS)+ indIb*REAL(n_datIb)/REAL(n_dat + n_datIb + n_datIs + n_datS) + indIs*REAL(n_datIs)/REAL(n_dat + n_datIb + n_datIs + n_datS) +indS*REAL(n_datS)/REAL(n_dat + n_datIb + n_datIs + n_datS)
            !ind=indT*REAL(n_dat)/REAL(n_dat + n_datIb + n_datIs)+ indIb*REAL(n_datIb)/REAL(n_dat + n_datIb + n_datIs) + indIs*REAL(n_datIs)/REAL(n_dat + n_datIb + n_datIs)
        ELSE
            ind=indT
        END IF            
    ELSEIF  (fun_obj=='KGE') THEN
        mean_mod=0.0d0
        DO i=1,n_dat
            mean_mod=mean_mod+Twat_mod_agg(I_inf(i,3))
        END DO
        mean_mod=mean_mod/REAL(n_dat)

        covar_mod=0.0d0
        TSS_mod=0.0d0
        DO i=1,n_dat
            TSS_mod=TSS_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)**2
            covar_mod=covar_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)*(Twat_obs_agg(I_inf(i,3))-mean_obs)
        END DO
        std_mod=DSQRT(TSS_mod/REAL(n_dat-1))
        covar_mod=covar_mod/REAL(n_dat-1)
        ind=1.0d0-DSQRT((std_mod/std_obs-1.0d0)**2 + (mean_mod/mean_obs-1.0d0)**2 + (covar_mod/(std_mod*std_obs)-1.0d0)**2)
    ELSEIF  (fun_obj=='RMS') THEN
        TSS=0.d0
        DO i=1,n_dat
            TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
        END DO
        ind=-DSQRT(TSS/n_dat)        ! sign - --> the calibration procedure maximizes ind.
    ELSEIF  (fun_obj=='XXX') THEN
        WRITE(*,*) 'XXX non definito'
    ELSE
        WRITE(*,*) 'Errore nella scelta della f. obiettivo'
    END IF

    RETURN
    END
    !-------------------------------------------------------------------------------
    !				AGGREGATION
    !-------------------------------------------------------------------------------
    SUBROUTINE aggregation

    ! Dichiarazione delle variabili
    USE commondata
    IMPLICIT NONE

    INTEGER :: i, j, k, status, pp
    INTEGER :: n_inf, n_pos, n_units, n_days, count, pos_tmp, n_infI,n_infIb, n_infIs, n_infS, n_posI, n_posIb, n_posIs, n_posS, n_unitsI, n_unitsIb, n_unitsIs, n_unitsS
    INTEGER :: month, month_curr
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: A
    INTEGER, DIMENSION(:), ALLOCATABLE :: B
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: AI, AIb, AIs, AS
    INTEGER, DIMENSION(:), ALLOCATABLE :: BI, BIb, BIs, BS
    REAL(KIND=8) :: tmp

    pp=LEN_TRIM(time_res)
    IF (pp==2) THEN
        unit=time_res(2:)
        READ(time_res,'(i1)') qty
    ELSEIF (pp==3) THEN
        unit=time_res(3:)
        READ(time_res,'(i2)') qty
    END IF

    ALLOCATE(I_pos(n_tot),stat=status)
    ALLOCATE(I_posI(n_tot),stat=status)
    ALLOCATE(I_posIb(n_tot),stat=status)
    ALLOCATE(I_posIs(n_tot),stat=status)
    ALLOCATE(I_posS(n_tot),stat=status)
    I_pos=-999
    I_posI=-999
    I_posIb=-999
    I_posIs=-999
    I_posS=-999
    Twat_obs_agg=-999
    hice_obs_agg=-999
    hbice_obs_agg=-999
    hsice_obs_agg=-999
    hsnow_obs_agg=-999
    n_inf=1
    n_pos=1
    n_infI=1
    n_infIb=1
    n_infIs=1
    n_infS=1
    n_posI=1
    n_posIb=1
    n_posIs=1
    n_posS=1
    IF (time_res=='1d') THEN    ! daily resolution
        n_units=n_tot-365
        n_unitsI=n_tot-365
        n_unitsIb=n_tot-365
        n_unitsIs=n_tot-365
        n_unitsS=n_tot-365
        ALLOCATE(I_inf(n_units,3))
        I_inf=-999
        DO i=366,n_tot    ! 1st year = warm-up period
            IF (Twat_obs(i) .ne. -999) THEN
                I_inf(n_inf,2)=n_pos
                I_inf(n_inf,3)=i
                I_pos(n_pos)=i
                Twat_obs_agg(I_inf(n_inf,3))=Twat_obs(i)
                n_inf=n_inf+1
                n_pos=n_pos+1
            END IF
        END DO
        ALLOCATE(I_infI(n_unitsI,3))
        I_infI=-999
        DO i=366,n_tot    ! 1st year = warm-up period
            IF (hice_obs(i) .ne. -999) THEN
                I_infI(n_infI,2)=n_posI
                I_infI(n_infI,3)=i
                I_posI(n_posI)=i
                hice_obs_agg(I_infI(n_infI,3))=hice_obs(i)
                n_infI=n_infI+1
                n_posI=n_posI+1
            END IF
        END DO
        ALLOCATE(I_infIb(n_unitsIb,3))
        I_infIb=-999
        DO i=366,n_tot    ! 1st year = warm-up period
            IF (hbice_obs(i) .ne. -999) THEN
                I_infIb(n_infIb,2)=n_posIb
                I_infIb(n_infIb,3)=i
                I_posIb(n_posIb)=i
                hbice_obs_agg(I_infIb(n_infIb,3))=hbice_obs(i)
                n_infIb=n_infIb+1
                n_posIb=n_posIb+1
            END IF
        END DO
        ALLOCATE(I_infIs(n_unitsIs,3))
        I_infIs=-999
        DO i=366,n_tot    ! 1st year = warm-up period
            IF (hsice_obs(i) .ne. -999) THEN
                I_infIs(n_infIs,2)=n_posIs
                I_infIs(n_infIs,3)=i
                I_posIs(n_posIs)=i
                hsice_obs_agg(I_infIs(n_infIs,3))=hsice_obs(i)
                n_infIs=n_infIs+1
                n_posIs=n_posIs+1
            END IF
        END DO
        ALLOCATE(I_infS(n_unitsS,3))
        I_infS=-999
        DO i=366,n_tot    ! 1st year = warm-up period
            IF (hsnow_obs(i) .ne. -999) THEN
                I_infS(n_infS,2)=n_posS
                I_infS(n_infS,3)=i
                I_posS(n_posS)=i
                hsnow_obs_agg(I_infS(n_infS,3))=hsnow_obs(i)
                n_infS=n_infS+1
                n_posS=n_posS+1
            END IF
        END DO
    ELSEIF (unit=='w') THEN      ! weekly resolution
        n_days=qty*7
        n_units=CEILING(REAL(n_tot-365)/REAL(n_days))
        ALLOCATE(I_inf(n_units,3))
        I_inf=-999
        DO i=366,n_tot,n_days
            tmp=0.0d0
            count=0
            pos_tmp=i+CEILING(0.5*n_days)-1
            DO j=0,n_days-1
                k=i+j
                IF (k .gt. n_tot) THEN
                    EXIT
                END IF
                IF (Twat_obs(k) .ne. -999) THEN
                    tmp=tmp+Twat_obs(k)
                    I_pos(n_pos)=k
                    n_pos=n_pos+1
                    count=count+1
                END IF
            END DO
            IF (count .ge. n_days*prc) THEN
                I_inf(n_inf,2)=n_pos-1
                I_inf(n_inf,3)=pos_tmp
                Twat_obs_agg(I_inf(n_inf,3))=tmp/count
                n_inf=n_inf+1
            ELSE
                I_pos(n_pos-count:n_pos)=-999
                n_pos=n_pos-count
            END IF
        END DO
    ELSEIF (unit=='m') THEN      ! monthly resolution
        n_units=CEILING(n_tot/(30.5))
        ALLOCATE(I_inf(n_units,3))
        I_inf=-999
        n_days=0
        month_curr=-999
        count=0
        DO i=366,n_tot
            month=date(i,2)
            IF (month .ne. month_curr) THEN
                IF (count .ge. n_days*prc .and. i .ne. 366) THEN
                    I_inf(n_inf,2)=n_pos-1
                    I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
                    Twat_obs_agg(I_inf(n_inf,3))=tmp/count
                    n_inf=n_inf+1
                ELSE
                    I_pos(n_pos-count:n_pos)=-999
                    n_pos=n_pos-count
                END IF
                month_curr=month
                count=0
                n_days=1
                tmp=0.0d0
            ELSE
                n_days=n_days+1
            END IF
            IF (Twat_obs(i) .ne. -999) THEN
                tmp=tmp+Twat_obs(i)
                I_pos(n_pos)=i
                n_pos=n_pos+1
                count=count+1
            END IF
        END DO
        ! Last month
        IF (count .ge. n_days*prc) THEN
            I_inf(n_inf,2)=n_pos-1
            I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
            Twat_obs_agg(I_inf(n_inf,3))=tmp/count
            n_inf=n_inf+1
        ELSE
            I_pos(n_pos-count:n_pos)=-999
            n_pos=n_pos-count
        END IF
    ELSE
        WRITE(*,*) 'Error: variable time_res'
    END IF

    n_dat=n_inf-1
    n_pos=n_pos-1
    n_datI=n_infI-1
    n_datIb=n_infIb-1
    n_datIs=n_infIs-1
    n_datS=n_infS-1
    n_posI=n_posI-1
    n_posIb=n_posIb-1
    n_posIs=n_posIs-1
    n_posS=n_posS-1

    I_inf(1,1)=1
    I_infI(1,1)=1
    I_infIb(1,1)=1
    I_infIs(1,1)=1
    I_infS(1,1)=1
    I_inf(2:n_dat,1)=I_inf(1:n_dat-1,2)+1
    I_infI(2:n_datI,1)=I_infI(1:n_datI-1,2)+1
    I_infIb(2:n_datIb,1)=I_infIb(1:n_datIb-1,2)+1
    I_infIs(2:n_datIs,1)=I_infIs(1:n_datIs-1,2)+1
    I_infS(2:n_datS,1)=I_infS(1:n_datS-1,2)+1
    ALLOCATE(A(n_units,3),B(n_tot))
    ALLOCATE(AI(n_unitsI,3),BI(n_tot))
    ALLOCATE(AIb(n_unitsIb,3),BIb(n_tot))
    ALLOCATE(AIs(n_unitsIs,3),BIs(n_tot))
    ALLOCATE(AS(n_unitsS,3),BS(n_tot))
    A=I_inf;
    AI=I_infI;
    AIb=I_infIb;
    AIs=I_infIs;
    AS=I_infS;
    B=I_pos;
    BI=I_posI;
    BIb=I_posIb;
    BIs=I_posIs;
    BS=I_posS;
    DEALLOCATE(I_inf,I_pos)
    DEALLOCATE(I_infI,I_posI)
    DEALLOCATE(I_infIb,I_posIb)
    DEALLOCATE(I_infIs,I_posIs)
    DEALLOCATE(I_infS,I_posS)
    ALLOCATE(I_inf(n_dat,3),I_pos(n_pos))
    ALLOCATE(I_infI(n_datI,3),I_posI(n_posI))
    ALLOCATE(I_infIb(n_datIb,3),I_posIb(n_posIb))
    ALLOCATE(I_infIs(n_datIs,3),I_posIs(n_posIs))
    ALLOCATE(I_infS(n_datS,3),I_posS(n_posS))
    I_inf=A(1:n_dat,:)
    I_infI=AI(1:n_datI,:)
    I_infIb=AIb(1:n_datIb,:)
    I_infIs=AIs(1:n_datIs,:)
    I_infS=AS(1:n_datS,:)
    I_pos=B(1:n_pos)
    I_posI=BI(1:n_posI)
    I_posIb=BIb(1:n_posIb)
    I_posIs=BIs(1:n_posIs)
    I_posS=BS(1:n_posS)
    DEALLOCATE(A,B)
    DEALLOCATE(AI,BI)
    DEALLOCATE(AIb,BIb)
    DEALLOCATE(AIs,BIs)
    DEALLOCATE(AS,BS)


    RETURN
    END
    !-------------------------------------------------------------------------------
    !				STATIS
    !-------------------------------------------------------------------------------
    SUBROUTINE statis
    ! Subroutine per il calcolo di media, somma degli scarti quadratici e std dei dati

    ! Dichiarazione delle variabili
    USE commondata
    IMPLICIT NONE

    INTEGER :: i, k, d, status

    mean_obs=0.d0
    meanI_obs=0.d0
    meanIb_obs=0.d0
    meanIs_obs=0.d0
    meanS_obs=0.d0
    TSS_obs=0.d0
    TSSI_obs=0.d0
    TSSIb_obs=0.d0
    TSSIs_obs=0.d0
    TSSS_obs=0.d0
    DO i=1,n_dat
        mean_obs=mean_obs+Twat_obs_agg(I_inf(i,3))
    END DO
    DO i=1,n_datI
        meanI_obs=meanI_obs+hice_obs_agg(I_infI(i,3))
    END DO
    DO i=1,n_datIb
        meanIb_obs=meanIb_obs+hbice_obs_agg(I_infIb(i,3))
    END DO
    DO i=1,n_datIs
        meanIs_obs=meanIs_obs+hsice_obs_agg(I_infIs(i,3))
    END DO
    DO i=1,n_datS
        meanS_obs=meanS_obs+hsnow_obs_agg(I_infS(i,3))
    END DO
    mean_obs=mean_obs/REAL(n_dat)
    meanI_obs=meanI_obs/REAL(n_datI)
    meanIb_obs=meanIb_obs/REAL(n_datIb)
    meanIs_obs=meanIs_obs/REAL(n_datIs)
    meanS_obs=meanS_obs/REAL(n_datS)

    DO i=1,n_dat
        TSS_obs=TSS_obs+(Twat_obs_agg(I_inf(i,3))-mean_obs)**2
    END DO
    DO i=1,n_datI
        TSSI_obs=TSSI_obs+(hice_obs_agg(I_infI(i,3))-meanI_obs)**2
    END DO
    DO i=1,n_datIb
        TSSIb_obs=TSSIb_obs+(hbice_obs_agg(I_infIb(i,3))-meanIb_obs)**2
    END DO
    DO i=1,n_datIs
        TSSIs_obs=TSSIs_obs+(hsice_obs_agg(I_infIs(i,3))-meanIs_obs)**2
    END DO
    DO i=1,n_datS
        TSSS_obs=TSSS_obs+(hsnow_obs_agg(I_infS(i,3))-meanS_obs)**2
    END DO

    std_obs=DSQRT(TSS_obs/REAL(n_dat-1))
    stdI_obs=DSQRT(TSSI_obs/REAL(n_datI-1))
    !stdS_obs=DSQRT(TSSS_obs/REAL(n_datS-1))

    RETURN
    END

    !-------------------------------------------------------------------------------
    !				FORWARD
    !-------------------------------------------------------------------------------
    SUBROUTINE forward
    USE commondata
    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: ei_check, ei,eiT,eiI,eiIb,eiIs,eiS

    par=par_best		! uso miglior set di paramteri

    CALL model

    ! Controllo: calcolo dell'indice di efficienza
    !CALL funcobj(ei_check)
    CALL funcobj(ei_check,eiT,eiI,eiIb,eiIs,eiS)
    WRITE(*,*) 'Indice efficienza calibrazione', ei_check

    IF (ABS(ei_check - finalfit) .gt. 0.0001) THEN
        WRITE(*,*) 'Errore efficienza in forward'
        WRITE(*,*) ei_check, finalfit
        PAUSE
    ELSE
        WRITE(*,*) 'Controllo superato'
    END IF

    WRITE(11,'(<n_par>(f10.6,1x))') (par_best(i),i=1,n_par)
    WRITE(11,'(f10.6)') ei_check

    OPEN(UNIT=12,FILE=TRIM(folder)//'/2_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'c_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
    DO i=1,n_tot
        !WRITE(12,1004) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i)
        ! WRITE(12,1005) (date(i,j),j=1,3),Tair(i), Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i),hbice_obs(i),hbice(i),hsice_obs(i),hsice(i),hsnow_obs(i),hsnow(i),prec(i) 
        WRITE(12,1006) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i),hbice_obs(i),hbice(i),hsice_obs(i),hsice(i),hsnow_obs(i),hsnow(i),prec(i),snowfall(i),Sf(i) 

    END DO
    CLOSE(12)

    OPEN(UNIT=14,FILE=TRIM(folder)//'/4_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'c_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
    DO i=1,n_tot
        WRITE(14,*) delta(i), hbice(i)
    END DO
    CLOSE(14)

    CALL read_validation

    IF (n_tot .lt. 365) THEN
        ei=-999
        GO TO 200
    END IF

    ! aggregate calibration data on the basis of time_res
    CALL aggregation

    CALL statis
    WRITE(*,*) 'mean, TSS and standard deviation (validation)'
    WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)
    WRITE(*,*)  SNGL(meanI_obs),SNGL(TSSI_obs),SNGL(stdI_obs)
    !WRITE(*,*)  SNGL(meanS_obs),SNGL(TSSS_obs),SNGL(stdS_obs)

    CALL model
    !CALL funcobj(ei)
    CALL funcobj(ei,eiT,eiI,eiIb,eiIs,eiS)
    WRITE(11,'(f10.6)') ei

    CLOSE(11)

    OPEN(UNIT=13,FILE=TRIM(folder)//'/3_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'v_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
    DO i=1,n_tot
        !WRITE(13,1004) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i)
        !WRITE(13,1005) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i),hbice_obs(i),hbice(i),hsice_obs(i),hsice(i),hsnow_obs(i),hsnow(i),prec(i) 
        WRITE(13,1006) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),hice_obs(i),hice(i),hbice_obs(i),hbice(i),hsice_obs(i),hsice(i),hsnow_obs(i),hsnow(i),prec(i),snowfall(i),Sf(i) 
    END DO
    CLOSE(13)

    OPEN(UNIT=15,FILE=TRIM(folder)//'/5_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'v_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
    DO i=1,n_tot
        WRITE(15,*) delta(i), hbice(i)
    END DO
    CLOSE(15)

1004 FORMAT(i4,1x,i4,1x,i4,1x,5(1x,f10.5))
1005 FORMAT(i4,1x,i4,1x,i4,1x,12(1x,f10.5))
1006 FORMAT(i4,1x,i4,1x,i4,1x,14(1x,f10.5))     

200 RETURN
    END

    !-------------------------------------------------------------------------------
    !				BEST
    !-------------------------------------------------------------------------------
    SUBROUTINE best(fit,part,foptim)
    USE commondata
    IMPLICIT NONE

    INTEGER,INTENT(OUT)::part
    INTEGER:: k
    REAL(KIND=8),INTENT(IN),DIMENSION(n_particles):: fit
    REAL(KIND=8),INTENT(OUT):: foptim

    foptim=-1e30				! valore molto piccolo
    DO k=1,n_particles
        IF(fit(k).gt.foptim) then
            foptim=fit(k)
            part=k
        END IF
    END DO

    RETURN
    END

    !-------------------------------------------------------------------------------
    !				LEAP YEAR
    !-------------------------------------------------------------------------------
    SUBROUTINE leap_year(Y,I)
    USE commondata
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: Y
    INTEGER,INTENT(OUT) :: I

    IF(MOD(Y,100).NE.0.AND.MOD(Y,4).EQ.0) THEN
        I=1
    ELSEIF(MOD(Y,400).EQ.0) THEN
        I=1
    ELSE
        I=0
    END IF

    RETURN
    END


    !-------------------------------------------------------------------------------
    !				LINEAR REGRESSION
    !-------------------------------------------------------------------------------
    SUBROUTINE linreg(n,X,Y,m,b,r2)
    IMPLICIT NONE

    ! adapted from http://www.pgccphy.net/Linreg/linreg_f90.txt
    ! Dr. David G. Simpson

    INTEGER :: i
    INTEGER, INTENT(IN):: n
    REAL(KIND=8), INTENT(IN):: X(n),Y(n)
    REAL(KIND=8), INTENT(OUT):: m, b, r2
    REAL(KIND=8):: sumx, sumx2, sumxy, sumy, sumy2

    sumx=0.0d0; sumx2=0.0d0; sumxy=0.0d0; sumy=0.0d0; sumy2=0.0d0;
    DO i=1,n
        sumx  = sumx + x(i)                     ! compute sum of x
        sumx2 = sumx2 + x(i) * x(i)             ! compute sum of x**2
        sumxy = sumxy + x(i) * y(i)             ! compute sum of x * y
        sumy  = sumy + y(i)                     ! compute sum of y
        sumy2 = sumy2 + y(i) * y(i)             ! compute sum of y**2
    END DO

    m = (n * sumxy  -  sumx * sumy) / (n * sumx2 - sumx**2)                          ! compute slope
    b = (sumy * sumx2  -  sumx * sumxy) / (n * sumx2  -  sumx**2)                    ! compute y-intercept
    r2 = (sumxy - sumx * sumy / n) /                                     &            ! compute correlation coefficient
        sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))
    r2=r2**2
    RETURN

    END
