PROGRAM MAIN
IMPLICIT NONE

integer, PARAMETER :: jmax=5000
REAL*8, dimension(jmax) :: r(jmax), rr(jmax), rho_dm(jmax), M_DM(jmax), M_DM_EXACT(jmax)
REAL*8, DIMENSION(jmax) :: vol_r(jmax), rho_g(jmax), rho_g_exact(jmax), M_S(jmax)
REAL*8, DIMENSION(jmax) :: bf_stars(jmax), M_GS(jmax), t_r(jmax), grad_T(jmax)

REAL*8, DIMENSION(jmax) :: rho_trash(jmax), M_trash(jmax)

REAL*8 :: pi, msol, kpc, years, mu, boltz, G, mp, zfesol, zfesn, snu , temp_K, eps_bar_frac, brayonic_fraction
REAL*8 :: rho_dm0, r_s, r_min, r_max, r_vir, b,a, time_start_debug, time_end_debug, time_scale_debug
REAL*8 :: x, rho_0, mBCG, SPL, SPR, rho_0_control, bf_spl, bf_spr, bf_c_control, rho_c, gyr, M_Fe_100kpc
REAL*8 :: fc, x500, r500, yy,y, M_fe_tot, output_time_control, delta_t_debug, SN_unit, M_Fe_100kpc_initial
REAL*8 :: M_Fe_100kpc_final, M_Fe_initial, zfeout, initial_time, kpc_check_mass, tau_diffusion
INTEGER :: i, random_approach_sn, index_100kpc, start_with_zero_iron, stop_debug_mass_index
CHARACTER(LEN=100) :: Format

LOGICAL :: COND_1

! Declaration of variables for the second part of the project

REAL*8 :: bg_zfe, v_turb, l_turb, param_D(1), C_param, dt, time_wanted, time_passed, age_of_supernovae
REAL*8, DIMENSION(jmax) :: grad_zfe(jmax), rho_g_mean(jmax), Z_Fe(jmax), rho_fe(jmax), n_e_paper(jmax), n_tot_paper(jmax)
REAL*8, DIMENSION(jmax) :: rho_paper(jmax), int_mass_sn(jmax), rho_stellar(jmax), source_term_rho(jmax)
REAL*8, DIMENSION(jmax) :: t_paper(jmax), t_profile(jmax), M_GSDM_T(jmax), rho_gas(jmax), M_fe(jmax)
INTEGER ::  timesteps, k_time_step


LOGICAL ::  end_of_time


brayonic_fraction = 0.16
fc = 1.0
pi=3.14159265359
msol = 1.989d33
kpc = 3.084e21
years=3.156d7
mu=0.61
boltz=1.38066e-16
G=6.6720e-8
mp=1.67265e-24
zfesol=1.8e-3
zfesn=0.744/1.4
snu=0.5

gyr = 3.156d16
rho_0 = 4.d-26

start_with_zero_iron = 0
! clear the terminal output with system library
call system("clear")
WRITE(*,*) "  _____________________________________________________________________"
WRITE(*,*) " |                                                                     | "
WRITE(*,*) " |                    IRON DIFFUSION IN CLUSTER                        |"
WRITE(*,*) " |_____________________________________________________________________|"
WRITE(*,*) 

r_s = 435.7 * kpc    ! radius of the cluster
rho_dm0 = 7.35d-26   ! density of dark matter g/cm^3
r_vir = 2800*kpc     ! virial radius of the cluster
r_min = 0.*kpc       ! minimum radius of the grid
r500 = r_vir/2
r_max = 3000 * kpc   ! maximum radius of cluster
temp_K = 8.9d7
b=8.*pi*G*mu*mp*rho_dm0*r_s**2./(27.*boltz*temp_K)
stop_debug_mass_index = 0
! GRID BUILDING !
DO i = 1, jmax
   r(i) = r_min + (i-1)*r_max/(jmax-1)
END DO
! Now let's make the grid so that it's denser in the center

! rmax_temp = r(jmax)**1.3

! write(*,*) "DEBUGGGGGGGGG", r(jmax), rmax_temp
! DO i = 1, jmax
!    ! r(i) = 1/(20*(r(i)/r_s)*(1+r(i)/r_s)**2)
!    ! r(i) = 1/(r(i)/rmax_temp/r_max)ù
!    r(i) = r(i)**1.3
!    r(i) = r(i)/rmax_temp*r_max
! END DO

FORMAT = "(e15.7,e15.7)"
OPEN(20, FILE="grid.dat")

DO i = 1, jmax-1
   rr(i)=r(i)+0.5*(r(i+1)-r(i))
   WRITE(20,FORMAT) r(i)/kpc, rr(i)/kpc
END DO
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
WRITE(20,FORMAT) r(jmax)/kpc, rr(jmax)/kpc
CLOSE(20)

! DARK MATTER DENSITY PROFILE !
DO i = 1, jmax
   x = rr(i)/r_s
   rho_dm(i) = rho_dm0 / (x*(1.+x)**2) ! Density of dark matter defined on the rr(i) grid
END DO

vol_r(1) = 4./3.*pi*r(1)**3
DO i = 2, jmax-1
   vol_r(i) = 4./3.*pi*(r(i+1)**3-r(i)**3) ! Volume between r(i+1) and r(i) ,DEFINED AT rr(i) or rr(i+/-1/2)
END DO
vol_r(jmax) = 4./3.*pi*(r(jmax)**3-r(jmax-1)**3)
                                          
! DM MASS INTEGRATION !
M_DM(1) = rho_dm0*vol_r(1)
M_DM_EXACT(1) = rho_dm0*vol_r(1)
x=rr(1)/r_s
DO i = 2, jmax
x = rr(i)/r_s
M_DM(i) = M_DM(i-1) + rho_dm(i) * vol_r(i)
M_DM_EXACT(i) = 4.*pi*rho_dm0*r_s**3*(log(1.+rr(i)/r_s)-rr(i)/(rr(i)+r_s))
END DO
                                 
FORMAT = "(e11.4,e11.4)"
OPEN(20, FILE="DM_EXACT_VS_NUMERIC.dat")
WRITE(20,*) "M_DM ", "M_DM_EXACT"
DO i = 1, jmax
   WRITE(20,FORMAT) M_DM(i)/msol, M_DM_EXACT(i)/msol
END DO
CLOSE(20)


FORMAT = "(3(e15.7))"
OPEN(20, FILE="rho_gas_only_DM.dat")
REWIND(20)
rho_g(1) = exp(log(rho_0))
rho_g_exact(1) = rho_0
WRITE(20, *) "rho_gas ", "rho_gas_exact ", "rho_dm"
WRITE(20, FORMAT) rho_g(1), rho_g_exact(1), rho_dm(1)
M_GS(1) = rho_g(1)*vol_r(1)
DO i = 2, jmax-1
   rho_g(i) = log(rho_g(i-1))-(G*M_DM(i)/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_g(i) = exp(rho_g(i))
   rho_g_exact(i) = rho_0*exp(-27.*b/2.)*(1.+r(i)/r_s)**(27.*b/(2.*r(i)/r_s))
   M_GS(i) = M_GS(i-1) + rho_g(i)*vol_r(i)
   WRITE(20, FORMAT) rho_g(i), rho_g_exact(i), rho_dm(i)
END DO
rho_g(jmax) = exp(log(rho_g(jmax-1))-(G*M_DM(jmax)/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
rho_g_exact(jmax) = rho_0*exp(-27*b/2)*(1+rr(jmax)/r_s)**(27*b/(2*rr(jmax)/r_s))
M_GS(jmax) = M_GS(jmax-1) + rho_g(jmax) * vol_r(jmax)
WRITE(20, FORMAT) rho_g(jmax), rho_g_exact(jmax), rho_dm(jmax)

CLOSE(20)


! Some useful variables for the analysis
! Check halving of the iron mass at XXX kpc
! WRITE(*,*) " At which radius do you want to check the halving of the iron mass?  [kpc]"
! READ(*,*) kpc_check_mass
kpc_check_mass = 100
DO i = 1, jmax
   IF(rr(i) .GE. kpc_check_mass*kpc) THEN
      index_100kpc = i
      EXIT
   END IF
END DO


mBCG = 10.**12 * msol ! Total mass of the central galaxy
a = 12.*kpc/(1.+sqrt(2.))
DO i = 1, jmax
   M_S(i) = mBCG*(r(i))**2/((r(i)+a)**2) ! Stars mass defined on the rr(i) grid
   bf_stars(i) = (M_S(i)+M_GS(i))/(M_DM(i)+M_S(i)+M_GS(i))
END DO
CLOSE(20)

WRITE(*,*) "  _____________________________________________________________________"
WRITE(*,*) " |                                                                     | "
WRITE(*,*) " |  ANALYSIS WITH JUST DARK MATTER WITH CONSTANT TEMPERATURE PROFILE   |"
WRITE(*,*) " |_____________________________________________________________________|"
WRITE(*,*) 

FORMAT = "(3(e15.7))"

!!! Bisection method to find the right density for the gas that makes the baryonic fraction to be 0.16

eps_bar_frac = 0.00001 ! epsilon for stopping the bisection method
COND_1 = .FALSE.
rho_0_control = 4.e-26
SPL = 1.E-22 ! Left guess for gas density, "very big"
SPR = 1.E-28 ! Right guess for gas density, "very small"

DO WHILE(COND_1 .EQV. .FALSE.)  ! Bisection method, while until the condition is true (right value for the gas density is found)
rho_trash(1) = exp(log(SPL)) 
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*M_DM(i)/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*M_DM(jmax)/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

! NOW I CALCULATE THE BAR FRACTION WITH SPL DENSITY, and save it in temp storage bf_spl.
! The virial radius of 2.8 Mpc is around index 4667 of the grid. But we et it to JMAX for semplicity
bf_spl = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

! Now I do the same for the right guess SPR
rho_trash(1) = exp(log(SPR)) 
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*M_DM(i)/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*M_DM(jmax)/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_spr = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

rho_c = (SPL + SPR)/2 ! Middle value between SPL and SPR

rho_trash(1) = exp(log(rho_c))
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*M_DM(i)/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*M_DM(jmax)/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_c_control =(M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

! NOW IF F(C) > 0 then SPL = C, AND SPR = OLD SPR
!		IF F(C) < 0 then SPL = OLD SPL, SPR = C
! UNTIL abs(F(c)) < epsilon 

IF(0.16 - bf_c_control .le. 0.) THEN
   SPL = rho_c
ELSE 
   SPR = rho_c
END IF

IF(abs(0.16 - bf_c_control) .le. eps_bar_frac) THEN
   COND_1 = .true.
   rho_trash(1) = rho_c
END IF


END DO ! End of the bisection method, now we have the right value for the gas density

M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*M_DM(i)/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*M_DM(jmax)/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

FORMAT = "(7(e15.7))"
OPEN(30, file="data_only_gas.dat")
WRITE(30,*) "r"," rr"," M_DM"," M_S"," M_G"," rho_dm", " rho_gas"
DO i = 1, jmax
   WRITE(30,FORMAT) r(i)/kpc, rr(i)/kpc, M_DM(i)/msol, M_S(i)/msol, M_trash(i)/msol, rho_dm(i), rho_trash(i)
END DO
CLOSE(30)
WRITE(*,*) "The right density is: ", rho_trash(1)
WRITE(*,*) "Baryonic fraction:     ", bf_c_control
WRITE(*,*) ""

WRITE(*,*) "  _____________________________________________________________________"
WRITE(*,*) " |                                                                     | "
WRITE(*,*) " |    ANALYSIS WITH STARS + DM, WITH CONSTANT TEMPERATURE PROFILE      | "
WRITE(*,*) " |_____________________________________________________________________| "
WRITE(*,*) 

COND_1 = .FALSE.
rho_0_control = 4.e-26
SPL = 1.E-22
SPR = 1.E-28

DO WHILE(COND_1 .EQV. .FALSE.)
rho_trash(1) = exp(log(SPL))
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_spl = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

rho_trash(1) = exp(log(SPR)) 
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_spr = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

rho_c = (SPL + SPR)/2

rho_trash(1) = exp(log(rho_c)) 
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/r(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc)
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_c_control =(M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

IF(0.170 - bf_c_control .le. 0.) THEN
   SPL = rho_c
ELSE 
   SPR = rho_c
END IF

IF(abs(0.170 - bf_c_control) .le. eps_bar_frac) THEN
   COND_1 = .true.
   rho_trash(1) = rho_c
END IF
END DO

FORMAT = "(7(e15.7))"
OPEN(30, file="data_DM_plus_STARS_T_constant.dat")
WRITE(30,*) "r"," rr"," M_DM"," M_S"," M_GSDM"," rho_dm", " rho_gas_DM_plus_STARS"
DO i = 1, jmax
   WRITE(30,FORMAT) r(i)/kpc, rr(i)/kpc, M_DM(i)/msol, M_S(i)/msol, M_trash(i)/msol, rho_dm(i), rho_trash(i)
END DO
CLOSE(30)
WRITE(*,*) "The right density is: ", rho_trash(1)
WRITE(*,*) "Baryonic fraction:     ", bf_c_control
WRITE(*,*) ""

WRITE(*,*) "  _____________________________________________________________________  "
WRITE(*,*) " |                                                                     | "
WRITE(*,*) " |    ANALYSIS WITH STARS + DM, WITH VARIABLE TEMPERATURE PROFILE      | "
WRITE(*,*) " |_____________________________________________________________________| "
WRITE(*,*) 

! NOW WE SHOULD DO THE ANALYSIS AGAIN USING A TEMPERATURE PROFILE NOT CONSTANT, SO THIS
! WILL ADD A GRADIENT TERM INTO THE EQUATION OF HYDROSTATIC EQUILIBRIUM

COND_1 = .FALSE.
rho_0_control = 4.e-25
SPL = 1.E-22
SPR = 1.E-28

! We first have to create the temperature profile in function of radius
OPEN(20, FILE='temperature.dat')
DO i = 1, jmax
   y=rr(i)/r500
   yy=y/0.045
   x500 = rr(i)/(1.4*1000*kpc)
   t_r(i) = temp_k*1.35*(yy**1.9+0.45)/(yy**1.9+1.)* &   !! this is for Perseus !!
        1./(1.+(y/0.6)**2)**0.45
   WRITE(20,"(e15.7)") t_r(i)
END DO
CLOSE(20)

grad_T(1) = 0.
DO i=2, jmax
   grad_T(i) = log(t_r(i))-log(t_r(i-1))
END DO

DO WHILE(COND_1 .EQV. .FALSE.)
rho_trash(1) = exp(log(SPL))
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   temp_K = 0.5*(t_r(i)+t_r(i-1))
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i)-grad_T(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc &
                      -grad_T(jmax))
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_spl = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

rho_trash(1) = exp(log(SPR))
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   temp_K = 0.5*(t_r(i)+t_r(i-1))
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/rr(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i)-grad_T(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
temp_K = 0.5*(t_r(jmax)+t_r(jmax-1))
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc &
                      -grad_T(i))
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_spr = (M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax-1)+M_DM(jmax-1)+M_S(jmax-1))

rho_c = (SPL + SPR)/2

rho_trash(1) = exp(log(rho_c))
M_trash(1) = rho_trash(1)*vol_r(1)
DO i = 2, jmax-1
   temp_K = 0.5*(t_r(i)+t_r(i-1))
   rho_trash(i) = log(rho_trash(i-1))-(G*(M_DM(i)+M_S(i))/r(i)**2)*mu*mp*(rr(i+1)-rr(i))/(boltz*temp_K)*fc
   rho_trash(i) = exp(rho_trash(i)-grad_T(i))
   M_trash(i) = M_trash(i-1) + rho_trash(i)*vol_r(i)
END DO
temp_K = 0.5*(t_r(jmax)+t_r(jmax-1))
rho_trash(jmax) = exp(log(rho_trash(jmax-1))-(G*(M_DM(jmax)+M_S(jmax))/r(jmax)**2)*mu*mp*(rr(jmax)-rr(jmax-1))/(boltz*temp_K)*fc & 
                      -grad_T(i))
M_trash(jmax) = M_trash(jmax-1) + rho_trash(jmax) * vol_r(jmax)

bf_c_control =(M_trash(jmax-1)+M_S(jmax-1))/(M_trash(jmax)+M_DM(jmax)+M_S(jmax))

IF(0.170 - bf_c_control .le. 0.) THEN
   SPL = rho_c
ELSE 
   SPR = rho_c
END IF

IF(abs(0.170 - bf_c_control) .le. eps_bar_frac) THEN
   COND_1 = .true.
   rho_trash(1) = rho_c
END IF

END DO

FORMAT = "(7(e15.7))"
OPEN(30, file="data_DM_plus_STARS_T_NOT_CONSTANT.dat")
WRITE(30,*) "r"," rr"," M_DM"," M_S"," M_GSDM_T"," rho_dm", " rho_gas_DM_plus_STARS_T"
DO i = 1, jmax
   WRITE(30,FORMAT) r(i)/kpc, rr(i)/kpc, M_DM(i)/msol, M_S(i)/msol, M_trash(i)/msol, rho_dm(i), rho_trash(i)
END DO

CLOSE(30)
WRITE(*,*) "The right density is: ", rho_trash(1)
WRITE(*,*) "Baryonic fraction:     ", bf_c_control
WRITE(*,*) ""
WRITE(*,*) ""
WRITE(*,*) " #################################################################### "
WRITE(*,*) " #                   IRON DIFFUSION IN THE ICM                      # "
WRITE(*,*) " #            START OF THE SECOND PART OF THE PROJECT               # "
WRITE(*,*) " #################################################################### "
WRITE(*,*) ""
WRITE(*,*) ""

! Variable used in the second part of the code
bg_zfe = 0.4*zfesol
v_turb = 260e5
l_turb = 35.*kpc
param_D(1) = 1.32e29

WRITE(*,*) "Insert diffusion parameter D in units of 10^29 cm^2/s ? "
READ(*,*) param_D(1)
param_D(1) = param_D(1)*1.e29 

WRITE(*,*) "Diffusion parameter D: ", param_D(1)
OPEN(66, FILE="rho_fe_2D.dat")
OPEN(67, FILE="M_Fe.dat")
OPEN(68, FILE="Z_Fe.dat")
OPEN(69, FILE="position_sn_explosion.dat")
OPEN(70, FILE="temperature_time.dat")
OPEN(78, FILE="gradz_t.dat")
OPEN(79, FILE="M_Fe_theo.dat")
output_time_control = 0.

DO i=1, jmax
   rho_stellar(i) = (10.**12/2./pi*msol)*(a/(rr(i)*(rr(i)+a)**3))
ENDDO

OPEN(30, FILE="data_DM_plus_STARS_T_NOT_CONSTANT.dat")
!FIle format: |, r  rr  M_DM  M_S  M_GSDM_T  rho_dm  rho_gas_DM_plus_STARS_T
FORMAT = "(7(e15.7))"
REWIND(30)
READ(30,*)


DO i=1, jmax
! Gas density previously calculated in the first part. CENTERED ON i + 1/2 
   READ(30,FORMAT) r(i), rr(i), M_DM(i), M_S(i), M_GSDM_T(i), rho_dm(i), rho_gas(i)
   r(i) = r(i) * kpc
   rr(i) = rr(i) * kpc
   M_DM(i) = M_DM(i) * msol
   M_S(i) = M_S(i) * msol
   M_GSDM_T(i) = M_GSDM_T(i) * msol
END DO

CLOSE(30)

! Defining the temperature and density of the paper
OPEN(20, FILE="temperature.dat")
OPEN(30, FILE="paper_profile.dat")
WRITE(30,*) "t_paper rho_paper t_profile rho_gas"
DO i=1, jmax
n_e_paper(i) = (4.6d-2/(1+(rr(i)/kpc/57)**2)**1.8)+(4.8d-3/(1+(rr(i)/kpc/200)**2)**0.87)
n_tot_paper(i) = 1.83*n_e_paper(i)
rho_paper(i) = 1.937d-24*n_e_paper(i)
READ(20,"(1(e15.7))") t_profile(i)
t_paper(i) = 7.*(1+(rr(i)/kpc/71.)**3)/(2.3+(rr(i)/kpc/71.)**3)*1.16d7

WRITE(30,"(4(e15.7))") t_paper(i), rho_paper(i), t_profile(i), rho_gas(i)
END DO
CLOSE(30)
CLOSE(20)

rho_g_mean(1) = rho_gas(1)
DO i=2, jmax-1
   rho_g_mean(i) = (rho_gas(i)+rho_gas(i-1))/2.
ENDDO
rho_g_mean(jmax) = (rho_gas(jmax)+rho_gas(jmax-1))/2.

zfeout = 0.42*zfesol

DO i = 1, jmax
   Z_Fe(i) = 1.4*0.3*((2.2+(rr(i)/kpc/80.)**3)/(1.+(rr(i)/kpc/80.)**3))*zfesol
   Z_Fe(i) = Z_Fe(i) - zfeout
   rho_fe(i) = Z_Fe(i)*rho_gas(i)/1.4
END DO

DO i=2, jmax-1
   grad_zfe(i) = (Z_Fe(i)-Z_Fe(i-1))/(rr(i)-rr(i-1))
ENDDO
grad_zfe(1) = 0.
grad_zfe(jmax) = 0.        ! this required by boundary condition so that it remains constant
open(77,file='grad_zfe.dat')
! output gradzfe to file
DO i=1, jmax
   WRITE(77,"(e15.7)") grad_zfe(i)
END DO
CLOSE(77)

open(29,file='temporary_data.dat')
do i = 1, jmax                                              ! RE-ENABLE FOR DEBUG
   write(29,"(e15.7)") Z_Fe(i)/zfesol
enddo
close(29)

! We first need to compute the gradient of the netallicity on r(i), which is just the difference between Z_fe(i)-Z_fe(i-1) divided
! by the interval (rr(i)-rr(i-1). 
! Also, we need to compute the density in r(i) which is just the mean value from rho_gas(i) and rho_gas(i-1)

C_param = 0.5 ! Courant parameter
end_of_time = .FALSE.
TIME_PASSED = 0.




dt=(r(2)-r(1))**2/(2*param_D(1))*C_param

random_approach_sn = 0
WRITE(*,*) "For how many years you want to evolve the system: [Gyr] (timestep =", dt/years, "years, or: ", dt/gyr, " Gyr"
READ(*,*) time_wanted
time_wanted = 8.
time_wanted = time_wanted * 3.1536d16
timesteps = INT(time_wanted / dt)
! ALLOCATE(rho_fe_2D(timesteps+1, jmax))
! ALLOCATE(M_Fe(timesteps+1, jmax))
WRITE(*,*) "How often in time you want to output the data to file?:  [Gyr] "
READ(*,*) delta_t_debug
delta_t_debug = 1.
delta_t_debug = delta_t_debug * 3.1536d16

k_time_step = 1

WRITE(*,*) "Do you want to start with NO IRON in the ICM? [1 = TRUE, 0 = FALSE]"
READ(*,*) start_with_zero_iron
start_with_zero_iron = 1

IF(start_with_zero_iron .eq. 1) THEN
WRITE(*,*) "How many supernovae per century you want to inject?  [N_SN/century]" 
READ(*,*) SN_unit
END IF
! WRITE(*,*) "Do you want to inject supernovae with a random approach? [1 = TRUE, 0 = FALSE]"
! READ(*,*) random_approach_sn

if (time_wanted .LE. 13.8*gyr) then
   initial_time = 13.8*gyr - time_wanted
ELSE
   initial_time = 8.0*gyr  ! Needed if we want to evolve the system to a time greater than 13.8 Gyr, considering that literature
                           ! suggest that Perseus cluster is 5-6 Gyr old
END IF

IF(start_with_zero_iron .eq. 1) THEN
   Z_Fe = 0.
   rho_fe = Z_Fe*rho_gas/1.4
END IF


! Just to know how much iron until 100kpc
! IRON MASS IN THE ICM


M_Fe = 0.
M_Fe(1) = rho_fe(1)*(vol_r(1))
DO i=2, jmax
   M_Fe(i) = (rho_fe(i) * vol_r(i)) ! Here volume is really a delta_V 
ENDDO
M_Fe_initial = sum(M_Fe)

DO i=1, index_100kpc
   M_Fe_100kpc_initial = M_Fe_100kpc_initial + M_Fe(i)
END DO


! START OF THE WHILE LOOP FOR THE IRON DIFFUSION
! The loop will stop when the time passed is equal to the time wanted by the user

DO WHILE(end_of_time .eqv. .FALSE.)
CALL CPU_TIME(time_start_debug)

! SOURCE TERM 
IF(start_with_zero_iron .eq. 1) THEN
   DO i = 1, jmax
   ! source_term_rho(i) = SN_unit * rho_stellar(i)*(4.7d-21)*((5.7*gyr + time_passed)/(13.8*gyr))**(-1.1)
   IF (time_wanted .LE. 13.8*gyr) THEN
   source_term_rho(i) = (3.13d-21 * rho_stellar(i) * SN_unit)*(((initial_time+time_passed)/(13.8*gyr))**(-1.26)) & 
                     + (4.7d-20*rho_stellar(i)*Z_Fe(i)/1.4 * ((initial_time+time_passed)/(13.8*gyr))**(-1.26))
   ELSE
   source_term_rho(i) = (3.13d-21 * rho_stellar(i) * SN_unit) & 
   + (4.7d-20*rho_stellar(i)*Z_Fe(i))
   END IF
   rho_fe(i) = rho_fe(i) + source_term_rho(i)*dt
   !Z_Fe(i) = rho_fe(i)/rho_gas(i)*1.4 ! Uncomment to disable the diffusion
   ENDDO
END IF

! ! DIFFUSION 
DO i=2, jmax-1                        
   rho_fe(i) = rho_fe(i) + (dt/1.4)*(r(i+1)**2*param_D(1)*rho_g_mean(i+1)*grad_zfe(i+1)- &
               r(i)**2*param_D(1)*rho_g_mean(i)*grad_zfe(i)) &
               /((r(i+1)**3-r(i)**3)/3.) ! COMMENT TO DISABLE DIFFUSION
   Z_Fe(i) = rho_fe(i)/rho_gas(i)*1.4 ! just updating the value at radius r(i) for the iron density, the diffusion is going on.
END DO                                           
Z_Fe(1) = Z_Fe(2)
Z_Fe(jmax) = Z_Fe(jmax-1)        
rho_fe(1)=rho_fe(2)              ! Boundary conditions
rho_fe(jmax) = rho_fe(jmax-1)
k_time_step = k_time_step+1

DO i=2, jmax-1
   grad_zfe(i) = (Z_Fe(i)-Z_Fe(i-1))/(rr(i)-rr(i-1)) ! COMMENT WHEN DIFFUSION DISABLE
ENDDO
grad_zfe(1) = 0.
grad_zfe(jmax) = 0.
! CONDITIONS FOR IRON CONSERVATION

! IRON MASS IN THE ICM
M_Fe = 0.
M_fe(1) = rho_fe(1)*(vol_r(1))
DO i=2, jmax
   M_Fe(i) = (rho_fe(i) * vol_r(i))
ENDDO
M_fe_tot = sum(M_Fe)
M_Fe_100kpc = 0.
M_Fe_100kpc_final = 0.
! Check the amount of iron until some kpc (Here we use 100 but user can change it at runtime)
DO i=1, index_100kpc
   M_Fe_100kpc = M_Fe_100kpc + M_Fe(i)
END DO
M_Fe_100kpc_final = M_Fe_100kpc
FORMAT = "(A32,F8.2,A23)"
! Check and print if the iron mass under some kpc has been halved
IF((M_Fe_100kpc .LE. M_Fe_100kpc_initial/2.) .AND. stop_debug_mass_index .EQ. 0 .AND. SN_unit .EQ. 0) THEN
   WRITE(*,*)
   WRITE(*,*) "--------------------------------------------------------------------- "
   WRITE(*,FORMAT) "The iron mass in the ICM under ", kpc_check_mass ," kpc has been halved"
   WRITE(*,*) "M_i: ", M_Fe_100kpc_initial/msol, " M_Fe(t): ", M_Fe_100kpc/msol, "tau: ", TIME_PASSED/3.15d16, "Gyr"
   WRITE(*,*) "--------------------------------------------------------------------- "
   WRITE(*,*)
   tau_diffusion = TIME_PASSED
   stop_debug_mass_index = 1
END IF

IF ((output_time_control .GE. delta_t_debug) .OR. (time_passed .EQ. dt) &
      .OR. (time_passed  .GE. time_wanted - dt)) THEN
   
   time_scale_debug = (time_end_debug - time_start_debug) * (timesteps-k_time_step)
   WRITE(*,*) "Fe <100 Kpc: ", M_Fe_100kpc/msol, &
         " Time passed: ", TIME_PASSED/3.15d16, "Gyr", "dt", dt/3.15d16, "Gyr"
   DO i=1, jmax
      WRITE(66,"(E15.7)", advance='no') rho_fe(i)
      WRITE(67,"(E15.7)", advance='no') M_Fe(i)/msol
      WRITE(68, "(E15.7)", advance='no') Z_Fe(i)/zfesol
      WRITE(70, "(E15.7)", advance='no') t_r(i)
      WRITE(78, "(E15.7)", advance='no') grad_zfe(i)
   END DO

   WRITE(79, "(E15.7)") (TIME_PASSED*(5d-22)*1d12)
   WRITE(67,"(3(E15.7))", advance='no') M_fe_tot/msol, M_Fe_100kpc/msol, TIME_PASSED
   WRITE(66,*)
   WRITE(67,*)
   WRITE(68,*)
   WRITE(70,*)
   WRITE(78,*)
   output_time_control = 0.
ELSE
   output_time_control = output_time_control + dt
END IF



CALL CPU_TIME(time_end_debug)

IF(TIME_PASSED .GE. time_wanted - dt) THEN
   end_of_time = .TRUE.
ELSE
   end_of_time = .FALSE.
END IF

TIME_PASSED = TIME_PASSED + dt

ENDDO ! SUPER ENDDO FINISH

! DEBUG INFORMATION TO FILE
! D dt rho_g_central M_Fe_initial M_Fe_final M_Fe_100kpc_i M_Fe_100kpc_f t_f SN_unit t_i kpc_diffusion tau_diffusion 

close(67)
close(66)
close(69)
close(68)
close(70)
close(78)
close(79)

OPEN(666, FILE="recap_data.dat", status="old", position="append")
WRITE(666,"(12(E11.4))") param_D(1), dt, rho_gas(1), M_Fe_initial/msol, &
M_fe_tot/msol, M_Fe_100kpc_initial/msol, &
M_Fe_100kpc_final/msol, TIME_PASSED/3.15d16, SN_unit, initial_time/3.15d16, &
tau_diffusion/3.15d16, kpc_check_mass
CLOSE(666)


CONTAINS

subroutine random_choice(arr, n, out_arr, probabilities)
   REAL*8, intent(in) :: arr(:)
   integer, intent(in) :: n
   REAL*8, intent(out) :: out_arr(n)
   REAL*8, intent(in) :: probabilities(:)
   integer :: i, j
   REAL*8 :: cum_sum, rand_num


   do i = 1, n
      cum_sum = 0.0
      call random_seed()
      call random_number(rand_num)
      
      do j = 1, size(arr)
           cum_sum = cum_sum + probabilities(j)
            if (rand_num <= cum_sum) then
               out_arr(i) = arr(j)
               exit
            end if
      end do
   end do
end subroutine random_choice

! subroutine NORMALIZE_GRID_BY_DENSITY(x,y)
!    REAL*8, intent(in) :: x
!    REAL*8, intent(out) :: y
!    integer :: i
!    REAL*8 :: xs, kpc
!    kpc = 3.084e21
!    xs = 435.7*kpc
!    y = 1/((x/xs)*(1+x/xs)**2)

! end subroutine

END PROGRAM
