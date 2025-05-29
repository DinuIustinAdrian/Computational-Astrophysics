MODULE DATA
   real*8 :: pi,pc,kpc,yr,kbol,mu,mp, c
   !INTEGER :: jmax
   
   
   parameter(pi=3.141592)
   parameter(pc=3.085d18)
   parameter(kpc=1000.*pc)
   parameter(yr=3.156d7)
   parameter(kbol=1.38d-16)
   parameter(mu=0.61)
   parameter(mp=1.67d-24)
   parameter(c=2.998d10)
   !parameter(jmax=1000)
   END MODULE DATA
   
   
   PROGRAM MAIN
   USE DATA
   IMPLICIT NONE
   
   LOGICAL :: integration_condition, DEBUG_POSITION, verbose_debug, cooling_OFF
   INTEGER :: i, ref_sys, ncycle, jmax, time_control_int, SPECIAL_PRINT, number_of_injected_cell
   CHARACTER(100) :: FORMAT
   REAL*8 :: x0, xmax, gamma, c2, cfl, t_f, t, cv, eval_time, E_tot, E_kin, cool
   REAL*8, ALLOCATABLE, DIMENSION(:) :: xb(:), xa(:), d(:), p(:), e(:), v(:), &
   s(:), temp(:), dxa(:), dxb(:), g2a(:), g31a(:), g2b(:), g31b(:), &
   dvl1a(:), dvl1b(:), q1(:), DIV_v(:), ms(:), dstar(:), F1(:), M(:), e_dstar(:), &
   F3(:), vstar(:), F2(:)
   REAL*8 :: dt_min, rho_0, temp_0, deltax, debug_time, debug_time_setup, L_x, cpu_time_0, cpu_time_1
   
   
   !E_0 = 10.e51 ! [erg/s] energy inject by SNe
   
   WRITE(*,*) "Initial density? [IN GR/CM^3] DEFAULT = 2*10^-24 : INSERT 0 FOR DEFAULT"
   read(*,*) rho_0
   rho_0 = rho_0 * 2.d-24
   if(rho_0 .eq. 0) then
      rho_0 = 2d-24
   end if
   

   write(*,*) "Initial temperature? [IN K] DEFAULT = 10^4 : INSERT 0 FOR DEFAULT"
   read(*,*) temp_0
   temp_0 = temp_0 * 10**4
   if(temp_0 .eq. 0) then
      temp_0 = 10**4
   end if
   

   cooling_OFF = .FALSE.
   WRITE(*,*) "Do you want cooling to be DISABLED? [DEFAULT = .FALSE.]"
   READ(*,*) cooling_OFF

   x0 = 0.
   xmax = 70* pc
   WRITE(*,*) "Maximum region radius? [IN PARSEC]"
   READ(*,*) xmax
   xmax = xmax * pc
   gamma = 1.4
    
                                       ! PROGRAM BEGINS !
   WRITE(*,*) "GRID POINTS?: "
   READ(*,*) jmax
   
   ! WRITE(*,*) " Long verbose debug to file at each timestep: Default = TRUE"
   ! READ(*,*) verbose_debug
   verbose_debug = .true.
   
   
   ALLOCATE(xa(jmax))     ! Grid A defined in 0,1,2,3
   ALLOCATE(xb(jmax))     ! Grid B defined in 1/2, 3/2, 5/2 etc
   ALLOCATE(d(jmax))      ! Density                defined in xb
   ALLOCATE(e(jmax))      ! Energy                 defined in xb
   ALLOCATE(v(jmax))      ! Velocity               defined in xb
   ALLOCATE(s(jmax))      ! ??                     defined in xb   
   ALLOCATE(temp(jmax))   ! Temperature            defined in xb
   
   ALLOCATE(p(jmax))
   ALLOCATE(q1(jmax))
   
   ALLOCATE(dxa(jmax))
   ALLOCATE(dxb(jmax))
   
   ALLOCATE(g2a(jmax))
   ALLOCATE(g31a(jmax))
   ALLOCATE(g2b(jmax))
   ALLOCATE(g31b(jmax))
   ALLOCATE(dvl1a(jmax))
   ALLOCATE(dvl1b(jmax), DIV_v(jmax), ms(jmax), dstar(jmax), F1(jmax))
   ALLOCATE(M(jmax), e_dstar(jmax), F3(jmax), F2(jmax), vstar(jmax))
   
   ! #####################################
   ! ########### GRID GENERATION #########
   ! #####################################

   
   xa(1) = x0+(xmax-x0)*(-1.)/(jmax-1.)
   !write(*,*) xa(1)
   DO i = 2, jmax
      xa(i) = xa(1)+(xmax-x0)*(i-1.)/(jmax-1.)
   END DO
   deltax = xa(3)-xa(2)
   DO i=1, jmax-1
      xb(i)=0.5*(xa(i)+xa(i+1))
   END DO
   xb(jmax)=xb(jmax-1)+(xb(jmax-1)-xb(jmax-2))
   
   ! COMPUTING GRADIENTS !
   DO i=2, jmax-1
      dxa(i) = xa(i+1) - xa(i)
      dxb(i) = xb(i) - xb(i-1)
   ENDDO
   dxa(1) = xa(2)-xa(1)
   dxa(jmax) = dxa(jmax-1)
   dxb(1) = dxb(2)
   dxb(jmax) = xb(jmax)-xb(jmax-1)
   
   FORMAT = "(4(e15.7))"
   OPEN(20, FILE='grid.dat')
   DO i=2, jmax
      WRITE(20,FORMAT) xa(i)/pc, xb(i)/pc, dxa(i), dxb(i)
      WRITE(*,FORMAT) xa(i), xb(i), dxa(i), dxb(i)
   ENDDO
   CLOSE(20)
   
   ref_sys = 2 ! Default spherical coordinates
                  ! UNCOMMENT FROM LINE 128 to 146 to enable the choice of the coordinate system

   ! DO WHILE((ref_sys .ne. 1) .and. (ref_sys .ne. 2) .and. (ref_sys .ne. 0))
   ! WRITE(*,*) "--------------------------------------------------------"
   ! WRITE(*,*) " Which coordinate system to use:"
   ! WRITE(*,*) "   [1]. Cartesian"
   ! WRITE(*,*) "   [2]. Spherical"
   ! WRITE(*,*) "   [0]. QUIT"
   ! WRITE(*,*) "--------------------------------------------------------"
   ! READ(*,*)  ref_sys
   ! IF((ref_sys .ne. 1) .and. (ref_sys .ne. 2) .and. (ref_sys .ne. 0)) THEN
   ! WRITE(*,*) "!_Value not valid_!"
   ! ELSE 
   
   ! END IF
   
   ! END DO
   !WRITE(*,*) ref_sys
   
   ! Now we compute the metric scale factor based
   ! on the coordinates system chosen.
   

   IF(ref_sys .eq. 1) THEN
      DO i = 1, jmax
         g2a(i) = 1.
         g31a(i) = 1.
         g2b(i) = 1.
         g31b(i) = 1.
      ENDDO
   
      DO i = 1, jmax-1
      dvl1a(i) = xa(i+1)-xa(i)
      ENDDO
      dvl1a(jmax) = dvl1a(jmax-1)
   
      DO i = 2, jmax
      dvl1b(i) = xb(i)-xb(i-1)
      ENDDO
      dvl1b(1) = dvl1b(2)
   END IF
   
   IF(ref_sys .eq. 2) THEN
      DO i = 1, jmax
         g2a(i) = xa(i)
         g31a(i) = xa(i)
         g2b(i) = xb(i)
         g31b(i) = xb(i)
      ENDDO
   
      DO i = 1, jmax-1
      dvl1a(i) = (xa(i+1)**3-xa(i)**3)/3.
      ENDDO
      dvl1a(jmax) = dvl1a(jmax-1)
   
      DO i = 2, jmax
      dvl1b(i) = (xb(i)**3-xb(i-1)**3)/3.
      ENDDO
      dvl1b(1) = dvl1b(2)
   END IF
   
   
   gamma = 1.4
   cv = 1.99d8
   
   t = 0.
   t_f = 5.1
   WRITE(*,*) "Insert simulation duration: ............ [ input * 10^5 years]"
   READ(*,*) t_f
   t_f = t_f * 10**5 * yr
   c2=3.
   
   WRITE(*,*) "Insert how often to print debug data (in 10^5 years): "
   READ(*,*) debug_time_setup
   debug_time_setup = debug_time_setup * 10**5 * yr
   debug_time = 0.
   
   
   ! DO i = 1, jmax                                  ! This is for the shock tube
   !    IF (xa(i) .LE. 0.5) THEN
   !       d(i) = 1.0
   !       p(i) = 1.0
   !    ELSE 
   !       d(i) = 0.125
   !       p(i) = 0.1
   !    END IF
   !    v(i) = 0.
   !    e(i)=p(i)/(gamma-1.) ! Energy
   ! ENDDO
   
   ncycle = 0
   

   time_control_int = 1.
   
   
   DO i = 2, jmax
      d(i) = rho_0
      v(i) = 0.
   END DO
   
   
   
   DO i=2, jmax
      temp(i) = temp_0
   END DO
   
   DO i=2, jmax
      e(i) = temp(i)*cv*d(i)
   END DO
   
   DO i=2, jmax
      p(i) = (gamma-1.)*e(i)
   END DO
   
   
   ! #####################################
   ! ########### ENERGY INJECTION #########
   ! #####################################

   number_of_injected_cell = 2 ! Default value of 2 cells, can be changed
   e(2) = 1.d51 / (4./3.*pi*xa(4)**3)
   p(2) = (gamma-1.)*e(2)
   temp(2) = e(2)/cv/d(2)

   DO i = 3, number_of_injected_cell+1
      e(i) = e(i-1)
      p(i) = p(i-1)
      temp(i) = (i-1)
   END DO
   
   
   ! DO i=0, jmax
   !    write(*,*) d(i), v(i), temp(i), e(i), p(i)
   ! END DO
   
   
     
 
   ! Now that we defined the initial condition the 
   ! code can start integrating over time

   ! #####################################
   ! ########### MAIN INTEGRATION #########
   ! #####################################
   
   dt_min = 1.d30
   cfl=0.01
   
   
      OPEN(30, FILE="velocity.dat")
      OPEN(31, FILE="pressure.dat")
      OPEN(32, FILE="energy.dat")
      OPEN(33, FILE="density.dat")
      OPEN(34, FILE="temperature.dat")
      OPEN(35, FILE="time.dat")
   
   
   open(36, file="luminosity_x.dat")
   open(37, file="shock_position.dat")
   open(38, file="total_energy.dat")
   
   SPECIAL_PRINT = 0
   DO WHILE(t < t_f)
      CALL cpu_time(cpu_time_0) ! for debug purposes
   ncycle = ncycle+1

   p=(gamma-1.)*e
   
   dt_min=1.d30
   
   
   !WRITE(*,*) cfl
   
   cfl=min(0.5,1.1*cfl) ! This will increase dt until possible
   DO i=2, jmax-1
      dt_min = min(dt_min,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gamma*p(i)/d(i))))
   END DO
   dt_min = cfl*dt_min
   t=t+dt_min
   
   ! First step is to update velocity due to the gradient of P
   !WRITE(*,*) v(0), v(1), v(jmax), v(jmax-1)
   DO i = 2, jmax-1
      v(i) = v(i) - dt_min*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i))

   END DO
   
   CALL boundary_conditions(v, jmax)
                                 
   DO i=2, jmax-1 ! Defining the q1 var. coeff for the artificial viscosity
      IF ((v(i+1)-v(i))<0.) THEN
         q1(i) = C2 * d(i)*(v(i+1)-v(i))**2 ! C2 goes from 1 to 3
      ELSE
         q1(i) = 0.
      END IF 
   END DO
   CALL BCb(v, jmax)
   
   ! Then cycle to subtract from v(i) the artificial viscosity term
   DO i=2, jmax-1
      v(i)=v(i)-dt_min*2.*(q1(i)-q1(i-1))/((d(i)+d(i-1))*dxb(i))
   END DO
   !v(1) = 0.
   CALL boundary_conditions(v, jmax)
   
   ! Energy update step
   DO i=2, jmax-1
      e(i)=e(i)-dt_min*q1(i)*(v(i+1)-v(i))/dxa(i)
   END DO
   CALL BCb(e, jmax)
   
   ! Then cycle to update the quantity considering compression heating
   
   DO i=2, jmax-1
      DIV_v(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
   END DO
   CALL boundary_conditions(DIV_v, jmax)
   
   DO i=2, jmax-1
   e(i)=e(i)*(1.-0.5*dt_min*(gamma-1.)*DIV_v(i))/(1.+0.5*dt_min*(gamma-1.)*DIV_v(i))
   END DO
   !DIV_v(1) = 0.
   CALL BCb(e, jmax)
   
   
   
      ! #####################################
      ! ########### COOLING FUNCTION #########   
      ! Radiation cooling function
   DO i=2, jmax-1
      e(i) = e(i) - dt_min * Cool(temp(i), d(i), cooling_OFF)
   END DO
   CALL BCb(e, jmax)
   
   DO i=2, jmax-1
      temp(i)=e(i)/cv/d(i)
      IF(temp(i)<10**4) THEN
         temp(i)=10**4
      END IF
      e(i)=cv*d(i)*temp(i)
   END DO
   CALL boundary_conditions(e, jmax)
   CALL boundary_conditions(temp, jmax)
   
   ! MOMENTUM DENSITY
   DO i=2, jmax-1   
         ms(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
   END DO
   !ms(1) = 0.
   CALL boundary_conditions(ms, jmax)
   
   ! Now we use the donor cell  method, which says that the quantity is constant
   ! over a cell volume. We define dstar(i) as the coefficient for the upwind method
   
   DO i=2, jmax-1
      IF(v(i)>0.) THEN
         dstar(i)=d(i-1) ! [Eq. 47] page 12 paper
      ELSE
         dstar(i)=d(i)
      END IF
   END DO
   dstar(jmax) = dstar(jmax-1)
   dstar(1) = dstar(3)
   
   DO i=2, jmax-1
      F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)     ! Centered at i [ Eq. 51] page 13
   END DO                                    ! We are using some simplified
                                             ! version. In relity(?), we should
                                             ! use an average of g2a and g31a
                                             ! but in 1D it's fine to proceed this way
   DO i=2, jmax-1
         M(i)=dstar(i)*v(i)
   END DO
   !M(1) = 0.
   CALL boundary_conditions(M, jmax)
   
   
   DO i=2, jmax-1
      IF (v(i)>0.) THEN
         e_dstar(i)=e(i-1)/d(i-1)   !! at i !!
      ELSE
         e_dstar(i)=e(i)/d(i)
      END IF
   END DO
      e_dstar(jmax)=e_dstar(jmax-1)
      e_dstar(1)=e_dstar(3)
   
   
                                          !Now we can update the density
   DO i=2, jmax-1
      d(i)=d(i)-dt_min*(F1(i+1)-F1(i))/dvl1a(i)
   END DO 
   CALL BCb(d, jmax)
      
   
   DO i=2, jmax-1
      F2(i)=e_dstar(i)*M(i)*g2a(i)*g31a(i)            
   END DO
   !F2(1) = 0.
   CALL boundary_conditions(F2, jmax)
   
   DO i=2, jmax-1
      e(i)=e(i)-dt_min*(F2(i+1)-F2(i))/dvl1a(i)
   END DO
   CALL BCb(e, jmax)
   
   
   DO i=2, jmax-1
      IF ((v(i-1)+v(i))*0.5>0) THEN
         vstar(i)=v(i-1)       
      ELSE
         vstar(i)=v(i)
      END IF
   END DO
   CALL BCb(vstar, jmax) 
   
   DO i=2, jmax-1
      F3(i)=vstar(i+1)*0.5*(M(i)+M(i+1))*g2b(i)*g31b(i)  
   END DO
   
   DO i=2, jmax-1
      ms(i)=ms(i)-dt_min/dvl1b(i)*(F3(i)-F3(i-1))
   END DO
   s(1) = 0.
   CALL boundary_conditions(s, jmax)
   
   DO i=2, jmax-1
      v(i)=2.*ms(i)/(d(i)+d(i-1))
   END DO
   v(1) = 0.
   CALL BCb(v, jmax)
   

   DO i=2, jmax-1
   temp(i) = e(i)/cv/d(i)
   END DO
   CALL boundary_conditions(temp, jmax)
   
   ! FORMAT = "(1A,I5,4A,E14.5,A2,E14.5,3A,E14.5,1A,E14.5,4A,E14.5,1A,E14.5,1A,E14.5)"
   
   IF (debug_time > debug_time_setup) THEN
      DEBUG_POSITION = .TRUE.
      debug_time = 0.
   ELSE IF(debug_time < debug_time_setup) THEN
      DEBUG_POSITION = .FALSE.
      debug_time = debug_time + dt_min
   END IF
   
   IF((t .le. 1.e4*yr + dt_min/2) .and. (t .ge. 1.e4*yr - dt_min/2)) THEN
      SPECIAL_PRINT = 1
   ELSE 
      SPECIAL_PRINT = 0
   END IF
   
       ! Calculate the luminosity in the x-ray band, basically we need to integrate over the volume the emissivity
   ! and then multiply by 4*pi*r^2, but we do it only when the time is right, based on the last if statement above
   ! controlling that the temperature is in the right range, above 10^6 K
   L_x=0.
   DO i=2, jmax
      IF(temp(i) .ge. 1.d6) THEN ! ONLY IF THE TEMPERATURE IS ABOVE 10^6 K
         L_x = L_x + 4./3.*pi*(xb(i)**3-xb(i-1)**3)*Cool(temp(i), d(i), .false.)
      ELSE
         L_x=L_x
      END IF
   END DO

   WRITE(36,*) L_x, t/yr

   IF(((debug_time .GE. debug_time_setup) .AND. (debug_time .LE. debug_time_setup+dt_min)) &
   .OR. (SPECIAL_PRINT .EQ. 1)) THEN
      IF(verbose_debug) THEN ! Writing to file only if verbose is set to TRUE
         write(30,*) (v(i), i=2,jmax)
         write(31,*) (p(i), i=2,jmax)
         write(32,*) (e(i), i=2,jmax)
         write(33,*) (d(i), i=2,jmax)
         write(34,*) (temp(i), i=2,jmax)
      END IF

   write(35,*) t/yr
   
   
   ! Find the radius of the shock at a certain time, considering the peak in density over the whole domain, at that time
   ! and then calculate the shock velocity

    DO i=2, jmax
      IF(d(i) == maxval(d)) THEN
         write(*,*) "SHOCK POSITION", xa(i)/pc
         write(*,*) "SHOCK VELOCITY", v(i)/100000.0  !in km/s
         write(37,*) t/yr, xb(i)/pc, v(i)/100000.0
      END IF
   END DO


! calculate the total energy as a function of time and write it to file
   E_tot = 0.
   E_kin = 0.
   DO i=2, jmax
      E_tot = E_tot + 4./3.*pi*(xa(i)**3-xa(i-1)**3)*e(i)
      E_kin = E_kin + 2./3.*pi*d(i)*(xa(i)**3-xa(i-1)**3)*v(i)*v(i)
   END DO
   write(38,*) t/yr, E_tot, E_kin
   END IF

   
   CALL cpu_time(cpu_time_1) ! for debug purposes
   !write(*,*) "CPU TIME", cpu_time_1-cpu_time_0
   !write(*,*) "Remaing time", int((t_f-t)/dt_min*(cpu_time_1-cpu_time_0)), " seconds"
   
   END DO ! END DO WHILE
   open(20,file='results.dat')
   
    do i=2,jmax  !! write the results in the file "results.dat"
       write (20,1000) xa(i),xb(i),d(i),v(i),e(i)/d(i),p(i)
    end do
   1000 format(6(1pe12.4))
   
   close(20)
   CLOSE(30)
   CLOSE(31)
   CLOSE(32)
   CLOSE(33)
   CLOSE(34)
   CLOSE(35)
   CLOSE(36)
   CLOSE(37)
   close(38)
   
   
   END PROGRAM
   
   
   
   
   SUBROUTINE boundary_conditions(z1, dim) !corrette BC per velocità e momento (riflessione)
   USE DATA
   IMPLICIT NONE
   INTEGER :: dim
   real*8, dimension (dim) :: z1
   
   
   z1(2)=0.
   z1(1)=z1(3)
   !z1(dim)=z1(dim-1)
   z1(1)=-z1(2)       !! ouflow !!
   z1(dim)=-z1(dim-2)
   z1(dim-1) = 0
   
   END SUBROUTINE boundary_conditions
   
   SUBROUTINE BCb(z2, dim) ! BC di outflow tradizionali
   USE DATA
   IMPLICIT NONE
   INTEGER :: dim
   real*8, dimension (dim) :: z2
   z2(1)=z2(2)
   z2(dim)=z2(dim-1)
   END SUBROUTINE BCb
   
   Real*8 FUNCTION Cool(Temp1, d1, disable_cooling)
   USE DATA
   IMPLICIT NONE
   Real*8:: Temp1, d1, kev
   LOGICAL:: disable_cooling
   kev = 1.16d7
   IF(disable_cooling) THEN
      cool = 0.
   ELSE

      IF(Temp1 > 2.32e5) THEN
         cool = 1.d-22*(8.6*1.e-3*(Temp1/kev)**(-1.7)+0.058*(Temp1/kev)**0.5+0.063)*(d1/2.17d-24)**2

      ELSE IF(Temp1 < 2.32e5 .AND. Temp1 > 0.0017235*kev) THEN

         COOL = (6.72d-22 * (Temp1/kev/0.02)**0.6)*(d1/2.17d-24)**2
      ELSE IF(Temp1 < 0.0017235*kev) THEN

         COOL = (1.544d-22*(Temp1/kev/0.0017235)**6)*(d1/2.17d-24)**2
      ELSE
      cool = 0.
      END IF
   END IF
   END FUNCTION Cool
   
   
   
   
   
   
    
   
   
   
   