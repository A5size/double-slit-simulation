module qm_lib
  use iso_c_binding
  Implicit none

    Double precision, save :: PI = 2.0d0*acos(0.0d0)
    Complex(kind(0d0)), save :: ei = (0.0d0, 1.0d0)
    Double precision, save :: time = 0.0d0
    Double precision, save :: dt, dx=0.2d0, dy=0.2d0

    Integer, save :: nx=400, ny=400 

    Complex(kind(0d0)), save, allocatable :: PSI2d(:,:)
    Complex(kind(0d0)), save, allocatable :: v2d(:,:)

    Double precision, save, allocatable :: field_x(:)
    Double precision, save, allocatable :: field_y(:)

  contains

    subroutine set_time(time_in) bind(c)
      Real(c_double), intent(in) :: time_in
      time = time_in
    end subroutine set_time

    subroutine get_time(time_out) bind(c)
      Real(c_double), intent(out) :: time_out
      time_out = time
    end subroutine get_time
    
    subroutine set_dx(dx_in) bind(c)
      Real(c_double), intent(in) :: dx_in
      dx = dx_in
    end subroutine set_dx

    subroutine get_dx(dx_out) bind(c)
      Real(c_double), intent(out) :: dx_out
      dx_out = dx
    end subroutine get_dx

    subroutine set_dy(dy_in) bind(c)
      Real(c_double), intent(in) :: dy_in
      dy = dy_in
    end subroutine set_dy

    subroutine get_dy(dy_out) bind(c)
      Real(c_double), intent(out) :: dy_out
      dy_out = dy
    end subroutine get_dy
    
    subroutine set_nx(nx_in) bind(c)
      Integer(c_int), intent(in) :: nx_in
      nx = nx_in
    end subroutine set_nx

    subroutine set_ny(ny_in) bind(c)
      Integer(c_int), intent(in) :: ny_in
      ny = ny_in
    end subroutine set_ny

    subroutine get_nx(nx_out) bind(c)
      Integer(c_int), intent(out) :: nx_out
      nx_out = nx
    end subroutine get_nx

    subroutine get_ny(ny_out) bind(c)
      Integer(c_int), intent(out) :: ny_out
      ny_out = ny
    end subroutine get_ny
    
    subroutine set_psi2d(psi2d_real_in, psi2d_imag_in) bind(c)
      real(c_double), intent(in) :: psi2d_real_in(ny, nx), psi2d_imag_in(ny, nx)
      PSI2d(:, :) = cmplx(transpose(psi2d_real_in(:, :)), transpose(psi2d_imag_in(:, :)), kind=kind(0.0d0))
    end subroutine set_psi2d

    subroutine get_psi2d(psi2d_real_out, psi2d_imag_out) bind(c)
      real(c_double), intent(out) :: psi2d_real_out(ny, nx), psi2d_imag_out(ny, nx)
      psi2d_real_out(:, :) = real(transpose(PSI2d(:, :)))
      psi2d_imag_out(:, :) = aimag(transpose(PSI2d(:, :)))
    end subroutine get_psi2d

    subroutine set_v2d(v2d_real_in, v2d_imag_in) bind(c)
      real(c_double), intent(in) :: v2d_real_in(ny, nx), v2d_imag_in(ny, nx)
      v2d(:, :) = cmplx(transpose(v2d_real_in(:, :)), transpose(v2d_imag_in(:, :)), kind=kind(0.0d0))
    end subroutine set_v2d

    subroutine get_v2d(v2d_real_out, v2d_imag_out) bind(c)
      real(c_double), intent(out) :: v2d_real_out(ny, nx), v2d_imag_out(ny, nx)
      v2d_real_out(:, :) = real(transpose(v2d(:, :)))
      v2d_imag_out(:, :) = aimag(transpose(v2d(:, :)))
    end subroutine get_v2d
    
    subroutine set_field_x(field_x_in) bind(c)
      real(c_double), intent(in) :: field_x_in(nx)
      field_x(:) = field_x_in(:)
    end subroutine set_field_x
    
    subroutine get_field_x(field_x_out) bind(c)
      real(c_double), intent(out) :: field_x_out(nx)
      field_x_out(:) = field_x(:)
    end subroutine get_field_x

    subroutine set_field_y(field_y_in) bind(c)
      real(c_double), intent(in) :: field_y_in(ny)
      field_y(:) = field_y_in(:)
    end subroutine set_field_y
    
    subroutine get_field_y(field_y_out) bind(c)
      real(c_double), intent(out) :: field_y_out(ny)
      field_y_out(:) = field_y(:)
    end subroutine get_field_y
    
    subroutine init_2d() bind(c)
      Integer :: i
      
      dt = (dx*dy)/2.0d0
      
      allocate(PSI2d(0:(nx-1),0:(ny-1)))
      allocate(v2d(0:(nx-1),0:(ny-1)))
      allocate(field_x(0:(nx-1)))
      allocate(field_y(0:(ny-1)))
      
      PSI2d(:,:) = (0.0d0, 0.0d0)
      v2d(:,:) = (0.0d0, 0.0d0)
      
      Do i=0, (nx-1)
         field_x(i) = dx*dble(i - nx/2)
      end Do
      
      Do i=0, (ny-1)
         field_y(i) = dy*dble(i - ny/2)
      end Do
      
    end subroutine init_2d

    subroutine put_gaussian(x, y, kx, ky, cc) bind(c)
      Real(c_double) :: x, y
      Real(c_double) :: kx, ky
      Real(c_double) :: cc

      Integer :: i, j
      Real(c_double) :: D
      Real(c_double) :: ux, uy
      Real(c_double) :: cx, cy
      
      D = sqrt((2.0d0*cc)/PI)
      
      Do i=0, (nx-1)
         Do j=0, (ny-1)
            cx = field_x(i) - x
            cy = field_y(j) - y
            ux = kx*field_x(i)
            uy = ky*field_y(j)
            PSI2d(i,j) = D*exp(-cc*(cx**2 + cy**2) + ei*(ux + uy))
         end Do
      end Do

    end subroutine put_gaussian

    subroutine set_v2d_box(sx, sy, ex, ey, Re_c, Im_c) bind(c)
      Real(c_double) :: sx, sy, ex, ey, Re_c, Im_c 
      Real(c_double) :: tem
      Integer :: sxi, syi, exi, eyi
      
      If(sx > ex) then
         tem = sx
         sx = ex
         ex = tem
      end If
      
      If(sy > ey) then
         tem = sy
         sy = ey
         ey = tem
      end If
      
      sxi = nx/2 + anint(sx/dx)
      exi = nx/2 + anint(ex/dx)
      
      syi = ny/2 + anint(sy/dy)
      eyi = ny/2 + anint(ey/dy)
      
      If(sxi<0) then
         sxi = 0
      end If
      
      If((nx-1)<exi) then
         exi = nx-1
      end If
      
      If(syi<0) then
         syi = 0
      end If
      
      If((ny-1)<eyi) then
         eyi = ny-1
      end If
      
      v2d(sxi:exi, syi:eyi) = Re_c + ei*Im_c
      
    end subroutine set_v2d_box
    
    subroutine step_2d(n) bind(c)
      Integer(c_int) :: n
      Integer :: i
      
      Do i=1, n
         call step_2d_sub()
      end Do

    end subroutine step_2d

    subroutine step_2d_sub()
      Integer :: i, j 
      
      Complex(kind(0d0)) :: k1(0:(nx-1),0:(ny-1))
      Complex(kind(0d0)) :: k2(0:(nx-1),0:(ny-1))
      Complex(kind(0d0)) :: k3(0:(nx-1),0:(ny-1))
      Complex(kind(0d0)) :: k4(0:(nx-1),0:(ny-1))
      
      Complex(kind(0d0)) :: tem1
      Complex(kind(0d0)) :: tem2
      Complex(kind(0d0)) :: tem3
      
      !Double precision :: cal_e = 2.0d-50
      
      k1(:,:) = (0.0d0, 0.0d0)
      k2(:,:) = (0.0d0, 0.0d0)
      k3(:,:) = (0.0d0, 0.0d0)
      k4(:,:) = (0.0d0, 0.0d0)
      
      Do j=1, (ny-2)
         Do i=1, (nx-2)
            
            !If(abs(PSI2d(i,j)) < cal_e) then
            !   cycle
            !end If

            !----------------------x--------------------------
            tem1 = PSI2d(i+1,j)
            tem2 = PSI2d(i,j)
            tem3 = PSI2d(i-1,j)
            k1(i,j) = (tem1 - 2.0d0*tem2 + tem3)/(dx**2)
            
            !----------------------y--------------------------
            tem1 = PSI2d(i,j+1)
            tem2 = PSI2d(i,j)
            tem3 = PSI2d(i,j-1)
            k1(i,j) = k1(i,j) + (tem1 - 2.0d0*tem2 + tem3)/(dy**2)

            !----------------------k--------------------------
            k1(i,j) = (ei/2.0d0)*k1(i,j)

            !----------------------v--------------------------
            k1(i,j) = k1(i,j) - ei*v2d(i,j)*PSI2d(i,j)

         end Do
      end Do

      Do j=1, (ny-2)
         Do i=1, (nx-2)
            
            !If(abs(PSI2d(i,j)) < cal_e) then
            !   cycle
            !end If
            
            !----------------------x--------------------------
            tem1 = PSI2d(i+1,j) + k1(i+1,j)*(dt/2.0d0)
            tem2 = PSI2d(i,j) + k1(i,j)*(dt/2.0d0)
            tem3 = PSI2d(i-1,j) + k1(i-1,j)*(dt/2.0d0)
            k2(i,j) = (tem1 - 2.0d0*tem2 + tem3)/(dx**2)

            !----------------------y--------------------------
            tem1 = PSI2d(i,j+1) + k1(i,j+1)*(dt/2.0d0)
            tem2 = PSI2d(i,j) + k1(i,j)*(dt/2.0d0)
            tem3 = PSI2d(i,j-1) + k1(i,j-1)*(dt/2.0d0)
            k2(i,j) = k2(i,j) + (tem1 - 2.0d0*tem2 + tem3)/(dy**2)

            !----------------------k--------------------------
            k2(i,j) = (ei/2.0d0)*k2(i,j)

            !----------------------v--------------------------
            k2(i,j) = k2(i,j) - ei*v2d(i,j)*(PSI2d(i,j) + k1(i,j)*(dt/2.0d0))

         end Do
      end Do

      Do j=1, (ny-2)
         Do i=1, (nx-2)
            
            !If(abs(PSI2d(i,j)) < cal_e) then
            !   cycle
            !end If

            !----------------------x--------------------------
            tem1 = PSI2d(i+1,j) + k2(i+1,j)*(dt/2.0d0)
            tem2 = PSI2d(i,j) + k2(i,j)*(dt/2.0d0)
            tem3 = PSI2d(i-1,j) + k2(i-1,j)*(dt/2.0d0)
            k3(i,j) = (tem1 - 2.0d0*tem2 + tem3)/(dx**2)

            !----------------------y--------------------------
            tem1 = PSI2d(i,j+1) + k2(i,j+1)*(dt/2.0d0)
            tem2 = PSI2d(i,j) + k2(i,j)*(dt/2.0d0)
            tem3 = PSI2d(i,j-1) + k2(i,j-1)*(dt/2.0d0)
            k3(i,j) = k3(i,j) + (tem1 - 2.0d0*tem2 + tem3)/(dy**2)

            !----------------------k--------------------------
            k3(i,j) = (ei/2.0d0)*k3(i,j)

            !----------------------v--------------------------
            k3(i,j) = k3(i,j) - ei*v2d(i,j)*(PSI2d(i,j) + k2(i,j)*(dt/2.0d0))

         end Do
      end Do

      Do j=1, (ny-2)
         Do i=1, (nx-2)

            !If(abs(PSI2d(i,j)) < cal_e) then
            !   cycle
            !end If

            !----------------------x--------------------------
            tem1 = PSI2d(i+1,j) + k3(i+1,j)*dt
            tem2 = PSI2d(i,j) + k3(i,j)*dt
            tem3 = PSI2d(i-1,j) + k3(i-1,j)*dt
            k4(i,j) = (tem1 - 2.0d0*tem2 + tem3)/(dx**2)

            !----------------------y--------------------------
            tem1 = PSI2d(i,j+1) + k3(i,j+1)*dt
            tem2 = PSI2d(i,j) + k3(i,j)*dt
            tem3 = PSI2d(i,j-1) + k3(i,j-1)*dt
            k4(i,j) = k4(i,j) + (tem1 - 2.0d0*tem2 + tem3)/(dy**2)

            !----------------------k--------------------------
            k4(i,j) = (ei/2.0d0)*k4(i,j)

            !----------------------v--------------------------
            k4(i,j) = k4(i,j) - ei*v2d(i,j)*(PSI2d(i,j) + k3(i,j)*dt)

         end Do
      end Do

      Do j=1, (ny-2)
         Do i=1, (nx-2)

            !If(abs(PSI2d(i,j)) < cal_e) then
            !   cycle
            !end If

            PSI2d(i,j) = PSI2d(i,j) + ((k1(i,j) + 2.0d0*k2(i,j) + 2.0d0*k3(i,j) + k4(i,j))/6.0d0)*dt
            
         end Do
      end Do
      time = time + dt

    end subroutine step_2d_sub
    
end module qm_lib
