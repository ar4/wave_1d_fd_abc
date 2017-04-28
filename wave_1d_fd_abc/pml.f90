module pml

  implicit none

contains

  subroutine step(f1, f2, lphi1, rphi1, lphi2, rphi2, sigma, sigma_x,  &
      model_padded, dt, dx, sources, sources_x, num_steps, nx_padded,  &
      num_sources, source_len, abc_width)

    integer, intent (in) :: nx_padded
    integer, intent (in) :: num_sources
    integer, intent (in) :: source_len
    integer, intent (in) :: abc_width
    real, intent (in out), dimension (nx_padded) :: f1
    real, intent (in out), dimension (nx_padded) :: f2
    real, intent (in out), dimension (abc_width) :: lphi1
    real, intent (in out), dimension (abc_width) :: rphi1
    real, intent (in out), dimension (abc_width) :: lphi2
    real, intent (in out), dimension (abc_width) :: rphi2
    real, intent (in), dimension (abc_width) :: sigma
    real, intent (in), dimension (abc_width) :: sigma_x
    real, intent (in), dimension (nx_padded) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    real, intent (in), dimension (num_sources, source_len) :: sources
    integer, intent (in), dimension (num_sources) :: sources_x
    integer, intent (in) :: num_steps

    integer :: step_idx
    logical :: even

    do step_idx = 1, num_steps
    even = (mod (step_idx, 2) == 0)
    if (even) then
      call one_step(f2, f1, lphi2, rphi2, lphi1, rphi1, sigma, sigma_x,&
        model_padded, dt, dx, sources, sources_x, step_idx, nx_padded, &
        num_sources, source_len, abc_width)
    else
      call one_step(f1, f2, lphi1, rphi1, lphi2, rphi2, sigma, sigma_x,&
        model_padded, dt, dx, sources, sources_x, step_idx, nx_padded, &
        num_sources, source_len, abc_width)
    end if
    end do

  end subroutine step


  subroutine one_step(f, fp, lphi, rphi, lphip, rphip, sigma, sigma_x, &
      model_padded, dt, dx, sources, sources_x, step_idx, nx_padded,   &
      num_sources, source_len, abc_width)

    integer, intent (in) :: nx_padded
    integer, intent (in) :: num_sources
    integer, intent (in) :: source_len
    integer, intent (in) :: abc_width
    real, intent (in out), dimension (nx_padded) :: f
    real, intent (in out), dimension (nx_padded) :: fp
    real, intent (in), dimension (abc_width) :: lphi
    real, intent (in), dimension (abc_width) :: rphi
    real, intent (in out), dimension (abc_width) :: lphip
    real, intent (in out), dimension (abc_width) :: rphip
    real, intent (in), dimension (abc_width) :: sigma
    real, intent (in), dimension (abc_width) :: sigma_x
    real, intent (in), dimension (nx_padded) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    real, intent (in), dimension (num_sources, source_len)  :: sources
    integer, intent (in), dimension (num_sources) :: sources_x
    integer, intent (in) :: step_idx

    integer :: i

    ! left PML
    do i = 5, abc_width
    call fd_pml(f, fp, lphi, lphip, sigma, sigma_x, model_padded,      &
      dt, dx, nx_padded, abc_width, i, i)
    end do

    ! interior (no PML)
    do i = abc_width + 1, nx_padded - abc_width
    call fd_interior(f, fp, model_padded, dt, dx, nx_padded, i)
    end do

    ! right PML
    do i = nx_padded - abc_width + 1, nx_padded - 4
    call fd_pml(f, fp, rphi, rphip, sigma, sigma_x, model_padded,      &
      dt, dx, nx_padded, abc_width, i, nx_padded - i + 1)
    end do

    ! source term
    do i = 1, num_sources
    call add_source(fp, model_padded, dt, sources(i, step_idx),        &
      sources_x(i), nx_padded, abc_width)
    end do

  end subroutine one_step


  subroutine fd_interior(f, fp, model_padded, dt, dx, nx_padded, i)

    integer, intent (in) :: nx_padded
    real, intent (in), dimension (nx_padded) :: f
    real, intent (in out), dimension (nx_padded) :: fp
    real, intent (in), dimension (nx_padded) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    integer, intent (in) :: i

    real :: f_xx

    f_xx = second_x_deriv(f, i, dx, nx_padded)
    fp(i) = (model_padded(i)**2 * dt**2 * f_xx + 2 * f(i) - fp(i))

  end subroutine fd_interior

  
  subroutine fd_pml(f, fp, phi, phip, sigma, sigma_x, model_padded,    &
      dt, dx, nx_padded, abc_width, i, j)

    integer, intent (in) :: nx_padded
    integer, intent (in) :: abc_width
    real, intent (in), dimension (nx_padded) :: f
    real, intent (in out), dimension (nx_padded) :: fp
    real, intent (in), dimension (abc_width) :: phi
    real, intent (in out), dimension (abc_width) :: phip
    real, intent (in), dimension (abc_width) :: sigma
    real, intent (in), dimension (abc_width) :: sigma_x
    real, intent (in), dimension (nx_padded) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    integer, intent (in) :: i
    integer, intent (in) :: j

    ! Based on PML explanation by Steven G. Johnson
    ! http://math.mit.edu/~stevenj/18.369/pml.pdf
    !
    ! Scalar wave equation: u_xx = u_tt/c^2  (1)
    !
    ! Split (1) into:
    ! phi_t = u_x  (2a)
    ! u_t = c^2 * phi_x  (2b)
    !
    ! Verification:
    ! x partial derivative of (2a): phi_(x,t) = u_xx  (3a)
    ! t partial derivative of (2b): u_tt = c^2 * phi_(x,t)  (3b)
    ! substitute (3a) in (3b): u_tt = c^2 * u_xx, which is (1).
    !
    ! For the PML region, we use the wavefield at
    ! x' = x * (1 + i *  sigma / omega), as this will exponentially
    ! damp the wavefield. Changing coordinates back to x,
    ! u_x -> u_x / (1 + i * sigma / omega). Applying this to (2):
    ! phi_t = u_x / (1 + i * sigma / omega)  (4a)
    ! u_t = c^2 * phi_x / (1 + i * sigma / omega)  (4b).
    !
    ! Assuming that u and phi can be written in the form
    ! A * exp(i(kx - omega*t)), u_t = -i*omega*u. Applying
    ! this to (4) and multiplying both sides by the denominator
    ! of the right hand sides:
    ! (1 + i * sigma / omega) * (-i*omega) * phi = u_x  (5a)
    ! (1 + i * sigma / omega) * (-i*omega) * u = c^2 * phi_x  (5a)
    !
    ! Multiplying,
    ! (1 + i * sigma / omega) * (-i*omega) = -i*omega + sigma.
    ! Using this in (5), and noting that -i*omega*u = u_t,
    ! phi_t = u_x - sigma * phi  (6a)
    ! u_t = c^2 * phi_x - sigma * u  (6b)
    !
    ! Taking the x partial derivative of (6a), the t partial derivative
    ! of (6b), and substituting, as in (3), we get
    ! u_tt = c^2 * (u_xx - (sigma * phi)_x) - sigma * u_t (7)
    ! Expanding the derivative (sigma * phi)_x using the product rule:
    ! u_tt = c^2 * (u_xx - (sigma_x * phi + sigma * phi_x))
    !        - sigma * u_t (8)
    !
    ! Using central differences for u_tt and u_t:
    ! (u(t+1) - 2*u(t) + u(t-1))/dt^2 = c^2 * (u_xx - (sigma_x * phi
    !                                                  + sigma * phi_x))
    !                                   - sigma * (u(t+1) - u(t-1))/(2*dt)
    ! Rearranging:
    ! u(t+1) = c^2 * dt^2 / (1 + dt * sigma / 2)
    !          * (u_xx - (sigma_x * phi + sigma * phi_x))
    !          + dt * sigma / (2 + dt * sigma) * u(t-1)
    !          + 1 / (1 + dt * sigma / 2) * (2 * u(t) - u(t-1))  (9)
    !
    ! Using a forward difference for u_t in (6a):
    ! (phi(t+1) - phi(t))/dt = u_x - sigma * phi(t),
    ! Rearranging:
    ! phi(t+1) = dt * u_x + phi(t) * (1 - dt * sigma)  (10)
    ! 
    ! 
    ! (9) shows how to update u, and (10) shows how to update phi.

    real :: f_xx
    real :: f_x
    real :: phi_x

    f_xx = second_x_deriv(f, i, dx, nx_padded)
    f_x = first_x_deriv(f, i, dx, nx_padded)
    phi_x = first_x_deriv(phi, j, dx, nx_padded)

    ! (9)
    fp(i) = model_padded(i)**2 * dt**2 / (1 + dt * sigma(j) / 2)       &
            * (f_xx - (sigma_x(j) * phi(j) + sigma(j) * phi_x))        &
            + dt * sigma(j) / (2 + dt * sigma(j)) * fp(i)              &
            + 1 / (1 + dt * sigma(j) / 2) * (2 * f(i) - fp(i))

    ! (10)
    phip(j) = dt * f_x + phi(j) * (1 - dt * sigma(j))

  end subroutine fd_pml


  subroutine add_source(fp, model_padded, dt, source, source_x,        &
      nx_padded, abc_width)

    integer, intent (in) :: nx_padded
    real, intent (in out), dimension (nx_padded) :: fp
    real, intent (in), dimension (nx_padded) :: model_padded
    real, intent (in) :: dt
    real, intent (in)  :: source
    integer, intent (in) :: source_x
    integer, intent (in) :: abc_width

    integer :: sx

    sx = source_x + abc_width + 1;
    fp(sx) = fp(sx) + (model_padded(sx)**2 * dt**2 * source)

  end subroutine add_source


  pure function first_x_deriv(f, i, dx, nx_padded)

    integer, intent (in) :: nx_padded
    real, intent (in), dimension (nx_padded) :: f
    integer, intent (in) :: i
    real, intent (in) :: dx

    real :: first_x_deriv

    !f_x = (                                                           &
    !  5*f(i-6)-72*f(i-5)                                              &
    !  +495*f(i-4)-2200*f(i-3)                                         &
    !  +7425*f(i-2)-23760*f(i-1)                                       &
    !  +23760*f(i+1)-7425*f(i+2)                                       &
    !  +2200*f(i+3)-495*f(i+4)                                         &
    !  +72*f(i+5)-5*f(i+6))/(27720*dx)
    first_x_deriv = (                                                  &
      1/280*f(i-4) - 4/105*f(i-3) + 1/5*f(i-2) - 4/5*f(i-1)            &
      -1/280*f(i+4) + 4/105*f(i+3) - 1/5*f(i+2) + 4/5*f(i+1)) / dx

  end function first_x_deriv


  pure function second_x_deriv(f, i, dx, nx_padded)

    integer, intent (in) :: nx_padded
    real, intent (in), dimension (nx_padded) :: f
    integer, intent (in) :: i
    real, intent (in) :: dx

    real :: second_x_deriv

    second_x_deriv = (                                                 &
      -735*f(i-8)+15360*f(i-7)                                         &
      -156800*f(i-6)+1053696*f(i-5)                                    & 
      -5350800*f(i-4)+22830080*f(i-3)                                  & 
      -94174080*f(i-2)+538137600*f(i-1)                                & 
      -924708642*f(i+0)                                                & 
      +538137600*f(i+1)-94174080*f(i+2)                                & 
      +22830080*f(i+3)-5350800*f(i+4)                                  & 
      +1053696*f(i+5)-156800*f(i+6)                                    & 
      +15360*f(i+7)-735*f(i+8))/(302702400*dx**2)


  end function second_x_deriv

end module pml
