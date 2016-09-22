module solvers
   implicit none
   integer, parameter :: wp = selected_real_kind(15,307)
   real(wp), parameter :: pi = 3.141592653d0
   real(wp), parameter :: nu = 0.005d0
   real(wp), parameter :: CFL = 0.3d0
   integer, parameter :: nx = 400
   real(wp), allocatable :: ts(:)
   real(wp) :: dx, dt, C, dxi, dtp, dta
   
 contains
   !> @brief initialize the distribution
   subroutine initialize(x,u)
     real(wp), dimension(0:nx+1), intent(out) :: x, u
     real(wp) :: nn
     integer :: i, N
     
     dx = 1d0/real(nx - 1)
     dxi = 1d0/dx
     x(1) = 0d0; x(0) = x(1) - dx
     do i=2,nx+1
        x(i) = x(i-1) + dx
     enddo
     
     do i=0,nx+1
        u(i) = sin(pi*x(i))
     enddo
     
     dt = CFL*dx*dx
     print *,"Parameters:"
     print *,"   dx:",dx
     print *,"   dt:",dt
     print *,"  cfl:",dt/dx/dx
     
! compute STS substeps
     dtp = dt
     dta = CFL*dx
     
     nn = FindRoot(1d0, dtp, dta)
     N = int(nn) + 1
     allocate(ts(N))
     dtp = CorrectTimeStep(N, dta)
     CALL ComputeSubSteps(dtp, ts)
   end subroutine initialize
   
   !> @brief Compute Crank-Nicolson step
   subroutine crank(t,x,u)
     real(wp), intent(in) :: t
     real(wp), dimension(0:nx+1), intent(in) :: x
     real(wp), dimension(0:nx+1), intent(inout) :: u
     real(wp), dimension(1:nx) :: a, b, c, d, un
     real(wp) :: r
     integer :: i
     
     r = dta*dxi*dxi;
     do i=1,nx
        a(i) = -r/2d0
        b(i) = 1d0 + r
        c(i) = -r/2d0
        d(i) = u(i)*(1d0 - r) + (0.5d0*r)*(u(i+1) + u(i-1))
     enddo
     
     call tridiag(a,b,c,d,un)
     u(1:nx) = un(1:nx)
     call bc(t+dta,x,u)
     
   end subroutine crank
   
   !> @brief Tridiagonal matrix solver algorithm
   subroutine tridiag(a,b,c,d,un)
      real(wp), dimension(1:nx), intent(in) :: a,b,c,d
      real(wp), dimension(1:nx), intent(out) :: un
      real(wp), dimension(1:nx)  :: dp, cp
      real(wp) :: m
      integer :: i
! initialize c' and d'
      cp(1) = c(1)/b(1)
      dp(1) = d(1)/b(1)
! upwards sweep
      do i=2,nx
         m = b(i) - cp(i-1)*a(i)
         cp(i) = c(i)/m
         dp(i) = (d(i) - dp(i-1)*a(i))/m
      enddo
! downwards sweep
      un(nx) = dp(nx)
      do i=nx-1,1,-1
         un(i) = dp(i) - cp(i)*un(i+1)
      enddo
   end subroutine tridiag

   !> @brief Compute heat equation using super time-stepping algorithm
   subroutine super(t,x,u)
     real(wp), intent(in) :: t
     real(wp), dimension(0:nx+1), intent(in) :: x
     real(wp), dimension(0:nx+1), intent(inout) :: u
     real(wp), dimension(0:nx+1) :: un, up, phi
     real(wp) :: tau, r, tstep
     integer :: i, m
     
     tstep = t
     phi = 0d0; un = u
     do m=1,N
        tau = ts(m)
        r = tau*dxi
        tstep = tstep + tau
! compute grad(u) using backwards difference
        do i=1,nx+1
           phi(i) = (up(i) - up(i-1))*dxi
        enddo
! compute u(n+1) using forwards difference
!   - this is effectively a Lax-Wendroff method for parabolic methods
        do i=1,nx
           un(i) = un(i) + r*(phi(i+1) - phi(i))
        enddo
        call bc(tstep,x,un)
     enddo
     u(1:150) = un(1:150)
   end subroutine super
   
   !> @brief Determines optimal number of STS substeps
   function FindRoot(x0, dtr, dta) result(NS)
     real(wp), intent(in) :: x0, dtr, dta
     real(wp) :: NS, a, b, c, N1, rt_nu
     
     NS = x0 + 1d0
     N1 = x0
     rt_nu = sqrt(nu)
     do while(abs(NS - N1) >= 1d-5)
        NS = N1
        a = (1d0 - rt_nu)/(1d0 + rt_nu)
        b = a**(2d0*NS)
        c = (1d0 - b)/(1d0 + b)
        N1 = NS + (dta - dtr*NS/(2d0*rt_nu)*c)/(dtr/(2d0*rt_nu)*(c - 2d0*NS*b*log(a)*(1d0+c)/(1d0+b)))
     enddo
   end function FindRoot
   
   !> @brief Corrects the parabolic time step for use in STS method
   pure function CorrectTimeStep(n0, dta) result(dtr)
     integer, intent(in) :: n0
     real(wp), intent(in) :: dta
     real(wp) :: a, b, c, dtr, rt_nu
     
     rt_nu = sqrt(nu)
     a = (1d0 - rt_nu)/(1d0 + rt_nu)
     b = a**(2*n0)
     c = (1d0 - b)/(1d0 + b)
     dtr = dta*2d0*rt_nu/(n0*c)
   end function CorrectTimeStep
   
   !> @brief Determines the substeps for STS method
   subroutine ComputeSubSteps(dtex, tau)
     real(wp), intent(in) :: dtex
     real(wp), intent(out) :: tau(:)
     real(wp) :: s
     integer :: j
     
     s = real(size(tau))
     do j=1,size(tau)
        tau(j) = dtex / ((-1d0 + nu)*cos(((2d0*j - 1d0)*pi)/(2d0*s))  +  1d0  +  nu);
     enddo
   end subroutine ComputeSubSteps
   
   !> @brief Adjusts the boundary conditions
   subroutine bc(t,x,u)
     real(wp), intent(in) :: t
     real(wp), dimension(0:nx+1), intent(in) :: x
     real(wp), dimension(0:nx+1), intent(inout) :: u
     
     u(0) = sin(pi*x(0))*exp(-pi*pi*t)
     u(nx+1) = sin(pi*x(nx+1))*exp(-pi*pi*t)
   end subroutine

end module solvers
