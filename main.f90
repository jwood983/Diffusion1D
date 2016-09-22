program diffusion_solvers
   use solvers
   real(wp), dimension(0:151) :: u, w, x
   real(wp), dimension(2) :: timers
   real(wp) :: t, tmax, t0, t1, t2
   integer :: i, j, k
   
   call initialize(x,u); w = u
   t = 0d0; tmax = 1d-1; timers=0d0
   call WriteInitial(x,u)
   j=0
   do
      j = j+1
! evolve via different methods
      call cpu_time(t0)
      call crank(t,x,u) ! crank-nicolson step
      call cpu_time(t1)
      call super(t,x,w) ! super time stepping step
      call cpu_time(t2)

! update timers
      timers(1) = timers(1) + t1 - t0
      timers(2) = timers(2) + t2 - t1
      t = t+dta
      if(t > tmax) exit
   enddo
   call WriteFinal(x,t,u,w)
   print *,"---"
   print *," CN:",timers(1)/real(j)/1d-6," microseconds/step"
   print *,"STS:",timers(2)/real(j)/1d-6," microseconds/step"
   print '(i0,a)',j," total steps employed"
   
 contains
   !> @brief output the distribution in curve format
   subroutine WriteInitial(x,u)
     real(wp), dimension(0:151), intent(in) :: x, u
     integer :: i, ierr
     
     open(unit=10,file='diffusion0000.curve',iostat=ierr)
     if(ierr/=0) then
        print *,"error: unable to open file: ierr=",ierr
        stop
     endif
     write(10,*) "# temp_cn"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), u(i)
     enddo
     write(10,*) "# temp_sts"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), u(i)
     enddo
     write(10,*) "# exact"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), u(i)
     enddo
     write(10,*) "# sts-cn"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), 0d0
     enddo
     close(10)
   end subroutine WriteInitial
   
   !> @brief write the distribution in curve format
   subroutine WriteFinal(x,t,u,w)
     real(wp), dimension(0:151), intent(in) :: x,u,w
     real(wp), intent(in) :: t
     real(wp) :: exact
     integer :: i, ierr
     
     open(unit=10,file='diffusion0001.curve',iostat=ierr)
     if(ierr/=0) then
        print *,"error: unable to open file: ierr=",ierr
        stop
     endif
     write(10,*) "# temp_cn"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), u(i)
     enddo
     write(10,*) "# temp_sts"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), w(i)
     enddo
     write(10,*) "# exact"
     do i=1,150
        exact = sin(pi*x(i))*exp(-pi*pi*t)
        write(10,'(f11.4,2x,es15.4)') x(i), exact
     enddo
     write(10,*) "# sts-cn"
     do i=1,150
        write(10,'(f11.4,2x,es15.4)') x(i), w(i)-u(i)
     enddo
     close(10)
   end subroutine WriteFinal
end program diffusion_solvers
