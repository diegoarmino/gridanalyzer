program kde
!----------------------------------------------------------------------------
! GAUSSIAN KERNEL DENSITY ESTIMATOR ANALYSIS
! This program performs KDE over a selection of a trajectory in order to obtain
! a 3D map of probability and the corresponding free energy estimation.
!----------------------------------------------------------------------------
implicit none

character(len=50) :: xyzfile, csigma, cnlines, cdelta
integer           :: nsteps, nlines, i, j, k, n, ierr, cnt
integer           :: ndx_i, ndx_j, ndx_k
double precision  :: xx, yy, zz
double precision  :: pi, pi2, sqrt_pi2, gnorm, sigma, delta

integer, dimension(3) :: grid_size, indx, gcentr_indx
integer, dimension(3) :: eval_strt, eval_end, eval_ndx
double precision,dimension(3)  :: gorig, gmax, xyz, gcentr_xyz, eval_xyz, gauss
double precision, allocatable, dimension(:,:,:) :: grid

call get_command_argument(1, xyzfile)
call get_command_argument(2, cnlines)
call get_command_argument(3, cdelta)
call get_command_argument(4, csigma)

read(cnlines,*) nlines
read(cdelta,*) delta
read(csigma,*) sigma

pi = 3.1415926535897931
pi2 = 2.0d0*pi
sqrt_pi2 = sqrt(pi2)
gnorm = 1/sqrt_pi2**3
gnorm = gnorm/sigma


open(unit=001,file=xyzfile, iostat=ierr)

read(001,*) gorig(1), gorig(2), gorig(3)
read(001,*) gmax(1),  gmax(2),  gmax(3)


!gorig = gorig - delta
!gmax  = gmax  + delta 

grid_size = int( (gmax - gorig)/delta )

allocate ( grid( grid_size(1), grid_size(2), grid_size(3) ) )

grid = 0d0

cnt=1
do n=1,nlines
   read(001,*) gcentr_xyz(1), gcentr_xyz(2), gcentr_xyz(3)
   !Determine the closest grid point to xx, yy, zz
   gcentr_indx = 1 + nint( (gcentr_xyz - gorig) / delta )
   !Determine grid points in which to evaluate the gaussian.
   eval_strt = gcentr_indx - 5
   eval_end  = gcentr_indx + 5
   !Check that evaluation bubble is inside the grid.
   if ( eval_strt(1) < 1 ) eval_strt(1) = 1
   if ( eval_strt(2) < 1 ) eval_strt(2) = 1
   if ( eval_strt(3) < 1 ) eval_strt(3) = 1
   if ( eval_end(1) > grid_size(1) ) eval_end(1) = grid_size(1)
   if ( eval_end(2) > grid_size(2) ) eval_end(2) = grid_size(2)
   if ( eval_end(3) > grid_size(3) ) eval_end(3) = grid_size(3)

   !Evaluate gaussian in grid bubble.
   do ndx_i=eval_strt(1),eval_end(1)
   do ndx_j=eval_strt(2),eval_end(2)
   do ndx_k=eval_strt(3),eval_end(3)
      ! Convert grid point indexes to cartesian coordinates to evalate gaussian.
      eval_ndx = (/ ndx_i,ndx_j,ndx_k /)
      eval_xyz = gorig + delta * dble(eval_ndx-1)
      ! Evaluate gaussian.
      gauss = exp( -0.5d0*((eval_xyz - gcentr_xyz)/sigma)**2 )
      ! Accumulate on grid.
      grid(ndx_i,ndx_j,ndx_k) = grid(ndx_i,ndx_j,ndx_k) + gnorm*gauss(1)*gauss(2)*gauss(3)
!      if (mod(cnt,1000)==0) write(*,*) grid(ndx_i,ndx_j,ndx_k)  ! DEBUG
   end do
   end do
   end do
   if (mod(cnt,1000)==0) write(*,'(A,I7,A,I7,A,F7.2,A)') "Finished step ",cnt, " / ", nlines," (",dble(cnt)/dble(nlines)*100d0,"%)"
   cnt=cnt+1
end do
close(unit=001)

! Normalize grid
grid = grid/cnt

! Print probability distribution.
open(unit=002,file="probability.dx",iostat=ierr)
write(002,'(A)') '#' 
write(002,'(A)') '#' 
write(002,'(A,3I4)') 'object 1 class gridpositions counts ',grid_size(1),grid_size(2),grid_size(3) 
write(002,'(A,3F8.4)') 'origin ',gorig(1),gorig(2),gorig(3)
write(002,'(A,3F4.1)') 'delta',delta,0d0,0d0
write(002,'(A,3F4.1)') 'delta',0d0,delta,0d0
write(002,'(A,3F4.1)') 'delta',0d0,0d0,delta
write(002,'(A,3I4)') 'object 2 class gridconnections counts ',grid_size(1),grid_size(2),grid_size(3) 
write(002,'(A,I8)') 'object 3 class array double rank 0 items ',grid_size(1)*grid_size(2)*grid_size(3) 

cnt=1
do i=1,grid_size(1)
do j=1,grid_size(2)
do k=1,grid_size(3)
   if (mod(cnt,3)==0) then
      if (grid(i,j,k) > 1.0d-20) then
         write(002,'(D11.4)') grid(i,j,k)
      else
         write(002,'(D11.4)') 0d0
      end if 
   else
      if (grid(i,j,k) > 1.0d-12) then
         write(002,'(D11.4)',advance='no') grid(i,j,k)
      else 
         write(002,'(D11.4)',advance='no') 0d0
      end if 
   end if
   cnt = cnt + 1
end do
end do
end do
write(002,*)
close(unit=002)

! Compute free energy grid and print to file.
open(unit=003,file="free_energy.dx",iostat=ierr)
write(003,'(A)') '#' 
write(003,'(A)') '#' 
write(003,'(A,3I4)') 'object 1 class gridpositions counts ',grid_size(1),grid_size(2),grid_size(3) 
write(003,'(A,3F8.4)') 'origin ',gorig(1),gorig(2),gorig(3)
write(003,'(A,3F4.1)') 'delta',delta,0d0,0d0
write(003,'(A,3F4.1)') 'delta',0d0,delta,0d0
write(003,'(A,3F4.1)') 'delta',0d0,0d0,delta
write(003,'(A,3I4)') 'object 2 class gridconnections counts ',grid_size(1),grid_size(2),grid_size(3) 
write(003,'(A,I8)') 'object 3 class array double rank 0 items ',grid_size(1)*grid_size(2)*grid_size(3) 

cnt=1
do i=1,grid_size(1)
do j=1,grid_size(2)
do k=1,grid_size(3)
   if (grid(i,j,k) > 1d-30) then
      grid(i,j,k) = -log(grid(i,j,k))
   else
      grid(i,j,k) = 60.0d0
   end if

   if (mod(cnt,3)==0) then
      write(003,'(F11.5)') grid(i,j,k)
   else
      write(003,'(F11.5)',advance='no') grid(i,j,k)
   end if
   cnt = cnt + 1
end do
end do
end do
write(003,*)
close(unit=003)

end program kde
