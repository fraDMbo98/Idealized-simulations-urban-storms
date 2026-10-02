program ppert
!
! ----------------------------------------------------------
!
implicit none
!
integer, parameter                      :: nx = 250
integer, parameter                      :: ny = 250
integer, parameter                      :: nz = 63
double precision, parameter             :: dx = 1.0e3
double precision, parameter             :: dy = 1.0e3

!
integer                                 :: itime, itnum
integer                                 :: irec
integer                                 :: i, j, k
double precision, parameter             :: error = 0.0001
double precision                        :: diff
double precision, dimension(nx, ny, nz) :: dn0, pi0, theta0
double precision, dimension(nx, ny, nz) :: qv0
double precision, dimension(nx, ny, nz) :: pi, qv
double precision, dimension(nx, ny, nz) :: rtotal, theta
double precision, dimension(nx, ny, nz) :: u, v, w, z
double precision, dimension(nx, ny, nz) :: lappd, pd
double precision, dimension(nx, ny, nz) :: lappb, buoy, pb
double precision, dimension(nx, ny, nz) :: dpd, dpb
double precision                        :: dz, dzl, dzr
double precision                        :: thv
double precision                        :: thv0, thv0m1
double precision                        :: thv0p1
double precision                        :: Dijk, oldpd
double precision			:: oldpb
character*13                            :: dirname
character*10                            :: dirname2
character*26                            :: binfile
character*23                            :: binfile2
!
! constants
double precision, parameter             :: g = 9.81
!
! ----------------------------------------------------------
!do itime = 2, 25
!
dirname = '/glade/u/home/fdemartin/'
dirname2 = 'pp.v-50u50'
!
! Read in all binary files 
! ------------------------
!
! Retrieve forcing term

binfile = 'dynamic.bin'
!
print*,'reading forcing'
open(unit=10, file = binfile, form = 'unformatted', &
    access = 'stream')
print*,' '

irec = 1
do k = 1,nz
   read(10) lappd(:,:,k)
   irec = irec + 1
enddo
close(10)

binfile = 'forcing.bin'
!
print*,'reading forcing'
open(unit=10, file = binfile, form = 'unformatted', &
    access = 'stream')
print*,' '

irec = 1
do k = 1,nz
   read(10) lappb(:,:,k)
   irec = irec + 1
enddo
close(10)

print*,maxval(lappb)

binfile = 'exner.bin'

print*,'reading exner'
open(unit=10, file = binfile, form = 'unformatted', &
    access = 'stream')
print*,' '

irec = 1
do k = 1,nz
   read(10) pi(:,:,k)
   irec = irec + 1
enddo
close(10)

binfile = 'exner0.bin'

print*,'reading exner0'
open(unit=10, file = binfile, form = 'unformatted', &
    access = 'stream')
print*,' '

irec = 1
do k = 1,nz
   read(10) pi0(:,:,k)
   irec = irec + 1
enddo
close(10)

binfile = 'z.bin'

print*,'reading z'
open(unit=10, file = binfile, form = 'unformatted', &
    access = 'stream')
print*,' '

irec = 1
do k = 1,nz
   read(10) z(:,:,k)
   irec = irec + 1
enddo
close(10)

!do k=1,nz
!        print*,z(2,2,k),k
!end do
!
! ----------------------------------------------------
! Diagnosis algorithm
! ----------------------------------------------------
!
print*,''
print*,'---------------------------------'
print*,' Computing pressure perturbations'
print*,'---------------------------------'
!
! Assign first guess for Laplacians
!
print*,'Assuming first guess for pressure perturbations'
print*,''
do i = 1,nx
   do j = 1,ny
      do k = 1,nz
         pd(i,j,k) = 0.0
         pb(i,j,k) = 0.0
         !if (k.eq.1) then
         !        pb(i,j,k)=100.
         !endif
      enddo
   enddo
enddo
!

open(unit=10, file='outputp.bin', form='unformatted', access='stream')
    write(10) pb
    close(10)

print*, ' '
print*, '--------------------------------'
print*, 'Gauss-Seidel iterations begining'
print*, '--------------------------------'
print*, ' '
itnum = 1
diff = error + 0.1
! ----------------------------------------------------------------
! Iterations
! ----------------------------------------------------------------
do while ((diff .gt. error).and.(itnum .le. 4000))
   print*,'iteration no. ', itnum, ' and diff = ',diff
   diff = error * 0.1
   do i = 2, nx-1
      do j = 2,ny-1
         do k = 2,nz-1
            dzl = (z(i,j,k+1) - z(i,j,k))
            dzr = (z(i,j,k) - z(i,j,k-1))
            !print*, z(i,j,k+1),z(i,j,k-1),k-1,k+1,i,j
            dz = 0.5*(z(i,j,k+1) - z(i,j,k-1))
            !print*,dz
            oldpd = pd(i,j,k)
            call gauss(pd(i,j,k),lappd(i,j,k),pd(i+1,j,k), &
                       pd(i-1,j,k),pd(i,j+1,k),pd(i,j-1,k), &
                       pd(i,j,k+1),pd(i,j,k-1), &
                       dx,dzl,dzr,dz)
            if (itnum .eq. 1) then
                   pb(i,j,k) = pi(i,j,k)**(3.5)*100000 - pi0(i,j,k)**(3.5)*100000 - &
                               pd(i,j,k)
            endif
            oldpb = pb(i,j,k)
            call gauss(pb(i,j,k),lappb(i,j,k),pb(i+1,j,k), &
                       pb(i-1,j,k),pb(i,j+1,k),pb(i,j-1,k), &
                       pb(i,j,k+1),pb(i,j,k-1), &
                       dx,dzl,dzr,dz)
            diff = max(diff,abs(oldpd-pd(i,j,k)))
            diff = max(diff,abs(oldpb-pb(i,j,k)))
         enddo
      enddo
   enddo
   !
   ! This adds zero-gradient boundary condition for dynamic pressure
   ! at the bottom keeping with Wilhelmson and Ogura (1972)
   !
   do i = 1,nx
      do j = 1,ny
         pb(i,j,1) = pi(i,j,1)**(3.5)*100000 - pi0(i,j,1)**(3.5)*100000 - pd(i,j,1)
      enddo
   enddo
   !
   ! Because convection is restricted to the center of our domain,
   !  we assume a zero gradient condition on the lateral boundary.
   !
   do i = 2, nx-1
      do k = 2,nz-1
        pd(i,1,k) = pd(i,2,k)   
        pb(i,1,k) = pb(i,2,k)
      enddo
   enddo
   do k =2,nz-1
      do j = 2,ny-1
        pd(1,j,k) = pd(2,j,k) 
        pb(1,j,k) = pb(2,j,k)
      enddo
   enddo
   !
   ! End iterations
   !
   itnum = itnum + 1
enddo
!
! ------------------------------------------------------
! Compute vertical derivatives of pressure perturbations
! ------------------------------------------------------
!
do i = 1, nx
   do j = 1, ny
        dpb(i,j,1) = -1e-3 * thv0 * &
             (pb(i,j,2)-pb(i,j,1))/(z(i, j, 2)-z(i, j, 1))
        dpd(i,j,1) = -1e-3*thv0*(pd(i,j,2) - pd(i,j,1)) / &
             (z(i,j,2)-z(i,j,1))
        do k = 2, nz
           dz = 1e3*(z(i,j,k+1)-z(i,j,k-1))
           thv0 = theta0(i, j, k) * (1.0 + 0.61 * qv(i, j, k))
           dpb(i,j,k) = -thv0 * (pb(i,j,k+1)-pb(i,j,k-1))/dz
           dpd(i,j,k) = -thv0 * (pd(i,j,k+1)-pd(i,j,k-1))/dz
        enddo
        thv0 = theta0(i, j, nz) * (1.0 + 0.61 * qv(i, j, nz))
        dpb(i,j,nz) = -1e-3 * thv0 * &
                      (pb(i,j,nz)-pb(i,j,nz-1)) / &
                      (z(i,j,nz)-z(i,j,nz-1))
        dpd(i,j,nz) = -1e-3 * thv0 * &
                    (pd(i,j,nz)-pd(i,j,nz-1)) / &
                    (z(i,j,nz)-z(i,j,nz-1))
    enddo
enddo
!
! ---------------------------------------------------------------
! Write out output
!
print*, ' '
print*, '-------------------------'
print*, 'Writing out final output'
!


open(unit=10, file='pb.bin', form='unformatted', access='stream')
    write(10) pb
    close(10)

open(unit=10, file='pd.bin', form='unformatted', access='stream')
    write(10) pd
    close(10)


!
end
!
subroutine gauss(val,forc,pip1,pim1,pjp1,pjm1, &
                pkp1,pkm1,dx,dzl,dzr,dz)
! ===============================================================
! Mass field inversion
! ===============================================================
double precision    :: pip1,pim1,pjp1,pjm1,pkp1,pkm1
double precision    :: Dijk,dx,dzl,dzr,dz,th0v
double precision    :: val, forc
!print*,forc,dx,dy,dzr,dzl,dz,pip1,pim1,pjp1,pjm1,pkp1,pkm1
val = (-forc*dx*dx*dzr*dzl*dz &
       + (pip1+pim1+pjp1+pjm1)*dz*dzr*dzl &
       +2.*pkp1*dzr*dx*dx+2.*pkm1*dzl*dx*dx &
       )/(4.0*dzr*dzl*dz+2.*(dzr+dzl)*dx*dx)
!print*,val
!
return
end subroutine gauss
