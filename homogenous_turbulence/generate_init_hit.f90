!Code for generating initial field for homogeneous isotropic turbulence
!Adapted from M. Shoeybi's code

!compilation
!export PATH=/usr/bin:$PATH
!gfortran -I/usr/local/Cellar/fftw/3.3.7/include/ -o generate_init_hit generate_init_hit.f90 -L/usr/local/Cellar/fftw/3.3.7/lib -lfftw3

! gfortran -I/usr/include/ -o generate_init_hit generate_init_hit.f90 -lfftw3

!ifort -I$WORKDIR/SRC/fftw/include -o generate_init_hit generate_init_hit.f90 -L$WORKDIR/SRC/fftw/lib -lfftw3

module rand_tools_m

contains

  subroutine randset(idum)
    
    implicit none

    integer :: idum

    real(8) :: r(97)
    integer :: ix1,ix2,ix3
    common /random_common/ r,ix1,ix2,ix3

    real(8) :: rm1,rm2
    integer :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
    parameter(m1=259200,ia1=7141,ic1=54773,rm1=1.d0/m1)
    parameter(m2=134456,ia2=8121,ic2=28441,rm2=1.d0/m2)
    parameter(m3=243000,ia3=4561,ic3=51349)

    integer :: j

    ix1 = mod(ic1 - idum , m1)
    ix1 = mod(ia1 * ix1 + ic1 , m1)
    ix2 = mod(ix1 , m2)
    ix1 = mod(ia1 * ix1 + ic1 , m1)
    ix3 = mod(ix1 , m3)

    do j = 1, 97
       ix1 = mod(ia1 * ix1 + ic1 , m1)
       ix2 = mod(ia2 * ix2 + ic2 , m2)
       r(j) = (dble(ix1) + dble(ix2) * rm2) * rm1
       if (r(j).gt.1.d0) r(j)=1.d0
       if (r(j).lt.0.d0) r(j)=0.d0
    enddo

  end subroutine randset


  subroutine randvec(ran,N)

    implicit none

    integer :: N
    real(8) :: ran(N)
    real(8) :: r(97)
    integer :: ix1,ix2,ix3
    common /random_common/ r,ix1,ix2,ix3

    real(8) :: rm1,rm2
    integer :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
    parameter(m1=259200,ia1=7141,ic1=54773,rm1=1.d0/m1)
    parameter(m2=134456,ia2=8121,ic2=28441,rm2=1.d0/m2)
    parameter(m3=243000,ia3=4561,ic3=51349)

    integer :: nombre,j

    do nombre=1,N
       ix1 = mod(ia1 * ix1 + ic1 , m1)
       ix2 = mod(ia2 * ix2 + ic2 , m2)
       ix3 = mod(ia3 * ix3 + ic3 , m3)

       j = 1 + mod(((97 * ix3) / m3),97)
       ran(nombre) = r(j)
       r(j) = (dble(ix1) + dble(ix2) * rm2) * rm1
       if (r(j).gt.1.d0) r(j)=1.d0
       if (r(j).lt.0.d0) r(j)=0.d0
    enddo

  end subroutine randvec  

end module rand_tools_m




program generate_init_hit

  use rand_tools_m

  implicit none

  !include fftw header file
  include "fftw3.f"

  integer, parameter :: N_box = 129

  real(8), parameter :: k0 = 12.0_8

  integer, parameter :: divergence_flag = 0
  !            0: use the noraml construction which gives 
  !               the divergecne free for spectral codes
  !            1: use modified wave numbers
  
  !constants
  real(8) :: pi, twopi

  ! grid information
  real(8) :: my_node_cc_min(3), my_node_cc_max(3), node_cc_min(3), node_cc_max(3)
  real(8) :: my_box_volume, box_volume
  real(8) :: length(3), dx(3)
  integer  :: N_total, num_node

  ! initial condition in Fourier space
  real(8) :: u0
  real(8) :: kx, ky, kz, kxyz, kmax, tolerance, kmin, kxy, kmax_m
  integer :: kx_max_i, ky_max_i, kz_max_i

  ! energy parameters
  real(8) :: Eq_n 
  complex(8) :: a_ran, b_ran
  complex, parameter :: ima = (0.0_8,1.0_8)

  ! random number
  real(8), allocatable :: phi1(:,:,:), phi2(:,:,:), phi3(:,:,:)
  real(8) :: max_phi12, max_phi13, max_phi23, min_dphi

  ! Fourier ans velocity vectors
  real(8), allocatable :: u(:,:,:), u3(:,:,:,:), u3_ext(:,:,:,:)
  complex(8), allocatable :: uk(:,:,:)!,uk_div(:,:,:,:)
  integer*8 :: fftw_plan 

  ! debug stuff
  real(8) :: my_max_u, my_min_u, my_mean_u(3)
  real(8) :: max_u, min_u, mean_u(3), max_divergence
  real(8) :: u_left(1:3), u_right(1:3)
  real(8), allocatable :: div(:,:,:)
     
  ! loop parameters
  integer :: i, j, k, i_u, ino, ierr, n_phi

  ! dump out file
  character(len=64) :: filename
  !-------------------------------------    

  write(*,*) '     Generate initial condition for HIT'

  pi = 4.0_8 * atan(1.0_8)
  twopi = 2.0_8 * pi

  if ( divergence_flag.eq.0 ) then
     write(*,*) '       using spectral divergence free field'
  elseif ( divergence_flag.eq.1 ) then
     write(*,*) '       using modified wave number to get finite differecne divergence free field'
  else
     write(*,*) 'Error in in calc_init_for_hom_iso_tur: undefined divergence_flag'
     stop
  endif

  ! total number of grid points
  N_total = N_box * N_box * N_box
  
   
  ! check even number of points has been passed
  num_node = N_box - 1
  if ( mod(num_node,2).ne.0 ) then
     write(*,*) 'Error in calc_init_for_hom_iso_tur, pass odd number of points '
     stop
  endif
    
  ! get the grid information
  ! box length  
  length(1:3) = twopi

  ! grid spacing
  dx(1:3) = length(1:3) / dble(N_box-1)
    
  ! write some data on the screen
  write(*,*) '*****************************************************************************'
  write(*,*) '  grid information'
  write(*,*) '     total number of grid points(cart)               =  ', N_total
  write(*,*) '     number of grid points in each direction(cart)   =  ', N_box
  write(*,*) '     length in x, y and z direction (should be 2*pi) =  ', length(1:3)
  write(*,*) '     grid spacing in x, y and z(cart)                =  ', dx(1:3)
  write(*,*) '*****************************************************************************'

  ! spectrum paprameters
  tolerance = 0.0001_8
  kmax = dble(N_box-1) / 2.0_8 - tolerance
  kmin = tolerance
  kx_max_i = idint( dble(N_box-1) / 2.0_8 + 0.1_8 ) + 1
  ky_max_i = N_box-1
  kz_max_i = N_box-1

  ! spectrum information
  ! write on the screen
  write(*,*) '*****************************************************************************'
  write(*,*) 'spectrum properties'
  write(*,*) '   wave number corresponding to most energetic wave =  ', k0
  write(*,*) '*****************************************************************************'      
  write(*,*) 'wave number information'
  write(*,*) '   maximun wave number(for homogeniety) =  ', kmax
  write(*,*) '   min wave number                      =  ', kmin
  write(*,*) '   max wave number in x direction       = ',kx_max_i
  write(*,*) '   max wave number in y direction       = ',ky_max_i
  write(*,*) '   max wave number in z direction       = ',kz_max_i
  write(*,*) '*****************************************************************************'      

  ! random numbers
  ! allocate memory
  allocate( phi1(kx_max_i,ky_max_i,kz_max_i) )
  allocate( phi2(kx_max_i,ky_max_i,kz_max_i) )
  allocate( phi3(kx_max_i,ky_max_i,kz_max_i) )

  ! call random number generator
  call randset(0)
  call randvec(phi1,kx_max_i*ky_max_i*kz_max_i)
  call randvec(phi2,kx_max_i*ky_max_i*kz_max_i)
  call randvec(phi3,kx_max_i*ky_max_i*kz_max_i)

  ! check the random numbers
  max_phi12 = maxval( abs( phi1(1:kx_max_i,1:ky_max_i,1:kz_max_i) -&
       phi2(1:kx_max_i,1:ky_max_i,1:kz_max_i) ) )
  max_phi13 = maxval( abs( phi1(1:kx_max_i,1:ky_max_i,1:kz_max_i) -&
       phi3(1:kx_max_i,1:ky_max_i,1:kz_max_i) ) )
  max_phi23 = maxval( abs( phi3(1:kx_max_i,1:ky_max_i,1:kz_max_i) -&
       phi2(1:kx_max_i,1:ky_max_i,1:kz_max_i) ) )
  
  min_dphi = min(max_phi12,max_phi13)
  min_dphi = min(min_dphi,max_phi23)

  if ( min_dphi.le.tolerance ) then
     write(*,*) 'Error in calc_init_cond_for_hom_iso_tu: random numbers are the same'
     stop
  endif
       
  ! allocate memory
  allocate( uk(kx_max_i,ky_max_i,kz_max_i) )
  allocate( u(ky_max_i,ky_max_i,kz_max_i) )
  allocate( u3(3,ky_max_i,ky_max_i,kz_max_i) )
  ! allocate( uk_div(3,kx_max_i,ky_max_i,kz_max_i) )

  ! modify the random numbers to get them in the 
  ! interval(-pi,pi)
  do k = 1,kz_max_i
     do j = 1,ky_max_i
        do i = 1,kx_max_i
           phi1(i,j,k) = twopi * (phi1(i,j,k)-0.5_8)
           phi2(i,j,k) = twopi * (phi2(i,j,k)-0.5_8)
           phi3(i,j,k) = twopi * (phi3(i,j,k)-0.5_8)
        enddo
     enddo
  enddo

  ! main loop for the fourier transforms
  do i_u = 1,3
     uk(1:kx_max_i,1:ky_max_i,1:kz_max_i) = ( 0.0_8, 0.0_8)
     u(1:ky_max_i,1:ky_max_i,1:kz_max_i) = 0.0_8
     do k = 1,kz_max_i
        do j = 1,ky_max_i
           do i = 1,kx_max_i
                
              ! modify the waves to get them symmetric
              kx = dble(i-1)
              ky = dble(j-1)
              if ( j.gt.kx_max_i ) ky = -dble(ky_max_i-j+1)
              kz = dble(k-1)
              if ( k.gt.kx_max_i ) kz = -dble(kz_max_i-k+1)

              kmax_m = dsqrt( kx**2 + ky**2 + kz**2 )

              kxyz = dsqrt( kx**2 + ky**2 + kz**2 )
              kxy = dsqrt( kx**2 + ky**2 )                

              ! enegry spectrum
              Eq_n = (kxyz/k0)**4 * dexp( -2.0_8*(kxyz/k0)**2 )
              
              a_ran = (0.0_8,0.0_8)
              b_ran = (0.0_8,0.0_8)
              
              if ( kxyz.gt.kmin .and. kxyz.le.kmax_m ) then
                 a_ran = dsqrt( Eq_n/(twopi* kxyz**2) ) *&
                      cdexp(ima*phi1(i,j,k)) * dcos(phi3(i,j,k))
                 b_ran = dsqrt( Eq_n/(twopi* kxyz**2) ) *&
                      cdexp(ima*phi2(i,j,k)) * dsin(phi3(i,j,k))
              endif
              
              if ( divergence_flag.eq.1 ) then
                 ! replace with modefied wave number
                 kx = sin( kx*twopi/dble(ky_max_i) ) / (twopi/dble(ky_max_i))
                 ky = sin( ky*twopi/dble(ky_max_i) ) / (twopi/dble(ky_max_i))
                 kz = sin( kz*twopi/dble(ky_max_i) ) / (twopi/dble(ky_max_i))
              endif
                
              kxyz = dsqrt( kx**2 + ky**2 + kz**2 )
              kxy = dsqrt( kx**2 + ky**2 )
              
              ! fourier coefficients
              uk(i,j,k) = (0.0_8,0.0_8)
              if ( kxyz.gt.kmin ) then
                 
                 if ( i_u.eq.1 ) then
                    
                    if ( kxy.gt.kmin ) then
                       uk(i,j,k) = (a_ran*kxyz*ky+b_ran*kx*kz)/(kxyz*kxy)
                    else
                       !uk(i,j,k) = a_ran
                       uk(i,j,k) = ( a_ran + b_ran ) / dsqrt(2.0_8)
                    endif
                    
                 elseif ( i_u.eq.2 ) then
                    
                    if ( kxy.gt.kmin ) then
                       uk(i,j,k) = (b_ran*ky*kz-a_ran*kxyz*kx)/(kxyz*kxy)
                    else
                       !uk(i,j,k) = b_ran
                       uk(i,j,k) = ( b_ran - a_ran ) / dsqrt(2.0_8)
                    endif

                 elseif ( i_u.eq.3 ) then

                    !if ( kxy.gt.kmin ) then
                    uk(i,j,k) = -b_ran * kxy / kxyz
                    !else
                    !   uk(i,j,k) = (0.0_8,0.0_8)
                    !endif
                    
                 endif
                 
              endif

           enddo
        enddo
     enddo

     ! at plane kx=0, the waves need to be complex conjugate
     do k = 2,kz_max_i
        do j = kx_max_i+1,ky_max_i
           uk(1,j,k) = dconjg(uk(1,ky_max_i-j+2,kz_max_i-k+2)) 
        enddo
     enddo
     do k = kx_max_i+1,kz_max_i
        uk(1,1,k) = dconjg(uk(1,1,kz_max_i-k+2))
     enddo
     do j = kx_max_i+1,ky_max_i
        uk(1,j,1) = dconjg(uk(1,ky_max_i-j+2,1))
     enddo
     
     ! remove odd-ball wave numbers
     uk(kx_max_i,1:ky_max_i,1:kz_max_i) = (0.0_8,0.0_8)
     uk(1:kx_max_i,kx_max_i,1:kz_max_i) = (0.0_8,0.0_8)
     uk(1:kx_max_i,1:ky_max_i,kx_max_i) = (0.0_8,0.0_8)
     
     ! remove mean
     uk(1,1,1) = (0.0_8,0.0_8)
     
     !take the fourier transform
     call dfftw_plan_dft_c2r_3d(fftw_plan, N_box-1, N_box-1, N_box-1, uk, u, 64)
     call dfftw_execute(fftw_plan)
     call dfftw_destroy_plan(fftw_plan)
     
     u3(i_u,1:N_box-1,1:N_box-1,1:N_box-1) = u(1:N_box-1,1:N_box-1,1:N_box-1)
     
  enddo
  
  ! deallocate memory
  deallocate(phi1); deallocate(phi2); deallocate(phi3)
  deallocate(uk); deallocate(u); 
  
  !add the periodic points to the grid   
  
  allocate( u3_ext(3,ky_max_i+1,ky_max_i+1,kz_max_i+1) )
  u3_ext(1:3,1:ky_max_i,1:ky_max_i,1:kz_max_i) = u3(1:3,1:ky_max_i,1:ky_max_i,1:kz_max_i)
  !periodic boundaries
  u3_ext(1:3,ky_max_i+1,1:ky_max_i,1:kz_max_i) = u3(1:3,1,1:ky_max_i,1:kz_max_i)
  u3_ext(1:3,1:ky_max_i,ky_max_i+1,1:kz_max_i) = u3(1:3,1:ky_max_i,1,1:kz_max_i)
  u3_ext(1:3,1:ky_max_i,1:ky_max_i,kz_max_i+1) = u3(1:3,1:ky_max_i,1:ky_max_i,1)
  
  u3_ext(1:3,ky_max_i+1,ky_max_i+1,1:kz_max_i) = u3(1:3,1,1,1:kz_max_i)
  u3_ext(1:3,ky_max_i+1,1:ky_max_i,kz_max_i+1) = u3(1:3,1,1:ky_max_i,1)
  u3_ext(1:3,1:ky_max_i,ky_max_i+1,kz_max_i+1) = u3(1:3,1:ky_max_i,1,1)
  u3_ext(1:3,ky_max_i+1,ky_max_i+1,kz_max_i+1) = u3(1:3,1,1,1)
  
  deallocate(u3)
  
  ! calculate divergence
  allocate(div(ky_max_i,ky_max_i,kz_max_i))

  do k = 1, kz_max_i
     do j = 1, ky_max_i
        do i = 1, kx_max_i

           if (i==1) then
              u_left(1) = u3_ext(1,ky_max_i,j,k)
           else
              u_left(1) = u3_ext(1,i-1,j,k)
           end if
           u_right(1) = u3_ext(1,i+1,j,k)

           if (j==1) then
              u_left(2) = u3_ext(2,i,ky_max_i,k)
           else
              u_left(2) = u3_ext(2,i,j-1,k)
           end if
           u_right(2) = u3_ext(2,i,j+1,k)

           if (k==1) then
              u_left(3) = u3_ext(3,i,j,kz_max_i)
           else
              u_left(3) = u3_ext(3,i,j,k-1)
           end if
           u_right(3) = u3_ext(3,i,j,k+1)

           div(i,j,k) = (u_right(1)-u_left(1))/2.0d0/dx(1) + &
                        (u_right(2)-u_left(2))/2.0d0/dx(2) + &
                        (u_right(3)-u_left(3))/2.0d0/dx(3)

        end do
     end do
  end do
  
  max_divergence = maxval(div)
  
  !dump the statistics on the screen
  write(*,*) '*****************************************************************************'
  write(*,*) 'initial condition'
  write(*,*) '   maximum mach number             =  ', max_u
  write(*,*) '   minimum mach number             =  ', min_u
  write(*,*) '   mean velocity                   =  ', mean_u(1:3)/box_volume
  write(*,*) '   mean velocity2                  =  ', mean_u(1:3)
  write(*,*) '   maximum divergence of velocity  =  ', max_divergence
  write(*,*) '*****************************************************************************'      
  

  !write the velocity field
  write(filename,'(a)') './init_field_hit.dat'
  open(10, file=filename, form='formatted', action='write')
  write(10,'(a)') 'VARIABLES = "X" "Y" "Z" "U" "V" "W"'
  write(10,'(a,i4.4,a,i4.4,a,i4.4)') "ZONE I=", ky_max_i+1, ", J=", ky_max_i+1, " K=", kz_max_i+1
  do k = 1, kz_max_i+1
     do j = 1, ky_max_i+1
        do i = 1, ky_max_i+1
           write(10,'(6e24.16)') dble(i-1)*dx(1), dble(j-1)*dx(2), dble(k-1)*dx(3), u3_ext(1:3,i,j,k)
        end do
     end do
  end do
  close(10)


  !deallocate memory
  deallocate(div)
  deallocate(u3_ext)
  
end program generate_init_hit
