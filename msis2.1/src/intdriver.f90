program msis2int_driver

  use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
  use msis_init, only: msisinit
  use netcdf

  implicit none (type, external)

  interface
    subroutine ghp8(day, utsec, z0, glat, glon, f107a, f107, ap, pres, alt, dn, tn)
      real, intent(in) :: day
      real, intent(in) :: utsec
      real, intent(in) :: z0
      real, intent(in) :: glat,glon
      real, intent(in) :: f107a,f107
      real, intent(in) :: ap(7)
      real, intent(in) :: pres
      real, intent(out) :: alt
      real, intent(out) :: dn(10)
      real, intent(out) :: tn
    end subroutine ghp8
    subroutine read_nc(fname, nx, ny, np, nt, lon, lat, z, t, o, o2, n2)
      character(len=*),                intent(in)  :: fname
      integer,                         intent(in)  :: nx, ny, np, nt
      real, dimension(nx),             intent(out) :: lon
      real, dimension(ny),             intent(out) :: lat
      real, dimension(nx, ny, np, nt), intent(out) :: z, t, o, o2, n2
    end subroutine read_nc
    subroutine write_nc(fname, nx, ny, nz, nt, lon, lat, hlevs, outdata, outhe, outo, outo2, outn2)
      character(len=*),                intent(in) :: fname
      integer,                         intent(in) :: nx, ny, nz, nt
      real, dimension(nx),             intent(in) :: lon
      real, dimension(ny),             intent(in) :: lat
      real, dimension(nz),             intent(in) :: hlevs
      real, dimension(nx, ny, nz, nt), intent(in) :: outdata, outhe, outo, outo2, outn2
    end subroutine write_nc
    subroutine check(istatus)
      integer, intent(in) :: istatus
    end subroutine check
  end interface

  integer, parameter :: nt = 1   ! number of time points
!  real    :: dt    ! change per timestep  (seconds)

  integer :: idoy  ! day of year
  real    :: ihour ! hour of the day
  character(len=255) :: input_file
  character(len=255) :: output_file
  real, dimension(nt)   :: f107a, f107
  real, dimension(7,nt) :: ap
  namelist /fixed_height/ idoy, ihour, input_file, output_file, f107, f107a, ap

  integer, parameter :: np    = 150 ! number of vertical input levels
  integer, parameter :: npref = 140 ! reference level at which to begin helium
  real,    parameter :: ppref = 3.22980e-08 ! WAM reference level pressure, Pa

  integer, parameter :: nz = 91   ! number of vertical output points
  integer, parameter :: ny = 91   ! number of latitudinal points  (input/output)
  integer, parameter :: nx = 90   ! number of longitudinal points (input/output)

  real, parameter    :: base_height = 100. ! output bottom vertical layer (km)

  real, parameter    :: odz = 10.        ! change per vertical layer in output (km)
  real, parameter    :: dy = 180./(ny-1) ! change per latitude  (degrees)
  real, parameter    :: dx = 360./nx     ! change per longitude (degrees)

  real, parameter :: Re = 6371008.7714 ! Earth radius, m
  real, parameter :: o_m  = 16.0    ! atomic oxygen atomic mass, u
  real, parameter :: o2_m = 32.0    ! molecular oxygen atomic mass, u
  real, parameter :: n2_m = 28.01   ! molecular nitrogen atomic mass, u
  real, parameter :: he_m = 4.0026  ! helium atomic mass, u
  real, parameter :: amu = 1.66e-27 ! atomic mass unit, kg
  real, parameter :: g = 9.80665    ! gravity, m/s^2
  real, parameter :: boltz = 1.38062e-23 ! Boltzmann constant, m^2*kg/(K*s^2)

  real, dimension(nx) :: lon
  real, dimension(ny) :: lat
  real, dimension(nz) :: hlevs, th
  real, dimension(:), allocatable :: utsec
  real, dimension(:), allocatable :: doy

  real, dimension(:,:,:), allocatable :: den_pref

  real :: zkm, dz, heh_pref, g_pref, thref

  character(255) :: argv

  real :: switch_legacy(1:25)

  integer :: ix, iy, iz, it, ip, izd

  ! msis variables
  real :: den(10)
  real :: temp      ! output that is ignored

  ! input/output variables
  real, dimension(:,:,:,:), allocatable :: dnh, z, t, o, o2, n2, oh, heh, o2h, n2h

  ! namelist init
  read(unit=5, nml=fixed_height)

  ! msis initialization
  switch_legacy(:) = 1   ! turn on all switches
  switch_legacy(9) = -1  ! turn on storm mode
  call msisinit(parmfile="msis21.parm", switch_legacy=switch_legacy)

  allocate(z(nx, ny, np, nt), t(nx, ny, np, nt), o(nx, ny, np, nt), &
           o2(nx, ny, np, nt), n2(nx, ny, np, nt))
  allocate(dnh(nx, ny, nz, nt), heh(nx, ny, nz, nt), oh(nx, ny, nz, nt), o2h(nx, ny, nz, nt), &
           n2h(nx, ny, nz, nt))
  allocate(den_pref(nx, ny, nt))
  allocate(doy(nt), utsec(nt))

  call read_nc(trim(input_file), nx, ny, np, nt, lon, lat, z, t, o, o2, n2)

  den_pref = o(:,:,npref,:) + o2(:,:,npref,:) + n2(:,:,npref,:)
  o( :,2:ny-1,:,:) = log(o( :,2:ny-1,:,:))
  o2(:,2:ny-1,:,:) = log(o2(:,2:ny-1,:,:))
  n2(:,2:ny-1,:,:) = log(n2(:,2:ny-1,:,:))

  ! --- execute program

  doy(1)   = real(idoy)
  utsec(1) = ihour*3600

!  do it = 2, nt
!    doy(it) = doy(it-1)
!    utsec(it) = modulo(ihour*3600 + (it-1) * dt, 24.*3600)
!    if (utsec(it) .lt. utsec(it-1)) then
!      doy(it) = doy(it) + 1
!    end if
!  end do
  do iz = 1, nz
    hlevs(iz) = (base_height + (iz-1) * odz) * 1000
  end do

  do it = 1, nt
    do iy = 2, ny-1
      do ix = 1, nx
        do iz = 1, nz
          ip = 1
          ! find pressure level at which plev is below hlev and plev+1 is above hlev
          do while(ip .le. np-1 .and. &
                   .not. (hlevs(iz) .ge. z(ix,iy,ip,it) .and. &
                          hlevs(iz) .le. z(ix,iy,ip+1,it)) )
            ip = ip + 1
            if (ip .eq. np) exit
          end do

          ! once hlev is above model top, no point in continuing in this loop
          if (ip .eq. np) exit

          ! (when appropriate, log-)interpolate values between pressure levels to fixed height grid
          dz = (hlevs(iz) - z(ix, iy, ip, it)) / (z(ix, iy, ip+1, it) - z(ix, iy, ip, it))
          th(iz)  =     dz *  t(ix, iy, ip+1, it) + (1-dz) *  t(ix, iy, ip, it)
          oh(ix, iy, iz, it)  = exp(dz *  o(ix, iy, ip+1, it) + (1-dz) *  o(ix, iy, ip, it))
          o2h(ix, iy, iz, it) = exp(dz * o2(ix, iy, ip+1, it) + (1-dz) * o2(ix, iy, ip, it))
          n2h(ix, iy, iz, it) = exp(dz * n2(ix, iy, ip+1, it) + (1-dz) * n2(ix, iy, ip, it))
          dnh(ix, iy, iz, it) = n2h(ix, iy, iz, it)*n2_m + o2h(ix, iy, iz, it)*o2_m + oh(ix, iy, iz, it)*o_m
        end do ! iz

        iz = 1
        ! find hlev such that hlev is just above the reference level
        do while(iz.lt.nz .and. hlevs(iz) .lt. z(ix, iy, npref, it))
          iz = iz + 1
        end do

        ! store iz as reference for downward extrapolation later
        izd = iz

        ! get reference level helium
        call ghp8(doy(it), utsec(it), z(ix, iy, npref, it)/1000, lat(iy), lon(ix), f107a(it), &
                  f107(it), ap(1:7,it), ppref, zkm, den, temp)

        heh_pref = den(5)/(den(2)+den(3)+den(4)) * den_pref(ix, iy, it)
        g_pref   = g*(re/(re+z(ix, iy, npref, it)))**2
        ! extrapolate up to hlev
        dz = z(ix, iy, npref, it) - hlevs(iz)
        heh(ix, iy, iz, it) = heh_pref * exp(dz*he_m*amu*g_pref/(boltz*t(ix, iy, npref, it)))
        dnh(ix, iy, iz, it) = dnh(ix, iy, iz, it) + heh(ix, iy, iz, it)*he_m

        ! sweep to model top, log extrapolating upwards, using previously calculated
        ! o, o2, n2 values below model top
        dz = -odz * 1000.
        do while(hlevs(iz) .le. z(ix, iy, np-1, it))
          g_pref = g*(re/(re+hlevs(iz)))**2

          heh(ix, iy, iz+1, it) = heh(ix, iy, iz, it) * exp(dz*he_m*amu*g_pref/(boltz*th(iz)))
          dnh(ix, iy, iz+1, it) = dnh(ix, iy, iz+1, it) + heh(ix, iy, iz+1, it)*he_m

          iz = iz + 1
        end do

        thref = th(iz-1)
        ! sweep from model top to output top, log extrapolating all values upwards
        do iz = iz, nz
          g_pref = g*(re/(re+hlevs(iz-1)))**2

          heh(ix, iy, iz, it) = heh(ix, iy, iz-1, it) * exp(dz*he_m*amu*g_pref/(boltz*thref))
          oh( ix, iy, iz, it) = oh( ix, iy, iz-1, it) * exp(dz* o_m*amu*g_pref/(boltz*thref))
          o2h(ix, iy, iz, it) = o2h(ix, iy, iz-1, it) * exp(dz*o2_m*amu*g_pref/(boltz*thref))
          n2h(ix, iy, iz, it) = n2h(ix, iy, iz-1, it) * exp(dz*n2_m*amu*g_pref/(boltz*thref))
          dnh(ix, iy, iz, it) = heh(ix, iy, iz, it)*he_m + n2h(ix ,iy, iz, it)*n2_m + o2h(ix, iy, iz, it)*o2_m + &
                                 oh(ix, iy, iz, it)*o_m
        end do

        ! sweep from reference level to base height, log extrapolating downwards
        ! and using previous values for o, o2, n2
        dz = odz * 1000.
        do iz = izd, 2, -1
          g_pref = g*(re/(re+hlevs(iz)))**2

          heh(ix, iy, iz-1, it) = heh(ix, iy, iz, it) * exp(dz*he_m*amu*g_pref/(boltz*th(iz)))
          dnh(ix, iy, iz-1, it) = dnh(ix, iy, iz-1, it) + heh(ix, iy, iz-1, it)*he_m
        end do
      end do !ix
    end do !iy
  end do !it

  dnh = dnh * amu

  call write_nc(trim(output_file), nx, ny, nz, nt, lon, lat, hlevs, dnh, heh, oh, o2h, n2h)

end program


subroutine ghp8(day, utsec, z0, glat, glon, f107a, f107, ap, pres, alt, dn, tn)

  use msis_constants, only: kB,Na,g0,rp
  use msis_init, only          : msisinit
  use msis_calc, only          : msiscalc
  use msis_utils, only         : alt2gph

  implicit none

  real(kind=rp),intent(in)    :: day
  real(kind=rp),intent(in)    :: utsec
  real(kind=rp),intent(in)    :: z0       !! first guess
  real(kind=rp),intent(in)    :: glat,glon
  real(kind=rp),intent(in)    :: f107a,f107
  real(kind=rp),intent(in)    :: ap(7)
  real(kind=rp),intent(in)    :: pres   !!! pressure in hPa
  real(kind=rp),intent(out)   :: alt
  real(kind=rp),intent(out)   :: dn(10)   !!! # density are now in different order then MSIS00
  real(kind=rp),intent(out)   :: tn
! Local variables
!  real(8), external           :: alt2gph
  real(kind=rp)               :: tex
  real(kind=rp)               :: plog,delta
  real(kind=rp)               :: zkm,pzkm
  real(kind=rp)               :: xn,gz,xmbar,scl
  integer                     :: n
  real(kind=rp),parameter     :: tol = 0.000043_rp
  integer,parameter           :: maxit = 30
  real(8)                     :: xlat,alt0,alt1

  plog = log10(pres*100.0_rp)
  zkm = z0
  delta = 1.0_rp

  n = 0

  do while ( abs(delta) .ge. tol .and. n .le. maxit)
    n = n + 1

    call msiscalc(day,utsec,zkm,glat,glon,f107a,f107,ap,tn,dn,tex)

    xn = sum(dn(2:8))
    pzkm = kB * xn * tn
    delta = plog - log10(pzkm)
    xmbar = dn(1) / xn / 1.66E-24_rp
    xlat = dble(glat)
    alt0 = dble(zkm)
    alt1 = alt0 + 1.0d0
    gz = real((alt2gph(xlat,alt1) - alt2gph(xlat,alt0)) * g0)
    scl = Na * kB * tn / (xmbar * gz)

    ! difference
    zkm = zkm - scl * delta / 1000.0_rp
  end do
  alt = zkm

end subroutine ghp8


subroutine read_nc(fname, nx, ny, np, nt, lon, lat, z, t, o, o2, n2)

  use netcdf

  implicit none

  character(len=*),                intent(in)  :: fname
  integer,                         intent(in)  :: nx, ny, np, nt
  real, dimension(nx),             intent(out) :: lon
  real, dimension(ny),             intent(out) :: lat
  real, dimension(nx, ny, np, nt), intent(out) :: z, t, o, o2, n2
! Local variables
  integer ncid, varid

  call check(nf90_open(trim(fname), nf90_nowrite, ncid))

  call check(nf90_inq_varid(ncid, "lon",          varid))
  call check(nf90_get_var(ncid, varid, lon))

  call check(nf90_inq_varid(ncid, "lat",          varid))
  call check(nf90_get_var(ncid, varid, lat))

  call check(nf90_inq_varid(ncid, "height",       varid))
  call check(nf90_get_var(ncid, varid, z))

  call check(nf90_inq_varid(ncid, "temp_neutral", varid))
  call check(nf90_get_var(ncid, varid, t))

  call check(nf90_inq_varid(ncid, "O_Density",    varid))
  call check(nf90_get_var(ncid, varid, o))

  call check(nf90_inq_varid(ncid, "O2_Density",   varid))
  call check(nf90_get_var(ncid, varid, o2))

  call check(nf90_inq_varid(ncid, "N2_Density",   varid))
  call check(nf90_get_var(ncid, varid, n2))

  call check(nf90_close(ncid))

end subroutine read_nc


subroutine write_nc(fname, nx, ny, nz, nt, lon, lat, hlevs, outdata, outhe, outo, outo2, outn2)

  use netcdf

  implicit none

  character(len=*),                intent(in) :: fname
  integer,                         intent(in) :: nx, ny, nz, nt
  real, dimension(nx),             intent(in) :: lon
  real, dimension(ny),             intent(in) :: lat
  real, dimension(nz),             intent(in) :: hlevs
  real, dimension(nx, ny, nz, nt), intent(in) :: outdata, outhe, outo, outo2, outn2
! Local variables
  integer :: ncid, x_dimid, x_varid, y_dimid, y_varid, z_dimid, z_varid, t_dimid, t_varid, den_varid
  integer :: he_varid, o_varid, o2_varid, n2_varid
  integer, dimension(4) :: dimids

  call check(nf90_create(trim(fname), nf90_clobber, ncid))

  call check(nf90_def_dim(ncid, "lon",   nx, x_dimid))
  call check(nf90_def_dim(ncid, "lat",   ny, y_dimid))
  call check(nf90_def_dim(ncid, "hlevs", nz, z_dimid))
  call check(nf90_def_dim(ncid, "time",  NF90_UNLIMITED, t_dimid))

  call check(nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, x_varid))
  call check(nf90_put_att(ncid, x_varid, "axis","X"))
  call check(nf90_put_att(ncid, x_varid, "long_name","longitude"))
  call check(nf90_put_att(ncid, x_varid, "units","degrees_east"))
  call check(nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, y_varid))
  call check(nf90_put_att(ncid, y_varid, "axis","Y"))
  call check(nf90_put_att(ncid, y_varid, "long_name","latitude"))
  call check(nf90_put_att(ncid, y_varid, "units","degrees_north"))
  call check(nf90_def_var(ncid, "hlevs", NF90_REAL, z_dimid, z_varid))
  call check(nf90_put_att(ncid, z_varid, "axis","Z"))
  call check(nf90_put_att(ncid, z_varid, "long_name","altitude"))
  call check(nf90_put_att(ncid, z_varid, "units","km"))
  call check(nf90_def_var(ncid, "time", NF90_DOUBLE, t_dimid, t_varid))
  call check(nf90_put_att(ncid, t_varid, "axis","T"))
  call check(nf90_put_att(ncid, t_varid, "long_name","time"))
  call check(nf90_put_att(ncid, t_varid, "units","days since 1970-01-01"))

  dimids = (/ x_dimid, y_dimid, z_dimid, t_dimid /)

  call check(nf90_def_var(ncid, "den", NF90_FLOAT, dimids, den_varid))
  call check(nf90_put_att(ncid, den_varid, "units","m^-3"))
  call check(nf90_put_att(ncid, den_varid, "standard_name","neutral density"))

!  call check(nf90_def_var(ncid, "he",  NF90_FLOAT, dimids, he_varid))
!  call check(nf90_put_att(ncid, den_varid, "units","m^-3"))
!  call check(nf90_put_att(ncid, den_varid, "standard_name","neutral density"))

!  call check(nf90_def_var(ncid, "o",   NF90_FLOAT, dimids, o_varid))
!  call check(nf90_put_att(ncid, den_varid, "units","m^-3"))
!  call check(nf90_put_att(ncid, den_varid, "standard_name","neutral density"))

!  call check(nf90_def_var(ncid, "o2",  NF90_FLOAT, dimids, o2_varid))
!  call check(nf90_put_att(ncid, den_varid, "units","m^-3"))
!  call check(nf90_put_att(ncid, den_varid, "standard_name","neutral density"))

!  call check(nf90_def_var(ncid, "n2",  NF90_FLOAT, dimids, n2_varid))
!  call check(nf90_put_att(ncid, den_varid, "units","m^-3"))
!  call check(nf90_put_att(ncid, den_varid, "standard_name","neutral density"))

  call check(nf90_enddef(ncid))

  call check(nf90_put_var(ncid, x_varid, lon))
  call check(nf90_put_var(ncid, y_varid, lat))
  call check(nf90_put_var(ncid, z_varid, hlevs / 1000))
  ! I am not smart enough to write a good function for adding time, so leave it to Python.
  call check(nf90_put_var(ncid, den_varid, outdata))
!  call check(nf90_put_var(ncid, he_varid,  outhe))
!  call check(nf90_put_var(ncid, o_varid,   outo))
!  call check(nf90_put_var(ncid, o2_varid,  outo2))
!  call check(nf90_put_var(ncid, n2_varid,  outn2))

  call check(nf90_close(ncid))

end subroutine write_nc


subroutine check(istatus)

  use netcdf

  implicit none

  integer, intent(in) :: istatus

  if (istatus /= nf90_noerr) then
    write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
  end if

end subroutine check
