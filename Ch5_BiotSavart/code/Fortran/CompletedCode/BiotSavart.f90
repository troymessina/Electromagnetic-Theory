! Example: Calculate the Magnetic Field around a Square Loop of Current-Carrying Wire with the Biot-Savart Law
! This code may be compiled with, for example, the gfortran compiler with
! > gfortran -o BiotSavart.x BiotSavart.f90
! Three configurations are set up in this file, via Fortran modules:
! - SquareLoop, which sets up the square loop of current-carrying wire.
! - CircleLoop, which sets up a circular loop of current-carrying wire.
! - StraightWire, which sets up a segment of straight current-carrying wire.
! In the module BiotSavartModule, the user may choose which wire to use
!    by uncommenting the module ``use`` statement corresponding to the
!    desired geometry, and commenting the other ``use`` statements.

module mesh
  implicit none
  ! The sizes of the calculation mesh
  integer, parameter :: NX=1, NY=50, NZ=50
  ! minimum and maximum values in the calculation space
  double precision, parameter :: minx=-2.0, maxx=2.0, &
       & miny = -2.0, maxy = 2.0, &
       & minz = -2.0, maxz = 2.0
  ! The spacing of the calculation mesh in each direction
  double precision, parameter :: dx=(maxx-minx)/NX, dy=(maxy-miny)/NY, dz=(maxz-minz)/NZ
  ! a tracker for the position in space
  double precision, dimension(3) :: r

  contains
  subroutine calc_position(ix, iy, iz)
    ! This subroutine takes a set of grid indices (ix,iy,iz) and converts it to a position,
    ! storing that value in the module variable r
    integer :: ix, iy, iz

    r(1) = 0.0d0 !minx + dx*ix
    r(2) = miny + dy*iy
    r(3) = minz + dz*iz
!!$    WRITE(*,*) "DEBUG: Calculating the position at r=", r(1), r(2), r(3)
  end subroutine calc_position

end module mesh

module Constants
  implicit none
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)
  double precision, parameter :: mu0 = 4.0d-7 * pi     ! T*m/A
  double precision, parameter :: I = 1.0d0             ! amps
end module Constants


module SquareLoop
  ! implement a square loop of current-carrying wire
  use Constants
  implicit none
  double precision, dimension(3, 4) :: vertices
  double precision, parameter :: length=1.0d0
  double precision, parameter :: field_at_origin = 2.0d0*sqrt(2.0d0)*mu0*I/(pi*length)

  contains

  subroutine setup()
    vertices(1,1) = 0.5d0*length
    vertices(2,1) = 0.5d0*length
    vertices(3,1) = 0.0d0

    vertices(1,2) = -0.5d0*length
    vertices(2,2) = 0.5d0*length
    vertices(3,2) = 0.0d0

    vertices(1,3) = -0.5d0*length
    vertices(2,3) = -0.5d0*length
    vertices(3,3) = 0.0d0

    vertices(1,4) = 0.5d0*length
    vertices(2,4) = -0.5d0*length
    vertices(3,4) = 0.0d0
  end subroutine setup

  subroutine find_dl(t, dl)
    double precision :: t
    double precision, dimension(3) :: dl
    integer :: vstart, vend

    if (t .lt. 0.25d0) then
       vstart = 1
       vend = 2
    else if (t .lt. 0.50d0) then
       vstart = 2
       vend = 3
    else if (t .lt. 0.75d0) then
       vstart = 3
       vend = 4
    else
       vstart = 4
       vend = 1
    end if

    ! each dl(i) = \frac{d rp(i)}{dt}
    dl(1) = 4.0d0*(vertices(1, vend) - vertices(1, vstart))
    dl(2) = 4.0d0*(vertices(2, vend) - vertices(2, vstart))
    dl(3) = 4.0d0*(vertices(3, vend) - vertices(3, vstart))
  end subroutine find_dl

  subroutine find_rp(t, rp)
    double precision :: t, tp
    double precision, dimension(3) :: rp
    integer :: vstart, vend

    if (t .lt. 0.25d0) then
       vstart = 1
       vend = 2
       tp = (t) / 0.25d0
    else if (t .lt. 0.50d0) then
       vstart = 2
       vend = 3
       tp = (t - 0.25d0) / 0.25d0
    else if (t .lt. 0.75d0) then
       vstart = 3
       vend = 4
       tp = (t - 0.50d0) / 0.25d0
    else
       vstart = 4
       vend = 1
       tp = (t - 0.75d0) / 0.25d0
    end if

    rp(1) = vertices(1, vstart) + tp*(vertices(1, vend) - vertices(1, vstart))
    rp(2) = vertices(2, vstart) + tp*(vertices(2, vend) - vertices(2, vstart))
    rp(3) = vertices(3, vstart) + tp*(vertices(3, vend) - vertices(3, vstart))
  end subroutine find_rp

  subroutine TestLoop()
    double precision :: t
    double precision, dimension(3) :: dl, rp
    integer :: N, it
    double precision :: dt

    ! This subroutine stores the shape and tangents to files.
    ! - the wire shape r'(t) is stored in the file fort.20
    ! - the tangent to the wire dl'(t) is stored in fort.21

    N = 10
    dt = 1.0d0 / N
    t = 0.0d0

    do it=1, N
       call find_dl(t, dl)
       call find_rp(t, rp)

       WRITE(20, *) t, rp(1), rp(2), rp(3)
       WRITE(21, *) t, dl(1), dl(2), dl(3)
       t = t + dt
    end do

  end subroutine TestLoop
end module SquareLoop

module CircleLoop
  ! implement a circular loop of current-carrying wire
  use Constants
  implicit none
  double precision, dimension(3, 4) :: vertices
  double precision, parameter :: radius=1.0d0
  double precision, parameter :: field_at_origin = mu0*I/(2.0d0 * radius)

  contains

  subroutine setup()

  end subroutine setup

  subroutine find_dl(t, dl)
    double precision :: t
    double precision, dimension(3) :: dl

    dl(1) = -2.0d0 * pi *radius*sin(2.0d0*pi*t)
    dl(2) = 2.0d0 * pi * radius*cos(2.0d0*pi*t)
    dl(3) = 0.0d0

  end subroutine find_dl

  subroutine find_rp(t, rp)
    double precision :: t
    double precision, dimension(3) :: rp

    rp(1) = radius*cos(2.0d0*pi*t)
    rp(2) = radius*sin(2.0d0*pi*t)
    rp(3) = 0.0d0

  end subroutine find_rp

  subroutine TestLoop()
    double precision :: t
    double precision, dimension(3) :: dl, rp
    integer :: N, it
    double precision :: dt

    ! This subroutine stores the shape and tangents to files.
    ! - the wire shape r'(t) is stored in the file fort.20
    ! - the tangent to the wire dl'(t) is stored in fort.21

    N = 50
    dt = 1.0d0 / N
    t = 0.0d0

    do it=1, N+1
       call find_dl(t, dl)
       call find_rp(t, rp)

       WRITE(20, *) t, rp(1), rp(2), rp(3)
       WRITE(21, *) t, dl(1), dl(2), dl(3)
       t = t + dt
    end do

  end subroutine TestLoop
end module CircleLoop

module StraightWire
  ! implement a long, straight current-carrying wire
  use Constants
  implicit none
  double precision, dimension(3, 2) :: vertices
  double precision, parameter :: distance = 0.5d0
  double precision, parameter :: length = 1.0d0
  double precision, parameter :: field_at_origin = mu0*I*length/(4.0d0*pi*distance * sqrt(distance**2 + (0.5d0*length)**2))
  ! NOTE: The wire itself is placed at y=distance, so that field_at_origin still provides an easy check

  contains

    subroutine setup()
      vertices(1,1) = 0.5d0 * length
      vertices(2,1) = distance
      vertices(3,1) = 0.0d0

      vertices(1,2) = -0.5d0 * length
      vertices(2,2) = distance
      vertices(3,2) = 0.0d0

    end subroutine setup

  subroutine find_dl(t, dl)
    double precision :: t
    double precision, dimension(3) :: dl
    integer, parameter :: vstart = 1, vend = 2

    ! calculate as dl(i) = \frac{d rp(i)}{dt}
    dl(1) = vertices(1, vend) - vertices(1, vstart)
    dl(2) = vertices(2, vend) - vertices(2, vstart)
    dl(3) = vertices(3, vend) - vertices(3, vstart)

  end subroutine find_dl

  subroutine find_rp(t, rp)
    double precision :: t
    double precision, dimension(3) :: rp
    integer, parameter :: vstart = 1, vend = 2

    rp(1) = vertices(1, vstart) + t*(vertices(1, vend) - vertices(1, vstart))
    rp(2) = vertices(2, vstart) + t*(vertices(2, vend) - vertices(2, vstart))
    rp(3) = vertices(3, vstart) + t*(vertices(3, vend) - vertices(3, vstart))

  end subroutine find_rp

  subroutine TestLoop()
    double precision :: t
    double precision, dimension(3) :: dl, rp
    integer :: N, it
    double precision :: dt

    ! This subroutine stores the shape and tangents to files.
    ! - the wire shape r'(t) is stored in the file fort.20
    ! - the tangent to the wire dl'(t) is stored in fort.21

    N = 50
    dt = 1.0d0 / N
    t = 0.0d0

    do it=1, N+1
       call find_dl(t, dl)
       call find_rp(t, rp)

       WRITE(20, *) t, rp(1), rp(2), rp(3)
       WRITE(21, *) t, dl(1), dl(2), dl(3)
       t = t + dt
    end do

  end subroutine TestLoop
end module StraightWire


module BiotSavartModule
  ! use the module that describes the grid points at which the magnetic field should be calculated
  use mesh
  ! use the module that describes the geometry of the desired current-carrying wire
  ! (comment the lines for the wire geometries you do not wish to work with presently)
  !use SquareLoop
  use StraightWire
  !use CircleLoop
  implicit none

  contains
  function IntegrandX(t) result(intx)
    double precision :: t, mag, intx, cross
    double precision, dimension(3) :: dl, rp

    call find_dl(t, dl)
    call find_rp(t, rp)

    cross = dl(2)*(r(3)-rp(3)) - dl(3)*(r(2) - rp(2))
    mag = sqrt(sum((r-rp)**2))

    if (mag .gt. 1.0d-5) then
       intx = cross / (mag**3)
    else
       intx = 0.0d0
    end if

  end function IntegrandX

  function IntegrandY(t) result(inty)
    double precision :: t, mag, inty, cross
    double precision, dimension(3) :: dl, rp

    call find_dl(t, dl)
    call find_rp(t, rp)

    cross = dl(3)*(r(1)-rp(1)) - dl(1)*(r(3) - rp(3))
    mag = sqrt(sum((r-rp)**2))

    if (mag .gt. 1.0d-5) then
       inty = cross / (mag**3)
    else
       inty = 0.0d0
    end if
  end function IntegrandY

  function IntegrandZ(t) result(intz)
    double precision :: t, mag, intz, cross
    double precision, dimension(3) :: dl, rp

    call find_dl(t, dl)
    call find_rp(t, rp)

    cross = dl(1)*(r(2)-rp(2)) - dl(2)*(r(1) - rp(1))
    mag = sqrt(sum((r-rp)**2))

    if (mag .gt. 1.0d-5) then
       intz = cross / (mag**3)
    else
       intz = 0.0d0
    end if
  end function IntegrandZ

  subroutine BiotSavartIntegral(Br)
    implicit none
    double precision, dimension(3) :: Br
    ! Things needed for Simpson Rule integration
    integer, parameter :: Nt = 200
    integer :: it
    double precision :: t, dt

!!$    WRITE(*,*) "DEBUG: Before calling integrals."

    ! perform each integral
    ! Use the Simpson Rule
    Br = 0.0d0
    dt = 1.0d0 / Nt
    t = 0.0d0
    Br(1) = IntegrandX(t)
    Br(2) = IntegrandY(t)
    Br(3) = IntegrandZ(t)
    do it=1, Nt/2-1
       ! odd terms
       t = dt * (2*it - 1)
       Br(1) = Br(1) + 4.0d0*IntegrandX(t)
       Br(2) = Br(2) + 4.0d0*IntegrandY(t)
       Br(3) = Br(3) + 4.0d0*IntegrandZ(t)
       ! even terms
       t = dt * (2*it)
       Br(1) = Br(1) + 2.0d0*IntegrandX(t)
       Br(2) = Br(2) + 2.0d0*IntegrandY(t)
       Br(3) = Br(3) + 2.0d0*IntegrandZ(t)
    end do
    ! pick up the last terms
    t = dt * (Nt - 1)
    Br(1) = Br(1) + 4.0d0*IntegrandX(t)
    Br(2) = Br(2) + 4.0d0*IntegrandY(t)
    Br(3) = Br(3) + 4.0d0*IntegrandZ(t)
    t = dt * (Nt)   ! that is, t = 1.0d0, the endpoint
    Br(1) = Br(1) + IntegrandX(t)
    Br(2) = Br(2) + IntegrandY(t)
    Br(3) = Br(3) + IntegrandZ(t)

    ! final factor
    Br = Br * mu0*dt / (3.0d0 * 4.0d0 * pi)

  end subroutine BiotSavartIntegral

end module BiotSavartModule


program BiotSavart
  use BiotSavartModule
  implicit none
  ! The magnetic field array
  !double precision, dimension(3, NX, NY, NZ) :: B
  double precision, dimension(3) :: B
  ! counters for the x, y, and z directions
  integer :: ix, iy, iz
  ! IO placeholder
  integer :: noutput = 60

  call setup()                  ! Constructs the current-carrying loop
  call TestLoop()               ! Prints out a sample of points in the loop

  ! a test at the center of the loop
  r(1)=0.0d0
  r(2)=0.0d0
  r(3)=0.0d0
  call BiotSavartIntegral(B)
  WRITE(*,*) "At center of loop, Bz=", B(3)
  WRITE(*,*) "Compare to ", field_at_origin
  WRITE(*,*) "    Percent error: ", 100.0d0*(abs(B(3))-field_at_origin)/field_at_origin

  ! open the file for output
  open(unit=noutput, file="MagField.txt", status="replace", form="formatted")

  do iz=1, NZ
     do iy=1, NY
        do ix=1, NX
           ! establish the vector r
           call calc_position(ix, iy, iz)
!!$           WRITE(*,*) "DEBUG: At r=", r(1), r(2), r(3)
           ! calculate B(:,ix,iy,iz) by the Biot-Savart integral
           call BiotSavartIntegral(B)
           write(noutput, *) r(1), r(2), r(3), B(1), B(2), B(3)
        end do ! ix
     end do ! iy
  end do ! iz

  ! close output file
  close(noutput)
  ! The data in the file "MagField.txt" may be plotted in gnuplot with the following commands:
  ! > enorm(x,y) = sqrt(x**2 + y**2)
  ! > plot "MagField.txt" using 2:3:(0.1*$5/enorm($5,$6)):(0.1*$6/enorm($5,$6)) w vectors

end program BiotSavart
