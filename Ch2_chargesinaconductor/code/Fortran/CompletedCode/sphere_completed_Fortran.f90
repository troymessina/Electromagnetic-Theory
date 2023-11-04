! Written by K. Roos
! Department of Mechanical Engineering
! Bradley University
! 
! NOTE: This implementation of the Charges in a Spherical Conductor
! (originally produced by Larry Engelhardt in python) contains the basic
! structure for calculating the dynamical motion of charges that interact
! via the Coulomb interation, using a typical Molecular Dynamics pair
! potential approach.  The x,y, and z coordinates of all the charges are dumped
! to a data file "positions.dat" after the specified number of time steps has
! been iterated. The components for the electric field due to all the single charges,
! as well as the magnitude of the electric field vector, are calculated at
! every time step, and dumped, along with the time, in a data file "field.dat"

program charge_sphere
implicit none

! variable and array declarations
integer i,n_charges,n,t_steps,j
real(8) Q,m,R,k,qi,dt,r1,r2,dx,dy,dz,r_mag,a_mag,d
real(8) fij,fijx,fijy,fijz,Px,Py,Pz,r_vector
real(8) dx_P,dy_P,dz_P,E_field
real(8), allocatable :: x(:),y(:),z(:),vx(:),vy(:),vz(:),time(:)
real(8), allocatable :: ax(:),ay(:),az(:),fx(:),fy(:),fz(:)
real(8), allocatable :: Ex(:),Ey(:),Ez(:),E_mag(:)

! prepare data files for dumping
open(1,FILE='positions.dat')
open(2,FILE='field.dat')

CALL RANDOM_SEED

! set physical parameter values
n_charges=100 ! number of individual charges
Q=5.0e-6 ! net charge in Coulombs
m=1.0e-3 ! mass of each charge in kg
R=0.1 ! radius of conducting sphere in meters
k=8.99e9 ! Coulomb constant
qi=Q/real(n_charges) ! charge on each individual charge

! algorithmic parameters
dt=0.001 ! time step in seconds
t_steps=100 ! total number of time steps

! Coordinates of field point P -- Electric field components due 
! to all the charges will be calculated at this point
Px=-0.05
Py=0.043
Pz=-0.01

! Preallocate arrays to be used for calculations
allocate(x(n_charges),y(n_charges),z(n_charges),vx(n_charges),vy(n_charges),vz(n_charges))
allocate(ax(n_charges),ay(n_charges),az(n_charges),fx(n_charges),fy(n_charges),fz(n_charges))
allocate(time(t_steps),Ex(t_steps),Ey(t_steps),Ez(t_steps),E_mag(t_steps))

x=0; y=0; z=0; vx=0; vy=0; vz=0; ax=0; ay=0; az=0; fx=0; fy=0; fz=0
Ex=0; Ey=0; Ez=0; E_mag=0

! create charges at random initial positions (all initially at rest)
! The center of the spherical conductor of radius R is at the origin.
do i=1,n_charges
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    x(i)=(-1)**int(2.0*r1)*R/sqrt(3.0)*r2
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    y(i)=(-1)**int(2.0*r1)*R/sqrt(3.0)*r2
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    z(i)=(-1)**int(2.0*r1)*R/sqrt(3.0)*r2
end do        

! Main loop
do n=1,t_steps 
     time(n)=real(n)*dt
!    ! reset net force in x,y,z directions to zero for all charges   
     fx=0; fy=0; fz=0
     ! Calculate forces on each pair of charges by considering i-j pairs.
     ! The structure of the nested for loops prevents double counting of
     ! pairs.
     do i=1,n_charges-1
        do j=i+1,n_charges
            dx=x(i)-x(j)
            dy=y(i)-y(j)
            dz=z(i)-z(j)
            r_mag=sqrt(dx*dx+dy*dy+dz*dz)
            ! Magnitude of mutual Coulomb force between charges i and j
            fij=k*qi*qi/(r_mag*r_mag)
            ! The three spatial components of the force between charges i and j
            fijx=fij*dx/r_mag
            fijy=fij*dy/r_mag
            fijz=fij*dz/r_mag
            ! Accumulation of the componenets of the net force acting on charges i and j
            fx(i)=fx(i)+fijx
            fy(i)=fy(i)+fijy
            fz(i)=fz(i)+fijz
            fx(j)=fx(j)-fijx
            fy(j)=fy(j)-fijy
            fz(j)=fz(j)-fijz
         end do
      end do
      ! velocities and positions updated for all charges  
      do i=1,n_charges
         ax(i)=fx(i)/m
         ay(i)=fy(i)/m
         az(i)=fz(i)/m
         a_mag=sqrt(ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i))  
         ! rescale acceleration in the case of a huge calculated acceleration
         if(a_mag>1000)then
            ax(i)=1000.0*ax(i)/a_mag
            ay(i)=1000.0*ay(i)/a_mag
            az(i)=1000.0*az(i)/a_mag
         end if
         ! Euler-Cromer algorithm for updating velocities and positions
         vx(i)=vx(i)+ax(i)*dt
         vy(i)=vy(i)+ay(i)*dt
         vz(i)=vz(i)+az(i)*dt
         x(i)=x(i)+vx(i)*dt
         y(i)=y(i)+vy(i)*dt
         z(i)=z(i)+vz(i)*dt
         ! d is the distance from center of sphere to charge
         d=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         ! if charge moved beyond surface of conuctor, bring it back to surface
         if(d>R)then
            x(i)=x(i)*R/d
            y(i)=y(i)*R/d
            z(i)=z(i)*R/d
         end if
        ! Calculation of Electric Field components at point P
        dx_P=Px-x(i)
        dy_P=Py-y(i)
        dz_P=Pz-z(i)
        r_vector=sqrt(dx_P**2+dy_P**2+dz_P**2)
        ! Each E-field component, as well as the magnitude of the
        ! electric field vector is calculated. Units of micro-N/C for
        ! the electric field are used for ease in plotting.
        Ex(n)=Ex(n)+1.0e-6*k*qi/r_vector**2*(dx_P/r_vector) 
        Ey(n)=Ey(n)+1.0e-6*k*qi/r_vector**2*(dy_P/r_vector)
        Ez(n)=Ez(n)+1.0e-6*k*qi/r_vector**2*(dz_P/r_vector)
        E_mag(n)=sqrt(Ex(n)**2+Ey(n)**2+Ez(n)**2)
      end do
end do ! end of Main loop
if (sqrt(Px**2+Py**2+Pz**2)>R) then
    ! If the field point P is outside the sphere, the electric field at P
    ! is equivalent to the field produced if the total charge were
    ! concentrated at the center of the sphere.  This value is calculated
    ! here for comparison with E_mag in the graph.
    print*, 'Field point is outside sphere'
    E_field=1e-6*k*Q/(Px**2+Py**2+Pz**2)
    write(*,*) 'For total charge concentrated at the center, E_field=', E_field;
else
    print*, 'Field point is inside sphere'
end if
! Dump coordinates of charges into "positions.dat"
do i=1,n_charges
    write(1,*) x(i),y(i),z(i)
end do 
! Dump time and  E-field components to data file "field.dat"
do i=1,t_steps
    write(2,"(5F7.3)") time(i),Ex(i),Ey(i),Ez(i),E_mag(i)
end do

close(1)
close(2)
         
end program charge_sphere