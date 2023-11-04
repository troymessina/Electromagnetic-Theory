! Written by K. Roos
! Department of Mechanical Engineering
! Bradley University
! 
! NOTE: This implementation of the Charges in a Cubic Conductor
! (originally produced by Larry Engelhardt in python) contains the basic
! structure for calculating the dynamical motion of charges that interact
! via the Coulomb interation, using a typical Molecular Dynamics pair
! potential approach.  The x,y, and z coordinates of all the charges are dumped
! to a data file "positions.dat" after the specificed number of time steps has
! been iterated. The components for the electric field due to all the single charges,
! as well as the magnitude of the electric field vector, are calculated at
! every time step, and dumped, along with the time, in a data file "field.dat"

program charge_cube
implicit none

! variable and array declarations
integer i,n_charges,n,t_steps,j
real(8) Q,m,L,k,qi,dt,r1,r2,dx,dy,dz,r_mag,f_mag,d
real(8) fij,fijx,fijy,fijz,Px,Py,Pz,r_vector
real(8) dx_P,dy_P,dz_P,E_field
real(8), allocatable :: x(:),y(:),z(:),vx(:),vy(:),vz(:),time(:)
real(8), allocatable :: fx(:),fy(:),fz(:)
real(8), allocatable :: Ex(:),Ey(:),Ez(:),E_mag(:)

! prepare data files for dumping
open(1,FILE='positions.dat')
open(2,FILE='field.dat')

CALL RANDOM_SEED

! set physical parameter values
n_charges=100 ! number of individual charges
Q=5.0e-6 ! net charge in Coulombs
m=1.0e-3 ! mass of each charge in kg
L=0.2 ! side length of conducting cube in meters
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
allocate(fx(n_charges),fy(n_charges),fz(n_charges))
allocate(time(t_steps),Ex(t_steps),Ey(t_steps),Ez(t_steps),E_mag(t_steps))

x=0; y=0; z=0; vx=0; vy=0; vz=0; fx=0; fy=0; fz=0
Ex=0; Ey=0; Ez=0; E_mag=0

! create charges at random initial positions (all initially at rest)
! The center of the cube is at the origin.
do i=1,n_charges
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    x(i)=(-1)**int(2.0*r1)*0.5*L*r2
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    y(i)=(-1)**int(2.0*r1)*0.5*L*r2
    CALL RANDOM_NUMBER(r1);CALL RANDOM_NUMBER(r2)
    z(i)=(-1)**int(2.0*r1)*0.5*L*r2
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
      ! For the algorithm for updating positions, use the large friction
      ! limit wherein the force is replaced by displacement. This approach 
      ! is the simplest stable algorithm for this geometry; thus, this
      ! simulation emphasizes the resultant configuration of excess charge
      ! after equilibration, rather than the correct physics governing the
      ! dynamical motion of the system.
      do i=1,n_charges
         f_mag=sqrt(fx(i)*fx(i)+fy(i)*fy(i)+fz(i)*fz(i))
         ! With force representing displacement, rescale displacement
         ! in the case of a large calculated displacement
         if(f_mag>L/100) then
            Fx(i)=0.01*L*Fx(i)/f_mag
            Fy(i)=0.01*L*Fy(i)/f_mag
            Fz(i)=0.01*L*Fz(i)/f_mag  
         end if
         ! Update positions with large friction approximation
         x(i)=x(i)+fx(i)
         y(i)=y(i)+fy(i);
         z(i)=z(i)+fz(i);
        
        ! if charge moved beyond surface of conductor, bring it back to surface
        if (x(i)>0.5*L) then 
            x(i)=0.5*L
        end if
        if (x(i)<-0.5*L) then 
            x(i)=-0.5*L
        end if
        if (y(i)>0.5*L) then
            y(i)=0.5*L
        end if
        if (y(i)<-0.5*L) then
            y(i)=-0.5*L
        end if
        if (z(i)>0.5*L) then
            z(i)=0.5*L
        end if
        if (z(i)<-0.5*L) then
            z(i)=-0.5*L
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
if (Px>0.5*L.or.Px<-0.5*L.or.Py>0.5*L.or.Py<-0.5*L.or.Pz>0.5*L.or.Pz<-0.5*L) then
    ! The field for all charge concentrated at the center of the cube is
    ! calculated for comparison to the graph. Due to cubic charge distribution
    ! E_mag on the graph should not match the value calculated here.
    print*, 'Field point is outside cube'
    E_field=1e-6*k*Q/(Px**2+Py**2+Pz**2)
    write(*,*) 'For total charge concentrated at the center, E_field=', E_field;
else
    print*, 'Field point is inside cube'
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
         
end program charge_cube