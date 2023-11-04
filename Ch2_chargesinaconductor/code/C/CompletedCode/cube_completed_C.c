/* 
Written by K. Roos
Department of Mechanical Engineering
Bradley University
 
NOTE: This implementation of the Charges in a Cubic Conductor
(originally produced by Larry Engelhardt in python) contains the basic
structure for calculating the dynamical motion of charges that interact
via the Coulomb interation, using a typical Molecular Dynamics pair
potential approach.  The x,y, and z coordinates of all the charges are dumped
to a data file "positions.dat" after the specified number of time steps has
been iterated. The components for the electric field due to all the single charges,
as well as the magnitude of the electric field vector, are calculated at
every time step, and dumped, along with the time, in a data file "field.dat"
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

int main () {
	// variable and array declarations 
  int i,n_charges,n,t_steps,j;
  double Q,m,L,k,qi,dt,r1,r2,dx,dy,dz,r_mag,f_mag,d;
  double fij,fijx,fijy,fijz,Px,Py,Pz,r_vector;
  double dx_P,dy_P,dz_P,E_field;
  double *x,*y,*z,*vx,*vy,*vz,*t,*fx,*fy,*fz;
  double *Ex,*Ey,*Ez,*E_mag;
  FILE *positions_file,*field_file;

  srand(time(NULL)); // initialize the random number generator

  // set physical parameter values 
  n_charges=100; // number of individual charges
  Q=5.0e-6; // net charge in Coulombs
  m=1.0e-3; // mass of each charge in kg
  L=0.2; // length of cube side in meters
  k=8.99e9; // Coulomb constant
  qi=Q/(double)n_charges; // charge on each individual charge
  
  // algorithmic parameters
  dt=0.001; // time step in seconds
  t_steps=100; // total number of time steps

  /* Coordinates of field point P -- Electric field components due 
  to all the charges will be calculated at this point */
  Px=-0.05;
  Py=0.043;
  Pz=-0.01;

  /* Allocate memory and initialize contents to zero */
  x=(double*)calloc(n_charges, sizeof(double));
  y=(double*)calloc(n_charges, sizeof(double));
  z=(double*)calloc(n_charges, sizeof(double));
  vx=(double*)calloc(n_charges, sizeof(double));
  vy=(double*)calloc(n_charges, sizeof(double));
  vz=(double*)calloc(n_charges, sizeof(double));
  fx=(double*)calloc(n_charges, sizeof(double));
  fy=(double*)calloc(n_charges, sizeof(double));
  fz=(double*)calloc(n_charges, sizeof(double));
  t=(double*)calloc(t_steps, sizeof(double)); // time array
  Ex=(double*)calloc(t_steps, sizeof(double));
  Ey=(double*)calloc(t_steps, sizeof(double));
  Ez=(double*)calloc(t_steps, sizeof(double));
  E_mag=(double*)calloc(t_steps, sizeof(double));

  /* Open the files used to store the results */
  positions_file=fopen("positions.dat","w");
  field_file=fopen("field.dat","w");

	/* create charges at random initial positions (all initially at rest)
    The center of the cube is at the origin.*/
	for(i=1; i<=n_charges; i++) {
		r1=int(2.0*(double)rand()/RAND_MAX);
		r2=(double)rand()/RAND_MAX;
		x[i]=pow(-1,r1)*0.5*L*r2;
		r1=int(2.0*(double)rand()/RAND_MAX);
		r2=(double)rand()/RAND_MAX;
		y[i]=pow(-1,r1)*0.5*L*r2;
		r1=int(2.0*(double)rand()/RAND_MAX);
		r2=(double)rand()/RAND_MAX;
		z[i]=pow(-1,r1)*0.5*L*r2;
	}

	// Main loop
	for(n=1; n<=t_steps; n++) {
		t[n]=(double)n*dt;
		// reset net force in x,y,z directions to zero for all charges
		for(i=1; i<=n_charges; i++) {
			fx[i]=0;
			fy[i]=0;
			fz[i]=0;
		}
		/* Calculate forces on each pair of charges by considering i-j pairs.
        The structure of the nested for loops prevents double counting of
        pairs.*/
		for(i=1; i<=n_charges-1; i++) {
			for(j=i+1; j<=n_charges; j++) {
				dx=x[i]-x[j];
				dy=y[i]-y[j];
				dz=z[i]-z[j];
				r_mag=sqrt(dx*dx+dy*dy+dz*dz);
				// Magnitude of mutual Coulomb force between charges i and j
                fij=k*qi*qi/(r_mag*r_mag);
				// The three spatial components of the force between charges i and j
                fijx=fij*dx/r_mag;
                fijy=fij*dy/r_mag;
                fijz=fij*dz/r_mag;
			    // Accumulation of the componenets of the net force acting on charges i and j
			    fx[i]=fx[i]+fijx;
                fy[i]=fy[i]+fijy;
                fz[i]=fz[i]+fijz;
                fx[j]=fx[j]-fijx;
                fy[j]=fy[j]-fijy;
                fz[j]=fz[j]-fijz;
			}
		}
		// velocities and positions updated for all charges
		/*  For the algorithm for updating positions, use the large friction
		limit wherein the force is replaced by displacement. This approach 
		is the simplest stable algorithm for this geometry; thus, this
		simulation emphasizes the resultant configuration of excess charge
		after equilibration, rather than the correct physics governing the
		dynamical motion of the system.
		*/
		for(i=1; i<=n_charges-1; i++) {
			f_mag=sqrt(fx[i]*fx[i]+fy[i]*fy[i]+fz[i]*fz[i]);  
            /* With force representing displacement, rescale displacement
            in the case of a large calculated displacement*/
			if(f_mag>L/1000) {
				fx[i]=0.01*L*fx[i]/f_mag;
                fy[i]=0.01*L*fy[i]/f_mag;
                fz[i]=0.01*L*fz[i]/f_mag;
			}
			// Update positions with large friction approximation
            x[i]=x[i]+fx[i];
            y[i]=y[i]+fy[i];
            z[i]=z[i]+fz[i];
			
			// if charge moved beyond surface of conuctor, bring it back to surface

			if (x[i]>0.5*L) {
				x[i]=0.5*L;
			}
			if (x[i]<-0.5*L) {
				x[i]=-0.5*L;
			}
			if (y[i]>0.5*L) {
				y[i]=0.5*L;
			}
			if (y[i]<-0.5*L) {
				y[i]=-0.5*L;
			}
			if (z[i]>0.5*L) {
				z[i]=0.5*L;
			}
			if (z[i]<-0.5*L) {
				z[i]=-0.5*L;
			}
			
			// Calculation of Electric Field components at point P
            dx_P=Px-x[i];
            dy_P=Py-y[i];
            dz_P=Pz-z[i];
            r_vector=sqrt(dx_P*dx_P+dy_P*dy_P+dz_P*dz_P);
			/* Each E-field component, as well as the magnitude of the
            electric field vector is calculated. Units of micro-N/C for
            the electric field are used for ease in plotting.*/
			Ex[n]=Ex[n]+1.0e-6*k*qi/pow(r_vector,2)*(dx_P/r_vector);
            Ey[n]=Ey[n]+1.0e-6*k*qi/pow(r_vector,2)*(dy_P/r_vector);
            Ez[n]=Ez[n]+1.0e-6*k*qi/pow(r_vector,2)*(dz_P/r_vector);
            E_mag[n]=sqrt(Ex[n]*Ex[n]+Ey[n]*Ey[n]+Ez[n]*Ez[n]);
		}
	}
	if(Px>0.5*L || Px<-0.5*L || Py>0.5*L || Py<-0.5*L || Pz>0.5*L || Pz<-0.5*L) {
        /* The field for all charge concentrated at the center of the cube is
		calculated for comparison to the graph. Due to cubic charge distribution
		E_mag on the graph should not match the value calculated here.*/
		cout<<"Field point is outside cube. ";
		E_field=1e-6*k*Q/(Px*Px+Py*Py+Pz*Pz);
		cout<<"For total charge concentrated at the center, E_field="<<E_field<<endl;
    } else {
        cout<<"Field point is inside cube";
    }
    // Dump coordinates of charges into "positions.dat"
	for(i=1; i<=n_charges; i++) {
		 fprintf(positions_file,"%3.3f %3.3f %3.3f\n",x[i],y[i],z[i]);
	}
	// Dump time and  E-field components to data file "field.dat"
	for(i=1; i<=t_steps; i++) {
		fprintf(field_file,"%3.3f %3.3f %3.3f %3.3f %3.3f\n",t[i],Ex[i],Ey[i],Ez[i],E_mag[i]);
	}
  
  /* Close the files */
	fclose(positions_file);
	fclose(field_file);

  return 0;
}

