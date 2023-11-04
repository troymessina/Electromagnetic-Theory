% Written by K. Roos
% Department of Mechanical Engineering
% Bradley University
% 
% NOTE: This implementation of the Charges in a Spherical Conductor
% (originally produced by Larry Engelhardt in python) contains the basic
% structure for calculating the dynamical motion of charges that interact
% via the Coulomb interation, using a typical Molecular Dynamics pair
% potential approach.  It does not attempt to animate the motion of the
% charges during the calculation (although there are several ways to
% accomplish this in MATLAB), but only plots the final position of the
% charges (including an image of the spherical conductor) 
% in a 3D scatter plot after the program has looped through the
% specified number of time steps.

clear; close all; clc;

% set physical parameter values
n_charges=100; % number of individual charges
Q=5e-6; % net charge in Coulombs
m=1e-3; % mass of each charge in kg
R=0.1; % radius of conducting sphere in meters
k=8.99e9; % Coulomb constant
qi=Q/n_charges; % charge on each individual charge

% Coordinates of field point P -- Electric field components due 
% to all the charges will be calculated at this point
Px=-0.05;
Py=0.043;
Pz=-0.01;

% algorithmic parameters
dt=0.001; % time step in seconds
t_steps=1000; % total number of time steps

% Preallocate arrays to be used for calculating electric Field components
% and time
Ex=zeros(1,t_steps);
Ey=zeros(1,t_steps);
Ez=zeros(1,t_steps);
E_mag=zeros(1,t_steps);
time=zeros(1,t_steps);

% Preallocate arrays for position, velocity, and acceleration of charges
% as well as the net force on each charge
x_pos=zeros(1,n_charges);
y_pos=zeros(1,n_charges);
z_pos=zeros(1,n_charges);

vx=zeros(1,n_charges);
vy=zeros(1,n_charges);
vz=zeros(1,n_charges);
v=zeros(1,n_charges);

ax=zeros(1,n_charges);
ay=zeros(1,n_charges);
az=zeros(1,n_charges);

Fx=zeros(1,n_charges);
Fy=zeros(1,n_charges);
Fz=zeros(1,n_charges);

% Create charges with random initial positions (all initially at rest).
% The center of the spherical conductor of radius R is at the origin.
for i=1:n_charges
        x_pos(i)=(-1)^floor(2*rand)*R/sqrt(3)*rand;
        y_pos(i)=(-1)^floor(2*rand)*R/sqrt(3)*rand;
        z_pos(i)=(-1)^floor(2*rand)*R/sqrt(3)*rand;    
end

% The next 7 commands create a figure with the image of the surface of the
% spherical conductor plotted.  The "hold on" command prepares the
% interpreter to include the scatter plot of the charges after the main
% loop over time steps is finished.
[x,y,z] = sphere(100);
figure('Position', [400, 100, 550, 600]);
h=surf(R*x,R*y,R*z);
set(h,'edgecolor','none','facecolor',[0.4 0.5 0.5]);
axis([-R R -R R -R R])
alpha 0.5;
hold on;

% Main loop
for n=1:t_steps % n is time index
    time(n)=n*dt; 
    % reset net force in x,y,z directions to zero for all particles
    for p=1:n_charges
        Fx(p)=0;
        Fy(p)=0;
        Fz(p)=0;
    end
    % Calculate forces on each pair of charges by considering i-j pairs.
    % The structure of the nested for loops prevents double counting of
    % pairs.
    for i=1:(n_charges-1)
        for j=(i+1):n_charges
            dx=x_pos(i)-x_pos(j);
            dy=y_pos(i)-y_pos(j);
            dz=z_pos(i)-z_pos(j);
            r_mag=sqrt(dx*dx+dy*dy+dz*dz);
            % Magnitude of mutual Coulomb force between charges i and j
            fij=k*qi*qi/(r_mag*r_mag); 
            % The three spatial components of the force between charges i
            % and j
            fijx=fij*dx/r_mag;
            fijy=fij*dy/r_mag;
            fijz=fij*dz/r_mag;
            % Accumulation of the componenets of the net force acting on
            % charges i and j
            Fx(i)=Fx(i)+fijx;
            Fy(i)=Fy(i)+fijy;
            Fz(i)=Fz(i)+fijz;
            Fx(j)=Fx(j)-fijx;
            Fy(j)=Fy(j)-fijy;
            Fz(j)=Fz(j)-fijz; 
        end
    end
    % velocities and positions updated for all charges
    for i=1:n_charges
        ax(i)=Fx(i)/m;
        ay(i)=Fy(i)/m;
        az(i)=Fz(i)/m;
        a_mag=sqrt(ax(i)*ax(i)+ay(i)*ay(i)+az(i)*az(i));
        % rescale acceleration in the case of a huge calculated
        % acceleration
        if(a_mag>1000)
            ax(i)=1000*ax(i)/a_mag;
            ay(i)=1000*ay(i)/a_mag;
            az(i)=1000*az(i)/a_mag;  
        end
        % Euler-Cromer algorithm for updating velocities and positions
        vx(i)=vx(i)+ax(i)*dt;
        vy(i)=vy(i)+ay(i)*dt;
        vz(i)=vz(i)+az(i)*dt;
        x_pos(i)=x_pos(i)+vx(i)*dt;
        y_pos(i)=y_pos(i)+vy(i)*dt;
        z_pos(i)=z_pos(i)+vz(i)*dt;
        
        % d is distance from center of sphere to charge i
        d=sqrt(x_pos(i)*x_pos(i)+y_pos(i)*y_pos(i)+z_pos(i)*z_pos(i));
        % if charge moved beyond surface of conuctor, bring it back to surface
        if(d>R) 
            x_pos(i)=x_pos(i)*R/d;
            y_pos(i)=y_pos(i)*R/d;
            z_pos(i)=z_pos(i)*R/d;
        end
        % Calculation of Electric Field components at point P
        %----------------------------------------------------------------
        %
        % THE CALCULATION OF THE ELECTRIC FIELD COMPONENTS AND MAGNITUDE 
        % GOES HERE
        %
        %----------------------------------------------------------------
    end
end % end of Main loop
if sqrt(Px^2+Py^2+Pz^2)>R
    % If the field point P is outside the sphere, the electric field at P
    % is equivalent to the field produced if the total charge were
    % concentrated at the center of the sphere.  This value is calculate
    % here for comparison wit E_mag in the graph.
    disp('Field point is outside sphere')
    E_field=1e-6*k*Q/(Px^2+Py^2+Pz^2);
    fprintf('For total charge concentrated at the center, E= %f \n', E_field);
else
    disp('Field point is inside sphere')
end
scatter3(x_pos,y_pos,z_pos,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
title('Charges in a Spherical Conductor');
figure;
plot(time,Ex,time,Ey,time,Ez,time,E_mag);
title('Electric Field components vs. Time');
xlabel('Time');
ylabel('Field Component or Magnitude (\mu N/C');
legend('Ex','Ey','Ez','E_{mag}');

