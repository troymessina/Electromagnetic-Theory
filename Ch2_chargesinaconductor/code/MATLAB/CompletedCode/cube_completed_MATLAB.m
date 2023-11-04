% Written by K. Roos
% Department of Mechanical Engineering
% Bradley University
% 
% NOTE: This implementation of the Charges in a Cubic Conductor
% (originally produced by Larry Engelhardt in python) contains the basic
% structure for calculating the dynamical motion of charges that interact
% via the Coulomb interation, using a typical Molecular Dynamics pair
% potential approach.  It does not attempt to animate the motion of the
% charges during the calculation (although there are several ways to
% accomplish this in MATLAB), but only plots the final position of the
% charges (including an image of the cubic conductor) 
% in a 3D scatter plot after the program has looped through the
% specified number of time steps.

clear; close all; clc;

% set physical parameter values
n_charges=100; % number of individual charges
Q=5e-6; % net charge in Coulombs
m=1e-3; % mass of each charge in kg
L=0.2; % side length of conducting cube in meters
k=8.99e9; % Coulomb constant
qi=Q/n_charges; % charge on each individual charge

% Coordinates of field point P -- Electric field components due 
% to all the charges will be calculated at this point
Px=0.07;
Py=-0.045;
Pz=0.09;

% algorithmic parameters
dt=0.001; % time step in seconds
t_steps=100; % total number of time steps

% Preallocate arrays to be used for calculating electric Field components
% and time
Ex=zeros(1,t_steps);
Ey=zeros(1,t_steps);
Ez=zeros(1,t_steps);
E_mag=zeros(1,t_steps);
time=zeros(1,t_steps);

% Preallocate arrays for position and velocity of charges
% as well as the net force on each charge
x_pos=zeros(1,n_charges);
y_pos=zeros(1,n_charges);
z_pos=zeros(1,n_charges);

vx=zeros(1,n_charges);
vy=zeros(1,n_charges);
vz=zeros(1,n_charges);
v=zeros(1,n_charges);

Fx=zeros(1,n_charges);
Fy=zeros(1,n_charges);
Fz=zeros(1,n_charges);

% Create charges with random initial positions (all initially at rest).
% The center of the cubic conductor of length L is at the origin.
for i=1:n_charges
        x_pos(i)=(-1)^floor(2*rand)*L/2*rand;
        y_pos(i)=(-1)^floor(2*rand)*L/2*rand;
        z_pos(i)=(-1)^floor(2*rand)*L/2*rand;    
end

% The next 7 commands create a figure with the image of the surface of the
% cubic conductor plotted with the aid of the function "plotcube" defined
% below.  The "hold on" command prepares the
% interpreter to include the scatter plot of the charges after the main
% loop over time steps is finished.
cube_origin = [0,0,0] ;   % center point of cube
cube_length = [L,L,L] ;  % your cube dimensions 
O = cube_origin-cube_length/2 ;       % Get the origin of cube so that P is at center
figure('Position', [400, 100, 550, 600]);
plotcube(cube_length,O,0.5,[0.4 0.5 0.5]);   % use function plotcube
axis([-L/2 L/2 -L/2 L/2 -L/2 L/2])
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
            % Accumulation of the components of the net force acting on
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
    % For the algorithm for updating positions, use the large friction
    % limit wherein the force is replaced by displacement. This approach 
    % is the simplest stable algorithm for this geometry; thus, this
    % simulation emphasizes the resultant configuration of excess charge
    % after equilibration, rather than the correct physics governing the
    % dynamical motion of the system.
    for i=1:n_charges
        f_mag=sqrt(Fx(i)*Fx(i)+Fy(i)*Fy(i)+Fz(i)*Fz(i));
        % With force representing displacement, rescale displacement
        % in the case of a large calculated displacement
        if(f_mag>L/100)
            Fx(i)=(L/100)*Fx(i)/f_mag;
            Fy(i)=(L/100)*Fy(i)/f_mag;
            Fz(i)=(L/100)*Fz(i)/f_mag;  
        end
        % Update positions with large friction approximation
        x_pos(i)=x_pos(i)+Fx(i);
        y_pos(i)=y_pos(i)+Fy(i);
        z_pos(i)=z_pos(i)+Fz(i);
        
        % if charge moved beyond surface of conductor, bring it back to surface
        if x_pos(i)>L/2 
            x_pos(i)=L/2;
        end
        if x_pos(i)<-L/2 
            x_pos(i)=-L/2;
        end
        if y_pos(i)>L/2 
            y_pos(i)=L/2;
        end
        if y_pos(i)<-L/2 
            y_pos(i)=-L/2;
        end
        if z_pos(i)>L/2 
            z_pos(i)=L/2;
        end
        if z_pos(i)<-L/2 
            z_pos(i)=-L/2;
        end
        
        % Calculation of Electric Field components at point P
        dx_P=Px-x_pos(i);
        dy_P=Py-y_pos(i);
        dz_P=Pz-z_pos(i);
        r_vector=sqrt(dx_P^2+dy_P^2+dz_P^2);
        % Each E-field component, as well as the magnitude of the
        % electric field vector is calculated. Units of micro-N/C for
        % the electric field are used for ease in plotting.
        Ex(n)=Ex(n)+1e-6*k*qi/r_vector^2*(dx_P/r_vector); 
        Ey(n)=Ey(n)+1e-6*k*qi/r_vector^2*(dy_P/r_vector);
        Ez(n)=Ez(n)+1e-6*k*qi/r_vector^2*(dz_P/r_vector);
        E_mag(n)=sqrt(Ex(n)^2+Ey(n)^2+Ez(n)^2);
    end
end % end of Main loop
if Px>L/2 || Px<-L/2 || Py>L/2 || Py<-L/2 || Pz>L/2 || Pz<-L/2
    disp('Field point is outside cube')
    % The field for all charge concentrated at the center of the cube is
    % calculated for comparison to the graph. Due to cubic charge distribution
    % E_mag on the graph should not match the value calculated here.
    E_field=1e-6*k*Q/(Px^2+Py^2+Pz^2);
    fprintf('For total charge concentrated at the center, E= %f \n', E_field);
else
    disp('Field point is inside cube')
end
% Plotting commands
scatter3(x_pos,y_pos,z_pos,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
title('Charges in a Cubic Conductor');
figure;
plot(time,Ex,time,Ey,time,Ez,time,E_mag);
title('Electric Field components vs. Time');
xlabel('Time');
ylabel('Field Component or Magnitude (\mu N/C');
legend('Ex','Ey','Ez','E_{mag}');

% function that produces cube image
function plotcube(varargin)
inArgs(1:nargin) = varargin;
% Create variables
[edges,origin,alpha,clr] = deal(inArgs{:});
XYZ = { ...
  [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
  [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
  [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
  [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
  [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
  [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
  };
XYZ = mat2cell(...
  cellfun( @(x,y,z) x*y+z , ...
    XYZ , ...
    repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
    repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
    'UniformOutput',false), ...
  6,[1 1 1]);
cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
  repmat({clr},6,1),...
  repmat({'FaceAlpha'},6,1),...
  repmat({alpha},6,1)...
  );
view(3);
end