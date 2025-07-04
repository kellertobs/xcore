 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  'testRNS_D02_d2_e1'; % run identifier  (D = 1e2; d0 = 1e-2; etam0 = 1e1)
restart   =  0;                   % restart from file (0: new run; <0: restart from last; >0: restart from specified frame)
nop       =  20;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file

% set model domain parameters
D         =  2;                   % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D/2.0;               % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
dt        =  0.01;                % initial time step [s]

% set physical parameters
xeq       =  0.01;                % equilibrium crystallinity of boundary layer [wt]
x0        =  xeq/1;              % initial background crystallinity [wt]
dx0       =  x0/10;               % background crystallinity perturbation [wt]
d0        =  1e-2;                % xtal size constant [m]
etam0     =  1e+1;                % melt viscosity constant [kg/m3]

% set boundary layer parameters
Da        =  0.1;                 % Dahmköhler number of boundary layer rate [s]
closed_bot=  0;                   % switch for closed bottom boundary to form cumulate pile

% set numerical model parameters
CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-5;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
alpha     =  0.9;                 % iterative step size parameter
Delta_cnv =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
Delta_sgr =  d0*10;               % correlation length for phase fluctuation diffusivity (multiple of d0, 10-20)
xi        =  1.0;                 % relative amplitude of random noise flux
gamma     =  1e-3;                % artificial horizontal inertia parameter (only applies if periodic)
bnd_w     =  max(Delta_sgr/2,2*h);    % width of boundary layer [m]

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

