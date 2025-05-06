 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  'D100_d1_e2';        % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file

% set model domain parameters
D         =  100;                 % chamber depth [m]
N         =  200;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
dt        =  0.1;                 % initial time step [s]

% set physical parameters
x0        =  0.01;                % initial background crystallinity [wt]
dx0       =  x0/3;                % background crystallinity perturbation [wt]
xb        =  0.00;                % initial boundary layer crystallinity [wt]
dxb       =  xb/3;                % boundary layer crystallinity perturbation [wt]

d0        =  0.01;                % xtal size constant [kg/m3]

% set boundary layer parameters
bnd_w     =  2*h;                 % width of boundary layer [m]
xeq       =  0.1;                 % equilibrium crystallinity of boundary layer [wt]
Da        =  0.1;                 % Dahmköhler number of boundary layer rate [s]

% set numerical model parameters
CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-5;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
alpha     =  0.75;                % iterative step size parameter
Delta_cnv =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
Delta_sgr =  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
Rec       =  1;                   % critical Reynolds number for ramping up eddy diffusivity
Scx       =  1;                   % xtal Schmidt number for applying eddy diffusivity to xtal diffusivity


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

