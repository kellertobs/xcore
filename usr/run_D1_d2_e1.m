% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  'D1_d2_e1';          % run identifier  (D = 1e2; d0 = 1e-2; etam0 = 1e1)
restart   =  0;                   % restart from file (0: new run; <0: restart from last; >0: restart from specified frame)
nop       =  50;                  % output frame plotted/saved every 'nop' time steps
nrh       =  5;                   % record metrics history every 'nrh' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
ndm_op    =  0;                   % plot nondimensionalised output 

% set model domain parameters
D         =  1e1;                 % chamber depth [m]
N         =  200;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D*1.5;               % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
t0end     =  1.5;                 % stop when dimensionless time is reached

% set crystallinity initial condition
xeq       =  0.01;                % equilibrium crystallinity of boundary layer [wt]
x0        =  xeq/10;              % initial background crystallinity [wt]

% set physical control parameters
d0        =  1e-2;                % xtal size constant [m]
etam0     =  1e+1;                % melt viscosity constant [kg/m3]
L0        =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
l0        =  d0*10;               % correlation length for phase fluctuation diffusivity (multiple of d0, 10-20)
R         =  0.5;                 % relative amplitude of crystallisation rate [s]
Xi        =  0.25;                % relative amplitude of random noise flux

% set numerical model parameters
CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-5;                % outer its relative tolerance
atol      =  1e-8;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
alpha     =  0.9;                 % iterative step size parameter
gamma     =  1e-3;                % artificial horizontal inertia parameter (only applies if periodic)


%*****  RUN XCORE MODEL  **************************************************
run('../src/main')
%**************************************************************************

