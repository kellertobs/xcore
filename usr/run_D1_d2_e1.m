% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  'D1_d2_e1';          % run identifier  (D = 1e2; d0 = 1e-2; etam0 = 1e1)
restart   =  0;                   % restart from file (0: new run; <0: restart from last; >0: restart from specified frame)
nop       =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file

% set model domain parameters
D         =  1e2;                 % chamber depth [m]
N         =  200;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D*1.5;               % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
t0end     =  1e+3;                % stop when dimensionless time is reached
xend      =  0.10;                % stop run when mean crystallinity reaches threshold
Dxend     =  1e-6;                % stop run when time-averaged change of crystallinity drops below threshold
DVend     =  1e-9;                % stop run when time-averaged change of convective speed drops below threshold

% set physical parameters
xeq       =  0.01;                % equilibrium crystallinity of boundary layer [wt]
x0        =  xeq/10;              % initial background crystallinity [wt]
dx0       =  x0/10;               % background crystallinity perturbation [wt]
xb        =  xeq;                 % initial boundary layer crystallinity [wt]
dxb       =  xb/10;               % boundary layer crystallinity perturbation [wt]

d0        =  1e-2;                % xtal size constant [m]
etam0     =  1e+1;                % melt viscosity constant [kg/m3]
R         =  1.0;                 % relative amplitude of crystallisation rate [s]
Xi        =  1.0;                 % relative amplitude of random noise flux
closed    =  1;                   % switch for closed bottom velocity boundary to form cumulate pile
open_sgr  =  1;                   % switch for open bottom boundary for crystal segregation

% set numerical model parameters
CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-5;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
alpha     =  0.9;                 % iterative step size parameter
elle      =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
ells      =  d0*10;               % correlation length for phase fluctuation diffusivity (multiple of d0, 10-20)
gamma     =  1e-3;                % artificial horizontal inertia parameter (only applies if periodic)
bnd_w     =  max(ells/2,h);       % width of boundary layer [m]


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

