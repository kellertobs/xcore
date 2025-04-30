 % prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID     =  'DEMO_2D';           % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  10;                   % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file
colourmap = 'lapaz';             % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set model domain parameters
D         =  10;                   % chamber depth [m]
N         =  100;                  % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  1e6;                 % number of time steps to take
tend      =  1*yr;                % end time for simulation [s]
dt        =  10;                  % initial time step [s]

% set initial thermo-chemical state
smth      =  25;
x0        =  0.01;               % initial temperature [deg C]
dx0       =  5e-3;
xb        =  0.0;
dxb       =  0e-3;
rhom0     =  2500;
rhox0     =  3000;
etam0     =  100;
etax0     =  1e18;
d0        =  0.03;

AA      =[ 0.5989, 0.1772; ...    % permission slopes
           0.0397, 0.1182 ];      % increases permission slopes away from step function 

BB      =[ 0.6870, 0.3130; ...  % permission step locations
           0.9998, 0.0002;];        % sets midpoint of step functions

CC      =[[0.9826, 0.0174]*9.1697; ... % permission step widths
          [0.1695, 0.8305]*4.2773;];   % factor increases width of step functions

% set thermo-chemical boundary parameters
bnd_w     =  D/100;
xeq       =  0.1;
Da        =  0.01;                   % wall cooling/assimilation time [s]
Ptop      =  1e5;                 % top pressure [Pa]

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
alpha     =  0.75;                % iterative step size parameter
rtol      =  1e-3;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
Delta_cnv =  h;                % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  10*dx0;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
etamin    =  1e2;
etacntr   =  1e8;
Rer       =  1;
Scx       =  1;
kmin      =  1e-12;

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

