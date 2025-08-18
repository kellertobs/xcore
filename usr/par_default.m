% set run parameters
runID     =  'default';           % run identifier
srcdir    =  '../src';            % output directory
outdir    =  '../out';            % output directory
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nrh       =  1;                   % record diagnostic history every 'nrh' time steps
nop       =  100;                 % output frame plotted/saved every 'nop' time steps
ndm_op    =  0;                   % plot nondimensionalised output
plot_op   =  1;                   % switch on to live plot results
save_op   =  0;                   % switch on to save output to file
plot_cv   =  0;                   % switch on to live plot iterative convergence
colourmap = 'lapaz';              % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set unit conversion parameters
hr        =  3600;                % conversion seconds to hours
yr        =  24*365.25*hr;        % conversion seconds to years
cm        =  0.01;                % conversion metre to centimetres
km        =  1000;                % conversion metre to kilometres

% set model domain parameters
D         =  10;                  % chamber depth [m]
N         =  100;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  1e6;                 % number of time steps to take
t0end     =  2;                   % stop when dimensionless time is reached
xend      =  1.00;                % stop run when mean crystallinity reaches threshold
tend      =  10*yr;               % end time for simulation [s]

% set initial phase fraction parameters
xeq       =  0.01;                % equilibrium crystallinity of boundary layer [wt]
x0        =  xeq/10;              % initial background crystallinity [wt]
dxg       =  0;                   % initial gaussian perturbation [wt] (for benchmarking)
seed      =  15;                  % random perturbation seed

% set buoyancy parameters
rhom0     =  2700;                % melt density constant [kg/m3]
rhox0     =  3200;                % xtal density constant [kg/m3]
d0        =  0.01;                % xtal size constant [m]
g0        =  10.;                 % gravity constant [m/s2]

% set rheological parameters
etam0     =  1e1;                 % melt viscosity constant [kg/m3]
etax0     =  1e18;                % xtal viscosity constant [kg/m3]
AA        = [ 0.5989, 0.1772; ...    % permission slopes
              0.0397, 0.1182 ];      % increases permission slopes away from step function 

BB        = [ 0.6870, 0.3130; ...    % permission step locations
              0.9998, 0.0002;];      % sets midpoint of step functions

CC        = [[0.9826, 0.0174]*9.1697; ... % permission step widths
             [0.1695, 0.8305]*4.2773;];   % factor increases width of step functions

% set physical control parameters
L0        =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
l0        =  d0*10;               % correlation length for phase fluctuation diffusivity (multiple of d0, 10-20)
R         =  1.0;                 % relative amplitude of crystallisation rate [s]
Xi        =  1.0;                 % relative amplitude of random noise flux
Ptop      =  1e5;                 % top boundary pressure [Pa]
open_cnv  =  0;                   % switch for open bottom boundary for crystal-driven convection
open_sgr  =  1;                   % switch for open bottom boundary for crystal segregation

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-5;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
alpha     =  0.75;                % iterative step size parameter
gamma     =  0e-3;                % artificial horizontal inertia parameter (only applies if periodic)
kmin      =  1e-16;               % minimum diffusivity
kmax      =  1e+16;               % maximum diffusivity
dtmax     =  1e32;                % maximum time step [s]
etacntr   =  1e+6;                % maximum viscosity contrast

% set other options
bnchm     =  0;                   % not a benchmark run
postprc   =  0;                   % not postprocessing mode
