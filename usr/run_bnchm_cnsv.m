% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% test decreasing time step
<<<<<<< Updated upstream
ATOL = [1e-3,1e-6,1e-9];
=======
ATOL = [1e-2,1e-4,1e-6];
>>>>>>> Stashed changes

for atol = ATOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop      =  10;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  1;                   % switch on to live plot iterative convergence

    % set model domain parameters
    D         =  1e1;                 % chamber depth [m]
    N         =  100;                  % number of grid points in z-direction
    h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
    L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

    % set model timing parameters
    Nt        =  5*nop;               % stop when dimensionless time is reached

    % set crystallinity initial condition
    xeq       =  0.01;                % equilibrium crystallinity of boundary layer [wt]
    x0        =  xeq/10;              % initial background crystallinity [wt]
    dxr       =  x0/10;                   % initial random perturbation [wt]
    dxg       =  0;                % initial gaussian perturbation [wt]

    % set physical control parameters
    d0        =  1e-2;                % xtal size constant [m]
    etam0     =  1e+1;                % melt viscosity constant [kg/m3]
    L0        =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
    l0        =  d0*10;               % correlation length for phase fluctuation diffusivity (multiple of d0, 10-20)
<<<<<<< Updated upstream
    R         =  1.0;                 % relative amplitude of crystallisation rate [s]
    Xi        =  1.0;                 % relative amplitude of random noise flux
    open_sgr  =  1;
    open_cnv  =  0;

    % set numerical model parameters
    CFL       =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
    rtol      =  atol/1e6;            % outer its relative tolerance
=======
    R         =  1;
    Xi        =  1;

    % set numerical model parameters
    TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL       =  1/2;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
    rtol      =  atol/1e6;            % outer its absolute tolerance
>>>>>>> Stashed changes
    maxit     =  100;                 % maximum outer its
    alpha     =  0.8;                 % iterative step size parameter
    gamma     =  1e-3;                % artificial horizontal inertia parameter (only applies if periodic)

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    EB = rms(diff(HST.EB(Nt/2:Nt))./diff(HST.time(Nt/2:Nt)));
    EM = rms(diff(HST.EM(Nt/2:Nt))./diff(HST.time(Nt/2:Nt)));
    EX = rms(diff(HST.EX(Nt/2:Nt))./diff(HST.time(Nt/2:Nt)));

    clist = [colororder;[0 0 0]];

    fh24 = figure(24);
    loglog(atol,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    loglog(atol,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Abs. residual tolerance [1]','Interpreter','latex','FontSize',15)
    ylabel('Rel. conservation error rate [1/s]','Interpreter','latex','FontSize',15)
    title('Global conservation with nonlinear convergence','Interpreter','latex','FontSize',18)

    if atol == ATOL(1)
        loglog(ATOL,eps.*ones(size(ATOL)),'k:' ,'LineWidth',2);  % plot trend for comparison
    end
    if atol == ATOL(end)
        legend('error $\bar{\rho}$','error $M$','error $X$','machine prec.','Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT,'_',ADVN];
print(fh24,name,'-dpng','-r300','-vector');