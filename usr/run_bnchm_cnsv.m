% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% test decreasing time step
ATOL = [1e-6,1e-9,1e-12];

for atol = ATOL

    % set run parameters
    runID    =  'bnchm_cnsv';        % run identifier
    restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop      =  10;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on to live plot of results
    plot_cv  =  1;                   % switch on to live plot iterative convergence

    % set model domain parameters
    N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)

    % set model timing parameters
    Nt       =  2*nop;              % number of time steps to take
    dt       =  0.1;                 % set initial time step

    % set initial crystallinity parameters
    x0       =  0.01;                % background crystallinity initial value [wt]
    dx0      =  0.00;                % background crystallinity random perturbation [wt]
    dxg      =  0.10;                % background crystallinity gaussian perturbation [wt]

    % set numerical model parameters
    TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL       =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
    rtol      =  atol/1e6;            % outer its absolute tolerance
    maxit     =  100;                 % maximum outer its
    alpha     =  0.75;                % iterative step size parameter
    Delta_cnv = h;

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    EB = norm(diff(hist.EB(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EB(Nt/2:Nt))         ));
    EM = norm(diff(hist.EM(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EM(Nt/2:Nt))         ));
    EX = norm(diff(hist.EX(Nt/2:Nt  ))./diff(hist.time(Nt/2:Nt)),'fro')./sqrt(length(diff(hist.EX(Nt/2:Nt))         ));

    clist = [colororder;[0 0 0]];

    fh21 = figure(21);
    loglog(atol,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    loglog(atol,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Abs. residual tolerance [1]','Interpreter','latex','FontSize',15)
    ylabel('Rel. conservation error rate [1/s]','Interpreter','latex','FontSize',15)
    title('Global conservation with nonlinear convergence','Interpreter','latex','FontSize',18)

    if atol == ATOL(1)
        loglog(ATOL,eps.*ones(size(ATOL)),'k:' ,'LineWidth',2);  % plot trend for comparison
        legend('error $\bar{\rho}$','error $M$','error $X$','machine prec.','Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT,'_',ADVN];
print(fh21,name,'-dpng','-r300','-vector');