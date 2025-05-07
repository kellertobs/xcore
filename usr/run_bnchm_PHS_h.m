% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_default')

% set run parameters
runID    =  'bnchm_PHS_h';        % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  0;                   % switch on to live plot of results

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  nop;                 % number of time steps to take
dt       =  1;                   % set initial time step

% set initial crystallinity parameters
x0       =  0.01;                % background crystallinity initial value [wt]
dx0      =  0.00;                % background crystallinity random perturbation [wt]
dxg      =  0.10;                % background crystallinity gaussian perturbation [wt]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
atol     =  1e-12;               % outer its absolute tolerance
rtol     =  atol/1e6;            % outer its absolute tolerance
maxit    =  100;                 % maximum outer its
alpha    =  0.80;                % iterative step size parameter


% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

NN    = 40*[1,2,4];
nshft = 1;

dt     =  D/NN(3)/100;
dtmax  =  D/NN(3)/100;

Nt     =  nshft*D/NN(1)/dt;           % number of time steps to take

for Ni = NN
    
    N     =  Ni;                % number of grid points in z-direction
    h     =  D/N;               % grid spacing

    Delta_cnv = h;

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0;  Wm(:) = 0;  Wx(:) = 0;  wx(:) = 0;  wm(:) = 0;
    U(:) = 0;  Um(:) = 1;  Ux(:) = 1;  upd_wx(:) = 0;
    P(:) = 0;

    % set diffusion parameters to zero to isolate advection
    kwx(:) = 0;  kx(:) = 0;  ke(:) = 0;

    % set parameters for non-dissipative, non-reactive flow
    rhoin = rho; rhoout = circshift(rho,Ni/NN(1)*nshft,2);
    Min   = M;   Mout   = circshift(M  ,Ni/NN(1)*nshft,2);
    Xin   = X;   Xout   = circshift(X  ,Ni/NN(1)*nshft,2);

    time  = 0;

    output;

    % physical time stepping loop
    while time <= tend && step <= Nt

        % time step info
        timing;

        % store previous solution
        store;

        % reset residuals and iteration count
        resnorm  = 1;
        resnorm0 = resnorm;
        iter     = 1;
        if frst; alpha = alpha/2; beta = beta/2; end

        % non-linear iteration loop
        while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit

            % solve thermo-chemical equations
            phsevo;

            % update non-linear parameters and auxiliary variables
            update;

            kwx(:) = 0;  kx(:) = 0;  ke(:) = 0;

            % report convergence
            report;

            iter = iter+1;
        end

        % print model diagnostics
        diagnose;

        % plot model results
        if ~mod(step,nop); output; end

        % increment time/step
        time = time+dt;
        step = step+1;
        if frst; alpha = alpha*2; beta = beta*2; frst=0; end

        figure(100); clf;
        plot(XX(ceil(Nz/2),:), Xout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/4),:), X(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
        xlabel('Distance [m]','Interpreter','latex','FontSize',16)
        ylabel('Crystallinity [wt]','Interpreter','latex','FontSize',16)
        drawnow;

    end

    % plot convergence
    EB = norm(rho-rhoout,'fro')./norm(rhoout,'fro');
    EM = norm(  M-  Mout,'fro')./norm(rhoout,'fro');
    EX = norm(  X-  Xout,'fro')./norm(rhoout,'fro');

    clist = [colororder;[0 0 0]];

    fh15 = figure(15);
    loglog(h,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    loglog(h,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    loglog(h,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Grid step [m]','Interpreter','latex','FontSize',15)
    ylabel('Rel. numerical error [1]','Interpreter','latex','FontSize',15)
    title('Numerical convergence in space','Interpreter','latex','FontSize',18)

    if Ni == NN(1)
        % loglog(D./NN,geomean([EB,EM,EX]).*(NN./NN(1)).^-1,'k-.','LineWidth',2);  % plot trend for comparison
        % loglog(D./NN,geomean([EB,EM,EX]).*(NN./NN(1)).^-3,'k--','LineWidth',2);  % plot trend for comparison
        loglog(D./NN,geomean([EB,EM,EX]).*(NN./NN(1)).^-5,'k-','LineWidth',2);  % plot trend for comparison
    end
    if Ni == NN(end)
        legend({'error $\bar{\rho}$','error $M$','error $X$','quintic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

    clear DD GG KP;

end

name = [opdir,'/',runID,'/',runID,'_',ADVN];
print(fh15,name,'-dpng','-r300','-vector');