%% *****  MODEL INITIALISATION  *******************************************

if ~postprc
% create output directory
if ~isfolder([outdir,'/',runID])
    mkdir([outdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1
    parfile = [outdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN REGC MODEL | %s  ********************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);
end

% load and process custom colormaps
switch colourmap
    case 'ocean'
        load ./colmap/ocean.mat
        colmap = ocean;   
    case 'batlow'
        load ./colmap/batlow.mat 
        colmap = batlow;
    case 'batlowW'
        load ./colmap/batlowW.mat
        colmap = batlowW;
    case 'batlowK'
        load ./colmap/batlowK.mat
        colmap = batlowK;
    case 'lipari'
        load ./colmap/lipari.mat
        colmap = lipari;  
    case 'glasgow'
        load ./colmap/glasgow.mat;
        colmap = glasgow; 
    case 'navia'
        load ./colmap/navia.mat;
        colmap = navia;  
    case 'lapaz'
        load ./colmap/lapaz.mat;
        colmap = lapaz;  
    case 'lajolla'
        load ./colmap/lajolla.mat;
        colmap = lajolla;  
end
nclmp = length(colmap);

BCA     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)
BCD     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

% get smoothed initialisation field
rng(seed);
smth = smth*Nx*Nz*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth)
    rp = rp + diffus(rp,1/8*ones(size(rp)),1,[1,2],BCD);
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for transient boundary layers
topshape = zeros(size(ZZ));
botshape = zeros(size(ZZ));
sdsshape = zeros(size(XX));
if ~any(bnd_h)
    switch bndmode
        case 0  % none
        case 1  % top only
            topshape = exp( ( -ZZ)/max(h,bnd_w));
        case 2  % bot only
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
        case 3  % top/bot only
            topshape = exp( ( -ZZ)/max(h,bnd_w));
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
        case 4 % all walls
            topshape = exp( ( -ZZ)/max(h,bnd_w));
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
            sdsshape = exp( ( -XX)/max(h,bnd_w)) ...
                     + exp(-(L-XX)/max(h,bnd_w));
        case 5 % only walls
            sdsshape = exp( ( -XX)/max(h,bnd_w)) ...
                     + exp(-(L-XX)/max(h,bnd_w));
    end
    sdsshape = max(0,sdsshape - topshape - botshape);
end

bnd_X = zeros(Nz,Nx);

% set specified boundaries to no slip, else to free slip
sds = -1;
top =  1;
bot = -1;

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [1,1:Nz,Nz];
ifx = [Nx,1:Nx+1,2];
ifz = [2,1:Nz+1,Nz];

% initialise crystallinity field
x   =  x0 + dx0.*rp + topshape.*(xb + dxb.*rp);  % potential temperature [C]
m   =  1-x;
xin =  x;

U   =  zeros(Nz+2,Nx+1);  UBG = U; upd_U = 0*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; wx = 0.*W; wx0 = 0.*W; wxo = wx; upd_W = 0*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*x; upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

Delta_cnv0 = Delta_cnv;
Re     = eps;  
Div_V  = 0.*x;  advn_rho = 0.*x;  advn_X = 0.*x; advn_M = 0.*x; drhodt = 0.*x;  drhodto = drhodt;
exx    = 0.*x;  ezz = 0.*x;  exz = zeros(Nz-1,Nx-1);  eII = 0.*x;  
txx    = 0.*x;  tzz = 0.*x;  txz = zeros(Nz-1,Nx-1);  tII = 0.*x; 
eta    = etam0 + zeros(Nz,Nx);
etamax = min(eta(:)) .* etacntr;
dV     = 0.*x; 
ke     = 0.*x;
Pref   = 1e5;
rho    = (x./rhox0 + m./rhom0).^-1;
Pt     = Ptop + rho.*g0.*ZZ;  Pl = Pt;  Pto = Pt; Ptoo = Pt;
rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhow.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;

% get volume fractions and bulk density
step    = 0;
FMtime  = 0;
XEtime  = 0;
UDtime  = 0;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

update;

X    = rho.*x;  Xo = X;  res_X = 0.*X;
M    = rho.*m;  Mo = M;  res_M = 0.*M;

Gx   =  Da.*(xeq-x).*rho./dt.*topshape;
Gm   =  -Gx;

xo   = x;
mo   = m;
Xo   = X;
Mo   = M;
rhoo = rho;
dto  = dt; 

% initialise correlation length for convective/turbulent regularisation
corrl;

% initialise auxiliary variables 
dwxdt   = 0.*wx; dwxdto = dwxdt;  advn_wx = 0.*wx;
dXdt    = 0.*x;  dXdto  = dXdt;
dMdt    = 0.*m;  dMdto  = dMdt;
upd_X   = 0.*X;
upd_M   = 0.*M;
upd_wx  = 0.*wx;
upd_rho = 0.*rho;

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumBdt = 0; dsumBdto = 0;
dsumMdt = 0; dsumMdto = 0;
dsumXdt = 0; dsumXdto = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','x','m','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','eII','tII','dt','time','step','dV','wx','dwxdt');
        name = [outdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);P(:)];
        RHO = X+M;

        update;
        corrl;
        update;
        store;
        fluidmech;
        update;
        output;

        time    = time+dt;
        step    = step+1;
        restart = 0;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        restart = 0;
        store;
        fluidmech;
        update;
        history;
        output;
        step = step+1;
    end
else
    % complete, plot, and save initial condition
    store;
    fluidmech;
    update;
    % history;
    output;
    step = step+1;
end

