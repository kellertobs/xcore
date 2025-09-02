%*****  MODEL INITIALISATION  *********************************************

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
fprintf('*****  RUN REGC MODEL | %s  ***************\n',datetime('now'));
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
lncls = colmap([5,135],:);

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
[XXc,ZZc] = meshgrid(Xf,Zf);

Nx = length(Xc);
Nz = length(Zc);

% get characteristic scales
scales;
if ~exist('dt','var'); dt = dt0/10; end


% Wavenumber grid
[kpx, kpz] = ndgrid( ...
    2*pi*ifftshift((0:Nz-1) - floor(Nz/2)) / (Nz*h), ...
    2*pi*ifftshift((0:Nx-1) - floor(Nx/2)) / (Nx*h)  );

padL0 = 4*ceil(L0/h);
padl0 = 4*ceil(l0/h);

[kpx_padL0, kpz_padL0] = ndgrid( ...
    2*pi*ifftshift((0:Nz-1+padL0) - floor((Nz+padL0)/2)) / ((Nz+padL0)*h), ...
    2*pi*ifftshift((0:Nx-1      ) - floor( Nx       /2)) / ( Nx       *h)  );

[kpx_padl0, kpz_padl0] = ndgrid( ...
    2*pi*ifftshift((0:Nz-1+padl0) - floor((Nz+padl0)/2)) / ((Nz+padl0)*h), ...
    2*pi*ifftshift((0:Nx-1      ) - floor( Nx       /2)) / ( Nx       *h)  );

kp2 = kpx.^2 + kpz.^2;
kp2_padL0 = kpx_padL0.^2 + kpz_padL0.^2;
kp2_padl0 = kpx_padl0.^2 + kpz_padl0.^2;

% Gaussian spatial filter kernels in Fourier space
Gkps = exp(-0.5 * ((l0h*sqrt(2))^2) * kp2);
Gkpe = exp(-0.5 * ((L0h*sqrt(2))^2) * kp2);
Gkps_padl0 = exp(-0.5 * ((l0h*sqrt(2))^2) * kp2_padl0);
Gkpe_padL0 = exp(-0.5 * ((L0h*sqrt(2))^2) * kp2_padL0);
Gkrp = exp(-0.5 * ((L0h+l0h    )^2) * kp2);

% initialise smooth random noise generation
rng(seed);

% Generate new white noise
rp  = randn(Nz, Nx);

% Filter white noise spatially
rp  = real(ifft2(Gkrp .* fft2(rp)));

% Rescale to unit standard deviation
rp  = (rp - mean(rp(:))) ./ std(rp(:));

% initialise noise flux potentials
psie = zeros(Nz+0, Nx+0);
psix = zeros(Nz+0, Nx+0);
psis = zeros(Nz+0, Nx+0);

% initialise noise flux components
xisw = zeros(Nz+1, Nx+2);
xisu = zeros(Nz+2, Nx+1);
xiew = zeros(Nz+1, Nx+2);
xieu = zeros(Nz+2, Nx+1);
xixw = zeros(Nz+1, Nx+2);
xixu = zeros(Nz+2, Nx+1);

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for initial and transient boundary layers
initshape = exp((-ZZ+h/2)/(L0h+l0h)/2); % width of initial boundary layer [m]
bndshape  = exp((-ZZ+h/2)/(L0h+l0h)/1); % width of crystal replenishing layer [m]

% set specified boundaries to no slip, else to free slip
sds        = -1;
top_cnv    =  1;
bot_cnv    =  1;
if open_cnv; bot_cnv = -1; end

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [1,1:Nz,Nz];
ifx = [Nx,1:Nx+1,2];
ifz = [2,1:Nz+1,Nz];

% initialise crystallinity field
gp  =  exp(-((XX-L/2)./(L/6)).^2) .* exp(-((ZZ-D/2)./(D/6)).^2);
if open_sgr
    xin =  initshape.*xeq + (1-initshape).*x0;
else
    xin = x0.*ones(Nz,Nx);
end
x   =  xin .* (1+dxr/x0.*rp+dxg/x0.*gp);
m   =  1-x;

U   =  zeros(Nz+2,Nx+1);  UBG = U; upd_U = 0*U; 
W   =  zeros(Nz+1,Nx+2);  WBG = W; wx = 0.*W; wm = 0.*W; wx0 = 0.*W; wxo = wx; upd_W = 0*W; Mx = 0*wx(:,2:end-1); 
P   =  zeros(Nz+2,Nx+2);  V   = 0.*x; vx = V; vxo = vx; upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux = U;  wx = W;  xiw = W;  qz_advn_X = W;  qz_advn_M = W;  qz_dffn_X = W;  qz_dffn_M = W;
Wm  = W;  Um = U;  ux = U;  xiu = U;  qx_advn_X = W;  qx_advn_M = W;  qx_dffn_X = W;  qx_dffn_M = W;

Re     = eps + 0.*x;  
Div_V  = 0.*x;  advn_rho = 0.*x;  advn_X = 0.*x; advn_M = 0.*x; drhodt = 0.*x;  drhodto = drhodt; 
xis = 0.*x;  xie = 0.*x;
exx    = 0.*x;  ezz = 0.*x;  exz = zeros(Nz-1,Nx-1);  eII = 0.*x;  
txx    = 0.*x;  tzz = 0.*x;  txz = zeros(Nz-1,Nx-1);  tII = 0.*x; 
eta    = etam0 + zeros(Nz,Nx);
etas   = etam0 + zeros(Nz,Nx);
MFS    = 0.*x; 
ke     = 0.*x;
Pref   = 1e5;
rho    = (x./rhox0 + m./rhom0).^-1;
Pt     = Ptop + rho.*g0.*ZZ;  Pl = Pt;  Pto = Pt; Ptoo = Pt;
rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhow.*W(:,2:end-1); advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); advn_mx = 0.*rhoUo;

% get volume fractions and bulk density
step    = 0;
FMtime  = 0;
XEtime  = 0;
UDtime  = 0;
dto     = dt;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

X = rho.*x;  Xo = X;  res_X = 0.*X;
M = rho.*m;  Mo = M;  res_M = 0.*M;

update;

Gx   =  0.*x;
Gm   =  -Gx;

xo   = x;
mo   = m;
Xo   = X;
Mo   = M;
Mxo  = Mx;
rhoo = rho;
dto  = dt; 

% initialise auxiliary variables 
dXdt    = 0.*x;   dXdto  = dXdt;
dMdt    = 0.*m;   dMdto  = dMdt;
drhodt  = 0.*rho; drhodto = drhodt;
upd_X   = 0.*X;
upd_MFS = 0.*MFS;

% initialise timing, history, and iterative parameters
step    = 0;
time    = 0;
iter    = 0;
HST     = [];
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
        fprintf('\n  restart from %s \n\n',name);
        load(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','dt','time','step','MFS','wx','wm','psie','psix','psis','xisw','xisu','xiew','xieu','xixw','xixu');
        name = [outdir,'/',runID,'/',runID,'_HST'];
        load(name,'HST');

        SOL = [W(:);U(:);P(:)];

        update;
        store;
        output;

        time    = time+dt;
        step    = step+1;
        restart = 0;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        restart = 0;
        store;
        history;
        output;
        step = step+1;
    end
else
    % complete, plot, and save initial condition
    store;
    history;
    output;
    step = step+1;
end

