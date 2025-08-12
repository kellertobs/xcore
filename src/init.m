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
open    =  1-closed;

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
dt = dt0/10;

% initialise smooth random noise generation
rng(seed);

% Wavenumber grid
[kwx, kwz] = ndgrid( ...
    2*pi*ifftshift((0:Nz+0) - floor((Nz+1)/2)) / ((Nz+1)*h), ...
    2*pi*ifftshift((0:Nx-1) - floor((Nx+0)/2)) / ((Nx+0)*h)  );

[kux, kuz] = ndgrid( ...
    2*pi*ifftshift((0:Nz-1) - floor((Nz+0)/2)) / ((Nz+0)*h), ...
    2*pi*ifftshift((0:Nx+0) - floor((Nx+1)/2)) / ((Nx+1)*h)  );

[kpx, kpz] = ndgrid( ...
    2*pi*ifftshift((0:Nz-1) - floor(Nz/2)) / (Nz*h), ...
    2*pi*ifftshift((0:Nx-1) - floor(Nx/2)) / (Nx*h)  );

[kcx, kcz] = ndgrid( ...
    2*pi*ifftshift((0:Nz-0) - floor((Nz+1)/2)) / ((Nz+1)*h), ...
    2*pi*ifftshift((0:Nx-0) - floor((Nx+1)/2)) / ((Nx+1)*h)  );

kw2 = kwx.^2 + kwz.^2;
ku2 = kux.^2 + kuz.^2;
kp2 = kpx.^2 + kpz.^2;
kc2 = kcx.^2 + kcz.^2;

% Gaussian spatial filter in Fourier space
Gkwe = exp(-0.5 * ((elle/2+h)^2) * kw2);
Gkue = exp(-0.5 * ((elle/2+h)^2) * ku2);
Gkws = exp(-0.5 * ((ells/2+h)^2) * kw2);
Gkus = exp(-0.5 * ((ells/2+h)^2) * ku2);
Gkps = exp(-0.5 * ((elle  +h)^2) * kc2);
Gkrp = exp(-0.5 * ((elle/2+ells/2+h)^2) * kp2);

% Generate new white noise
rp  = randn(Nz  , Nx  );
rp  = fft2(rp );

% Filter white noise spatially
rp  = real(ifft2(Gkrp  .* rp ));

% Rescale to unit standard deviation
rp  = (rp  - mean(rp (:))) / std(rp (:));

% initialise noise flux variables
xisw  = zeros(Nz+1, Nx+2);
xisu  = zeros(Nz+2, Nx+1);
xiew  = zeros(Nz+1, Nx+2);
xieu  = zeros(Nz+2, Nx+1);
xiexw = zeros(Nz+1, Nx+2);
xiexu = zeros(Nz+2, Nx+1);

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for transient boundary layers
bnd_w     =  max(ells,2*h);         % width of boundary layer [m]
initshape = exp((-ZZ+h/2)/bnd_w);%1./(1+exp((ZZ-1.5*bnd_w-h/2)/bnd_w*6));
bndshape  = exp((-ZZ+h/2)/h);

% set specified boundaries to no slip, else to free slip
sds = -1;
top =  1;
bot = -1;
if closed; bot = 1; bot_sgr = 1; end
if open_sgr; bot_sgr = -1; end

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [1,1:Nz,Nz];
ifx = [Nx,1:Nx+1,2];
ifz = [2,1:Nz+1,Nz];

% initialise crystallinity field
xin  = initshape.*xeq + (1-initshape).*x0;
x   =  xin .* (1+1/10.*rp);
m   =  1-x;

U   =  zeros(Nz+2,Nx+1);  UBG = U; upd_U = 0*U; 
W   =  zeros(Nz+1,Nx+2);  WBG = W; wx = 0.*W; wm = 0.*W; wx0 = 0.*W; wxo = wx; upd_W = 0*W; Mx = 0*wx(:,2:end-1); 
P   =  zeros(Nz+2,Nx+2);  V   = 0.*x; vx = V; vxo = vx; upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux = U;  wx = W;  xiw = W;
Wm  = W;  Um = U;  ux = U;  xiu = U;

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
rhoWo  = rhow.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;

% get volume fractions and bulk density
step    = 0;
FMtime  = 0;
XEtime  = 0;
UDtime  = 0;
dto     = dt;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

X    = rho.*x;  Xo = X;  res_X = 0.*X;
M    = rho.*m;  Mo = M;  res_M = 0.*M;

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
dwxdt   = 0.*wx; dwxdto = dwxdt;  advn_Mx = 0.*wx;
dXdt    = 0.*x;  dXdto  = dXdt;
dMdt    = 0.*m;  dMdto  = dMdt;
dMxdt   = 0.*Mx; dMxdto = dMxdt;
upd_X   = 0.*X;
upd_M   = 0.*M;
upd_Mx  = 0.*Mx;
upd_rho = 0.*rho;
tau_p   = 1;
relax   = 0;

% initialise timing and iterative parameters
frst    = 1;
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
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','dt','time','step','MFS','wx','xisw','xisu','xiew','xieu','xiexw','xiexu');
        name = [outdir,'/',runID,'/',runID,'_HST'];
        load(name,'HST');

        SOL = [W(:);U(:);P(:)];

        update;
        corrl;
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

