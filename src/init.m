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

% calculate and print characteristic scales
D0      =  D/10;
d0      =  d0;
L0      =  elle;
l0      =  ells;
h0      =  h;

rho0    =  rhom0;
Drho0   =  rhox0-rhom0;
Dchi0   =  xeq/2;
eta0    =  etam0;
C0      =  Dchi0*eta0/d0^2;

W0      =  1./(2.*L0.*rho0).*(sqrt(4.*Dchi0.*Drho0.*g0.*L0.*rho0.*D0.^2 + eta0.^2) - eta0);
w0      =  1./(2.*l0.*rho0).*(sqrt(4.*       Drho0.*g0.*l0.*rho0.*d0.^2 + eta0.^2) - eta0);
tW0     =  D0/W0;
tw0     =  D0/w0;
t0      =  D0/(W0 + w0);
p0      =  W0*eta0/D0;

ke0     =  W0*L0;
ks0     =  w0*l0;
dt0     =  min([(h0/2)^2/(ks0+ke0) , (h0/2)/(W0+w0)]);

xie0    =  Xi*sqrt(ke0.*(L0./(L0+h0)).^3./((L0+h0)/2/W0)+dt0);
xis0    =  Xi*sqrt(ks0.*(l0./(l0+h0)).^3./((L0+h0)/2/w0)+dt0);
xi0     =  (xie0 + xis0);

tau0    =  h0./(W0 + w0 + eps) + dt0;
G0      =  R.*Dchi0.*rho0./tau0;

etae0   =  ke0*rho0;
etas0   =  ks0*rho0;
Da0     =  G0/(rho0/t0);
Ns0     =  (xie0 + xis0)/(W0 + w0);
Rs0     =  w0/W0;
Ra0     =  W0*D0/(ks0 + ke0);
ReD0    =  W0*D0/((eta0+etae0)/rho0);
Red0    =  w0*d0/((eta0+etas0)/rho0);

fprintf(1,'\n  Scaled domain depth D0    = %1.0e [m]',D0);
fprintf(1,'\n  Crystal size        d0    = %1.0e [m]',d0);
fprintf(1,'\n  Eddy  corrl. length L0    = %1.0e [m]',L0);
fprintf(1,'\n  Segr. corrl. length l0    = %1.0e [m] \n',l0);

fprintf(1,'\n  Density             rho0  = %1.0f  [kg/m3]',rho0);
fprintf(1,'\n  Density contrast    Drho0 = %1.0f   [kg/m3]',Drho0);
fprintf(1,'\n  Cristal. contrast   Dchi0 = %1.3f [wt]',Dchi0);
fprintf(1,'\n  Viscosity           eta0  = %1.0e [Pas]',eta0);
fprintf(1,'\n  Drag coefficient    Cx0   = %1.0e [Pas/m2] \n',C0);

fprintf(1,'\n  Eddy  diffusivity   ke0   = %1.1e [m2/s]',ke0);
fprintf(1,'\n  Segr. diffusivity   ks0   = %1.1e [m2/s]',ks0);
fprintf(1,'\n  Eddy  viscosity     etae  = %1.1e [Pas]',etae0);
fprintf(1,'\n  Segr. viscosity     etas0 = %1.1e [Pas] \n',etas0);

fprintf(1,'\n  Eddy noise rate     xie0  = %1.2e [m/s]',xie0);
fprintf(1,'\n  Segr. noise rate    xis0  = %1.2e [m/s]\n',xis0);

fprintf(1,'\n  Reaction rate       G0    = %1.2e [kg/m3/s]\n',G0);

fprintf(1,'\n  Convection  speed   W0    = %1.2e [m/s]',W0);
fprintf(1,'\n  Segregation speed   w0    = %1.2e [m/s]\n',w0);

fprintf(1,'\n  Convection  time    tW0   = %1.2e [s]',tW0);
fprintf(1,'\n  Segregation time    tw0   = %1.2e [s] \n',tw0);

fprintf(1,'\n  Dahmköhler No       Da0   = %1.2e [1]',Da0);
fprintf(1,'\n  Noise No            Ns0   = %1.2e [1]',Ns0);
fprintf(1,'\n  Segregation No      Rs0   = %1.2e [1]',Rs0);
fprintf(1,'\n  Rayleigh No         Ra0   = %1.2e [1]',Ra0);
fprintf(1,'\n  Domain  Reynolds No ReD0  = %1.2e [1]',ReD0);
fprintf(1,'\n  Crystal Reynolds No Red0  = %1.2e [1]\n\n\n',Red0);


BCA     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)
BCD     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)
open    = 1-closed;

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
% smth = max(4,max((Delta_sgr/2/h)^2,(Delta_cnv/2/h)^2));   % random perturbation smoothness
% rng(seed);
% rp   = randn(Nz,Nx);
% for i = 1:ceil(smth)
%     ksmth = min(1,smth-i+1);
%     rp = rp + diffus(rp,ksmth/8*ones(size(rp)),1,[1,2],{'periodic','periodic'});
% end
% rp  = (rp-mean(rp(:)))./std(rp(:));
% 
% gp = exp(-(XX-L/2  ).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2) ...
%    + exp(-(XX-L/2+L).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2) ...
%    + exp(-(XX-L/2-L).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2);

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

kw2 = kwx.^2 + kwz.^2;
ku2 = kux.^2 + kuz.^2;
kp2 = kpx.^2 + kpz.^2;

% Gaussian spatial filter in Fourier space
Gkwe = exp(-0.5 * (((elle+h)/2)^2) * kw2);
Gkue = exp(-0.5 * (((elle+h)/2)^2) * ku2);
Gkws = exp(-0.5 * (((ells+h)/2)^2) * kw2);
Gkus = exp(-0.5 * (((ells+h)/2)^2) * ku2);
Gkp  = exp(-0.5 * (((elle+ells+h)/2)^2) * kp2);

% Generate new white noise
rwe = randn(Nz+1, Nx+0);
rue = randn(Nz+0, Nx+1);
rws = randn(Nz+1, Nx+0);
rus = randn(Nz+0, Nx+1);
rp  = randn(Nz  , Nx  );

rwe = fft2(rwe);
rue = fft2(rue);
rws = fft2(rws);
rus = fft2(rus);
rp  = fft2(rp );

% Filter white noise spatially
rwe = real(ifft2(Gkwe .* rwe));
rue = real(ifft2(Gkue .* rue));
rws = real(ifft2(Gkws .* rws));
rus = real(ifft2(Gkus .* rus));
rp  = real(ifft2(Gkp  .* rp ));

% Rescale to unit standard deviation
rwe = (rwe - mean(rwe(:))) / std(rwe(:)) .* (1-exp((-ZZw(:,2:end-1))/max(h,elle/2)) - closed.*exp(-(D-ZZw(:,2:end-1))/max(h,elle/2)));
rue = (rue - mean(rue(:))) / std(rue(:));
rws = (rws - mean(rws(:))) / std(rws(:)) .* (1-exp((-ZZw(:,2:end-1))/max(h,ells/2)) - closed.*exp(-(D-ZZw(:,2:end-1))/max(h,ells/2)));
rus = (rus - mean(rus(:))) / std(rus(:));
rp  = (rp  - mean(rp (:))) / std(rp (:));

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for transient boundary layers
topshape = exp( ( -ZZ+h/2)/bnd_w);
botshape = exp( (D-ZZ+h/2)/bnd_w);

% set specified boundaries to no slip, else to free slip
sds = -1;
top =  1;
bot = -1;
if closed; bot = 1; end

% set ghosted index arrays
icx = [Nx,1:Nx,1];
icz = [1,1:Nz,Nz];
ifx = [Nx,1:Nx+1,2];
ifz = [2,1:Nz+1,Nz];

% initialise crystallinity field
x   =  x0 + dx0.*rp + topshape.*(xb + dxb.*rp);
m   =  1-x;
xin =  x;

U   =  zeros(Nz+2,Nx+1);  UBG = U; upd_U = 0*U; 
W   =  zeros(Nz+1,Nx+2);  WBG = W; wx = 0.*W; wm = 0.*W; wx0 = 0.*W; wxo = wx; upd_W = 0*W; Mx = 0*wx(:,2:end-1); 
P   =  zeros(Nz+2,Nx+2);  V   = 0.*x; vx = V; vxo = vx; upd_P = 0*P;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wx  = W;  Ux = U;  wx = W;  xiw = W;
Wm  = W;  Um = U;  ux = U;  xiu = U;

Delta_cnv0 = Delta_cnv;
Re     = eps + 0.*x;  
Rex    = eps + 0.*x;
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

% initialise correlation length for convective/turbulent regularisation
corrl;

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
        load(name,'U','W','P','Pt','x','m','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','etas','eII','tII','ke','ks','kx','RaD','ReD','Rs','Red','dt','time','step','MFS','wm','wx','Mx','dMxdt');
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

