%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  *****************************

tic;

% update phase indicators
hasx = x >= eps^0.5;
hasm = m >= eps^0.5;

rho    = 1./(m./rhom0  + x./rhox0);

rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;

rhoref = mean(rhow,2);

Drhom  = rhom0 - rhow;
Drhox  = rhox0 - rhow;
Drho   = rhow  - rhoref;

rhoW   = rhow.*W(:,2:end-1);
rhoU   = rhou.*U(2:end-1,:);

% convert weight to volume fraction, update bulk density
chi    = max(eps,min(1-eps, x.*rho./rhox0));
mu     = max(eps,min(1-eps, m.*rho./rhom0));

chiw   = (chi(icz(1:end-1),:)+chi(icz(2:end),:))./2;
 muw   = ( mu(icz(1:end-1),:)+ mu(icz(2:end),:))./2;

xw     = (x(icz(1:end-1),:)+x(icz(2:end),:))./2;
mw     = (m(icz(1:end-1),:)+m(icz(2:end),:))./2;

xu     = (x(:,icx(1:end-1))+x(:,icx(2:end)))./2;
muu    = (m(:,icx(1:end-1))+m(:,icx(2:end)))./2;

Xw     = (X(icz(1:end-1),:)+X(icz(2:end),:))/2;
Mw     = (M(icz(1:end-1),:)+M(icz(2:end),:))/2;

% update lithostatic pressure
Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
Pt          = max(Ptop/100,Pl + P(2:end-1,2:end-1));

% get coefficient contrasts
kv = [etax0;etam0];
Mv = [etax0;etam0].'./[etax0;etam0];

% get permission weights
ff = permute(cat(3,chi,mu ),[3,1,2]);
FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
etai = kv.*thtv;

% get effective viscosity
etax   = squeeze(etai(1,:,:));
etamix = squeeze(sum(ff.*etai,1));

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V/3;                                 % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V/3;                                 % z-normal strain rate
exz = (diff(U,1,1)./h+diff(W,1,2)./h)/2;                                   % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% update velocity magnitudes
V   = sqrt(((W  (1:end-1,2:end-1)+W  (2:end,2:end-1))/2).^2 ...
         + ((U  (2:end-1,1:end-1)+U  (2:end-1,2:end))/2).^2);              % convection speed magnitude
bndtapers = (1 - (exp((-ZZ)/l0) + exp(-(D-ZZ)/l0)).*(1-open_sgr));
vx  = d0^2./etas.*(rhox0-rhom0).*g0.*bndtapers;                            % xtal segregation speed magnitude
vm  = vx.*x./(1-x);                                                        % melt segregation speed magnitude
xis = sqrt(((xisw (1:end-1,2:end-1)+xisw (2:end,2:end-1))/2).^2 ...
         + ((xisu (2:end-1,1:end-1)+xisu (2:end-1,2:end))/2).^2);          % settling noise flux magnitude 
xiex= sqrt(((xiexw(1:end-1,2:end-1)+xiexw(2:end,2:end-1))/2).^2 ...
         + ((xiexu(2:end-1,1:end-1)+xiexu(2:end-1,2:end))/2).^2);          % xtal eddy noise flux magnitude
xie = sqrt(((xiew (1:end-1,2:end-1)+xiew (2:end,2:end-1))/2).^2 ...
         + ((xieu (2:end-1,1:end-1)+xieu (2:end-1,2:end))/2).^2);          % eddy noise flux magnitude
xix = xis + xiex;

% update diffusion parameters
bndtapere = (1 - (exp((-ZZ)/L0) + exp(-(D-ZZ)/L0).*(1-open_cnv)));
ke   = eII.*L0.^2 .* bndtapere;                                            % turbulent eddy diffusivity
ks   = vx .*l0;                                                            % segregation diffusivity
kx   = (ks + fReL*ke);                                                     % regularised particle diffusivity 

% update viscosities
etae = fReL*ke.*rho;                                                       % eddy viscosity
eta  = (eta + etamix + etae)/2;                                            % effective viscosity

etat = fRel.*ks.*rho;                                                      % turbulent drag viscosity
etas = (etas + etax + etat)/2;                                             % effective drag viscosity   

% limit total viscosity contrast
etamax = geomean(eta(:)).*(etacntr/2);
etamin = geomean(eta(:))./(etacntr/2);
eta    = 1./(1./etamax + 1./eta) + etamin;

etamax = geomean(etas(:)).*(etacntr/2);
etamin = geomean(etas(:))./(etacntr/2);
etas   = 1./(1./etamax + 1./etas) + etamin;

% interpolate to staggered nodes
etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

etasw = (etas(icz(1:end-1),:).*etas(icz(2:end),:)).^0.5;

% update dimensionless numbers
ReD = V .*D0./(eta ./rho);                                                 % Reynolds number on scaled domain length
Red = vx.*d0./(etas./rho);                                                 % particle Reynolds number
Ra  = V .*D0./kx;                                                          % Rayleigh number on scale domain length 
Rc  = V./vx;                                                               % particle settling number
Ne  = xie./V;                                                              % eddy noise flux number
Ns  = xix./vx;                                                             % settling noise flux number

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% update time step
dtk = (h/2)^2/max(kx(:)); % diffusive time step size
dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:)]+eps));  % advective time step size
dt  =  min([1.5*dto,min([dtk,CFL*dta]),dtmax]); % time step size

% record timing
UDtime = UDtime + toc;
