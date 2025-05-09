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
chi    = max(0,min(1, x.*rho./rhox0));
mu     = max(0,min(1, m.*rho./rhom0));

chiw   = (chi(1:end-1,:)+chi(2:end,:))/2;
 muw   = ( mu(1:end-1,:)+ mu(2:end,:))/2;

xw     = (x(icz(1:end-1),icx)+x(icz(2:end),icx))./2;
mw     = (m(icz(1:end-1),icx)+m(icz(2:end),icx))./2;

Xw     = (X(icz(1:end-1),:)+X(icz(2:end),:))/2;
Mw     = (M(icz(1:end-1),:)+M(icz(2:end),:))/2;

% update lithostatic pressure
Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
Pt          = max(Ptop/100,Pl + P(2:end-1,2:end-1));

% update effective constituent sizes
dm = d0.*(1-mu ).^0.5;
dx = d0.*(1-chi).^0.5;

% get coefficient contrasts
kv = [etax0;etam0];
Mv = [etax0;etam0].'./[etax0;etam0];

% get permission weights
dd = max(eps^0.5,min(1-eps^0.5,permute(cat(3,dx ,dm ),[3,1,2])));
ff = max(eps^0.5,min(1-eps^0.5,permute(cat(3,chi,mu ),[3,1,2])));
FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
Kv   = ff.*kv.*thtv;
Cv   = Kv.*(1-ff)./dd.^2;

% get effective viscosity
eta0 = squeeze(sum(Kv,1));

% get segregation cofficients
Ksgr = ff./Cv;

Cvx  = squeeze(Cv(1,:,:));
Cvm  = squeeze(Cv(2,:,:));

Cvxw = (Cvx(icz(1:end-1),:).*Cvx(icz(2:end),:)).^0.5;
Cvmw = (Cvm(icz(1:end-1),:).*Cvm(icz(2:end),:)).^0.5;

Ksgr_x = squeeze(Ksgr(1,:,:)) + eps^2;
Ksgr_m = squeeze(Ksgr(2,:,:)) + eps^2;


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
V  = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
       + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);                    % convection speed magnitude
vx  = abs(wx(1:end-1,2:end-1)+wx(2:end,2:end-1))/2;                        % segregation speed magnitude

% update diffusion parameters
eII0   = eta0./rho./Delta_cnv0^2;
eIIe   = eII .* (1-exp(-eII./eII0)+eps);
ke     = (ke + eIIe.*Delta_cnv.^2)/2;                                      % turbulent eddy diffusivity
fRe    = (1-exp(-Re./Rec)+eps);
ke  = 1./(1./kmax + 1./ke) + kmin;

kwx = vx.*Delta_sgr.*hasx;                                                 % segregation diffusivity
kx  = (kwx + ke.*fRe/Scx).*chi;                                            % regularised solid fraction diffusion 
eta = ke.*rho.*fRe + eta0;                                                 % regularised momentum diffusion

etamax = etacntr.*max(min(eta(:)),etamin);
eta    = 1./(1./etamax + 1./eta) + etamin;

etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

% update dimensionless numbers
Re     = V .*Delta_cnv./(eta./rho);                                        % Reynolds number on correlation length scale
ReD    = V .*D        ./(eta./rho);                                        % Reynolds number on domain length scale
Rex    = vx.*d0       ./(eta./rho);                                        % particle Reynolds number
Ra     = V.*Delta_cnv./kx;                                                 % Rayleigh number on correlation length scale
RaD    = V.*D        ./kx;                                                 % Rayleigh number on domain length scale
Rux    = vx./V;                                                            % particle settling number

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

UDtime = UDtime + toc;
