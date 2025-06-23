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

chiw   = (chi(icz(1:end-1),:)+chi(icz(2:end),:))./2;
 muw   = ( mu(icz(1:end-1),:)+ mu(icz(2:end),:))./2;

xw     = (x(icz(1:end-1),:)+x(icz(2:end),:))./2;
mw     = (m(icz(1:end-1),:)+m(icz(2:end),:))./2;

Xw     = (X(icz(1:end-1),:)+X(icz(2:end),:))/2;
Mw     = (M(icz(1:end-1),:)+M(icz(2:end),:))/2;

% update lithostatic pressure
Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
Pt          = max(Ptop/100,Pl + P(2:end-1,2:end-1));

% % update effective constituent sizes
% dm = d0.*(1-mu ).^0.5;
% dx = d0.*(1-chi).^0.5;

% get coefficient contrasts
kv = [etax0;etam0];
Mv = [etax0;etam0].'./[etax0;etam0];

% get permission weights
ff = max(eps^0.5,min(1-eps^0.5,permute(cat(3,chi,mu ),[3,1,2])));
FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
Kv   = ff.*kv.*thtv;
Cv   = Kv./d0.^2;

% get effective viscosity
eta0 = squeeze(sum(Kv,1));

% get segregation drag cofficient
Cx0 = squeeze(Cv(1,:,:));

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
        + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);                   % convection speed magnitude
vx = abs(wx([2,2:end-1],2:end-1)+wx(2:end,2:end-1))/2;                         % segregation speed magnitude

% update diffusion parameters
ke   = eII.*Delta_cnv.^2;                                                  % turbulent eddy diffusivity

fRe  = (1-exp(-Re./Rec)+eps);                                              % ramp-up factor for eddy diffusivity
etae = ke.*rho;                                                            % eddy viscosity
eta  = (eta + eta0 + fRe.*etae)/2;                                         % effective viscosity

ks   = vx.*Delta_sgr.*hasx;                                                % segregation diffusivity
kx   = (ks + ke.*fRe/Scx);                                                 % regularised particle diffusivity 

fRex = (1-exp(-Rex./Rexc)+eps);                                            % ramp-up factor for turbulent drag coefficient
Cxt  = chi.*rho.*ks./d0^2;                                                 % turbulent drag coefficient
Cx   = (Cx + Cx0 + fRex.*Cxt)/2;                                           % effective drag coefficient   

% limit total contrast in Cx
Cxmax = geomean(Cx(:)).*(Cxcntr/2);
Cxmin = geomean(Cx(:))./(Cxcntr/2);
Cx    = 1./(1./Cxmax + 1./Cx) + Cxmin;

% interpolate to staggered nodes
Cxw   = (Cx(icz(1:end-1),:)+Cx(icz(2:end),:))/2;

% limit total contrast in eta
etamax = geomean(eta(:)).*(etacntr/2);
etamin = geomean(eta(:))./(etacntr/2);
eta    = 1./(1./etamax + 1./eta) + etamin;

% interpolate to staggered nodes
etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

% update dimensionless numbers
etasgr = Cx.*d0^2./chi;
Re  = V .*Delta_cnv./(eta   ./rho);                                        % Reynolds number on correlation length scale
ReD = V .*D/10     ./(eta   ./rho);                                        % Reynolds number on domain length scale
Rex = vx.*d0       ./(etasgr./rho);                                        % particle Reynolds number
Ra  = V .*Delta_cnv./kx;                                                   % Rayleigh number on correlation length scale
RaD = V .*D/10     ./kx;                                                   % Rayleigh number on domain length scale
Rux = vx./V;                                                               % particle settling number

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% record timing
UDtime = UDtime + toc;
