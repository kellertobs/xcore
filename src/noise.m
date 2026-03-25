% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter<=1
    % generate white noise
    re = randn(Nz+0, Nx+0);
    rx = randn(Nz+0, Nx+0);
    rs = randn(Nz+0, Nx+0);
end


% noise decorrelation time
taue = L0/2./(V +eps);
taus = l0/2./(vx+eps);
taux = sqrt(taus.*taue);

% particle generation noise relative to xix0, xis0
rg = xix0.*real(ifft2(Gkpx .* fft2(re))) ...
   + xis0.*real(ifft2(Gkps .* fft2(rs)));
rg = (rg - mean(rg(:))) ./ std(rg(:));

% noise flux amplitudes
sge = Xi * sqrt(          fReL.*ke./taue);             % eddy mixture noise speed
sgx = Xi * sqrt(chi.*sqrt(fReL.*ke./taue.*ks./taus));  % eddy crystal noise speed
sgs = Xi * sqrt(chi.*           ks./taus);             % settling noise speed

% Ornstein–Uhlenbeck time update for evolving noise
taue = (1./taue + 1./(1000*dt)).^-1 + dt;
taus = (1./taus + 1./(1000*dt)).^-1 + dt;
taux = (1./taux + 1./(1000*dt)).^-1 + dt;

Fte  = exp(-dt./taue);                                        % eddy noise time evolution factor
Fts  = exp(-dt./taus);                                        % settling noise time evolution factor
Ftx  = exp(-dt./taux);                                        % settling noise time evolution factor
Ftg  = exp(-dt./(xix0.*L0/2/W0+xis0.*l0/2/w0).*(xix0+xis0));  % particle generation noise time evolution factor                        % generation noise time evolution factor

fL   = 2/sqrt(1 - exp(-h^2/(2*L0h^2)));                       % scaling factor for potential field to noise component variance
fl   = 2/sqrt(1 - exp(-h^2/(2*l0h^2)));                       % scaling factor for potential field to noise component variance
fLl  = 2/sqrt(1 - exp(-h^2/(2*Llh^2)));                       % scaling factor for potential field to noise component variance

psie = Fte .* psieo + sqrt(1 - Fte.^2) .* sge.*fL .* re;      % update eddy noise stream function
psix = Ftx .* psixo + sqrt(1 - Ftx.^2) .* sgx.*fLl.* rx;      % update particle eddy noise potential
psis = Fts .* psiso + sqrt(1 - Fts.^2) .* sgs.*fl .* rs;      % update particle settling noise potential
psig = Ftg .* psigo + sqrt(1 - Ftg.^2)            .* rg;      % update particle generation noise

% filter noise amplitudes spatially to decorrelation length
% pad top/bot to avoid periodic bleed over on non-periodic boundaries
psie_padL0 = zeros(Nz+padL0,Nx); 
psix_padLl = zeros(Nz+padLl,Nx);
psis_padl0 = zeros(Nz+padl0,Nx);

psie_padL0(1+padL0/2:end-padL0/2,:) = psie;
psix_padLl(1+padLl/2:end-padLl/2,:) = psix;
psis_padl0(1+padl0/2:end-padl0/2,:) = psis;

psie_padL0 = real(ifft2(Gkpe_padL0 .* fft2(psie_padL0)));
psix_padLl = real(ifft2(Gkpx_padLl .* fft2(psix_padLl)));
psis_padl0 = real(ifft2(Gkps_padl0 .* fft2(psis_padl0)));

psie_flt = psie_padL0(1+padL0/2:end-padL0/2,:);
psix_flt = psix_padLl(1+padLl/2:end-padLl/2,:);
psis_flt = psis_padl0(1+padl0/2:end-padl0/2,:); 

% scale back to mean and std of raw noise amplitude fields
psie_flt = (psie_flt-mean(psie_flt(:))).*std(psie(:))./(std(psie_flt(:)) + eps) + mean(psie(:));
psix_flt = (psix_flt-mean(psix_flt(:))).*std(psix(:))./(std(psix_flt(:)) + eps) + mean(psix(:));
psis_flt = (psis_flt-mean(psis_flt(:))).*std(psis(:))./(std(psis_flt(:)) + eps) + mean(psis(:));

% taper noise potential fields towards closed boundaries and zero xtals
xtaper    = (1-exp(-chi./1e-4));
psix_flt  = psix_flt .* xtaper;
psis_flt  = psis_flt .* xtaper;


% take derivative of noise stream function and potentials to get flux components

% eddy noise components
psiec = (psie_flt(icz(1:end-1),icx(1:end-1))+psie_flt(icz(1:end-1),icx(2:end)) ...
      +  psie_flt(icz(2:end  ),icx(1:end-1))+psie_flt(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% particle eddy noise components
xixu = -ddx(psix_flt(icz,icx),1);
xixw = -ddz(psix_flt(icz,icx),1);

% particle settling noise components
xisu = -ddx(psis_flt(icz,icx),1);
xisw = -ddz(psis_flt(icz,icx),1);

% combine particle noise
xiwx = (xisw + xixw);
xiux = (xisu + xixu);

% get corresponding melt flux
xiwm = -xw(:,icx)./mw (:,icx).*xiwx;
xium = -xu(icz,:)./muu(icz,:).*xiux;



