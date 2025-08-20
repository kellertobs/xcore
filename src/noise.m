% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter==1
    % generate white noise
    re  = randn(Nz+0, Nx+0);
    rex = randn(Nz+0, Nx+0);
    rs  = randn(Nz+0, Nx+0);

    % filter white noise spatially to decorrelation length
    re  = real(ifft2(Gkpe .* fft2(re )));
    rex = real(ifft2(Gkpe .* fft2(rex)));
    rs  = real(ifft2(Gkps .* fft2(rs )));

    % rescale filtered noise to zero mean and unit variance
    re  = (re  - mean(re (:))) ./ std(re (:));
    rex = (rex - mean(rex(:))) ./ std(rex(:));
    rs  = (rs  - mean(rs (:))) ./ std(rs (:));
end

% noise decorrelation time
taue = L0h./(V +eps);
taus = l0h./(vx+eps);

% temporal evolution factor
Fe   = exp(-dt./taue);
Fs   = exp(-dt./taus);

% noise flux amplitudes
sgs   = Xi * sqrt(chi.*      ks./taus .* (l0h./(l0h+h)).^3);  % segregation noise speed
sgex  = Xi * sqrt(chi.*fReL.*ke./taue .* (L0h./(L0h+h)).^3);  % eddy crystal noise speed
sge   = Xi * sqrt(     fReL.*ke./taue .* (L0h./(L0h+h)).^3);  % eddy mixture noise speed

% Ornstein–Uhlenbeck time update
fL   =  h./sqrt(2 * (1 - exp(-(h./(L0h+h)).^2)));  % scaling factor for potential field to noise component variance
fl   =  h./sqrt(2 * (1 - exp(-(h./(l0h+h)).^2)));  % scaling factor for potential field to noise component variance

psie  =  Fe .* psieo  + sqrt(1 - Fe .^2) .* sge .*fL .* re;    % update eddy noise stream function
psiex =  Fe .* psiexo + sqrt(1 - Fe .^2) .* sgex.*fL .* rex;   % update particle eddy noise potential
psis  =  Fs .* psiso  + sqrt(1 - Fs .^2) .* sgs .*fl .* rs;    % update particle settling noise potential

% take derivative of noise stream function and potentials to get flux components

% eddy noise components
psiec  = (psie(icz(1:end-1),icx(1:end-1))+psie(icz(1:end-1),icx(2:end)) ...
       +  psie(icz(2:end  ),icx(1:end-1))+psie(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,h); xieu = xieu(icz,:);
xiew  = -ddx(psiec,h); xiew = xiew(:,icx);

% particle eddy noise components
xiexu =  ddx(psiex(icz,icx),h);
xiexw =  ddz(psiex(icz,icx),h);

% particle settling noise components
xisu  =  ddx(psis (icz,icx),h);
xisw  =  ddz(psis (icz,icx),h);

% combine particle noise and get corresponding melt flux
xiwx  = (xisw + xiexw);
xiux  = (xisu + xiexu);
xiwm  = -xw(:,icx)./mw (:,icx).*xiwx;
xium  = -xu(icz,:)./muu(icz,:).*xiux;



