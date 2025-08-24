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

% noise flux amplitudes
sge  = Xi * sqrt(     fReL.*ke./taue .* (L0h./(L0h+h)).^3);    % eddy mixture noise speed
sgex = Xi * sqrt(chi.*fReL.*ke./taue .* (L0h./(L0h+h)).^3);    % eddy crystal noise speed
sgs  = Xi * sqrt(chi.*      ks./taus .* (l0h./(l0h+h)).^3);    % segregation noise speed

% Ornstein–Uhlenbeck time update for evolving noise
Fte   = exp(-dt./taue);                                        % eddy noise time evolution factor
Fts   = exp(-dt./taus);                                        % settling noise time evolution factor

fL    =  2/sqrt(1 - exp(-h^2/(2*L0h^2)));                      % scaling factor for potential field to noise component variance
fl    =  2/sqrt(1 - exp(-h^2/(2*l0h^2)));                      % scaling factor for potential field to noise component variance

psie  =  Fte .* psieo  + sqrt(1 - Fte.^2) .* sge .*fL .* re;   % update eddy noise stream function
psiex =  Fte .* psiexo + sqrt(1 - Fte.^2) .* sgex.*fL .* rex;  % update particle eddy noise potential
psis  =  Fts .* psiso  + sqrt(1 - Fts.^2) .* sgs .*fl .* rs;   % update particle settling noise potential


% take derivative of noise stream function and potentials to get flux components

% eddy noise components
psiec  = (psie(icz(1:end-1),icx(1:end-1))+psie(icz(1:end-1),icx(2:end)) ...
       +  psie(icz(2:end  ),icx(1:end-1))+psie(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% particle eddy noise components
xiexu = -ddx(psiex(icz,icx),1);
xiexw = -ddz(psiex(icz,icx),1);

% particle settling noise components
xisu  = -ddx(psis (icz,icx),1);
xisw  = -ddz(psis (icz,icx),1);

% combine particle noise
xiwx  = (xisw + xiexw);
xiux  = (xisu + xiexu);

% get corresponding melt flux
xiwm  = -xw(:,icx)./mw (:,icx).*xiwx;
xium  = -xu(icz,:)./muu(icz,:).*xiux;



