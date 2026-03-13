% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter<=1
    % generate white noise
    re = randn(Nz+0, Nx+0);
    % rx = randn(Nz+0, Nx+0);
    rs = randn(Nz+0, Nx+0);

    % % filter white noise spatially to decorrelation length
    % re = real(ifft2(Gkpe .* fft2(re)));
    % rx = real(ifft2(Gkpe .* fft2(rx)));
    % rs = real(ifft2(Gkps .* fft2(rs)));
    % 
    % % rescale filtered noise to zero mean and unit variance
    % re = (re - mean(re(:))) ./ std(re(:));
    % rx = (rx - mean(rx(:))) ./ std(rx(:));
    % rs = (rs - mean(rs(:))) ./ std(rs(:));
end


% noise decorrelation time
taue = L0h./(V +eps);
taus = l0h./(vx+eps);

% noise flux amplitudes
sge = Xi * sqrt(          fReL.*ke     ./taue .* (L0h./(L0h+h)).^3);    % eddy mixture noise speed
sgx = Xi * sqrt(chi.*sqrt(fReL.*ke.*ks)./taue .* (L0h./(L0h+h)).^3);    % eddy crystal noise speed
sgs = Xi * sqrt(chi.*               ks ./taus .* (l0h./(l0h+h)).^3);    % settling noise speed


% Ornstein–Uhlenbeck time update for evolving noise
Fte  = exp(-dt./taue);                                        % eddy noise time evolution factor
Fts  = exp(-dt./taus);                                        % settling noise time evolution factor

fL   = 2/sqrt(1 - exp(-h^2/(2*L0h^2)));                       % scaling factor for potential field to noise component variance
fl   = 2/sqrt(1 - exp(-h^2/(2*l0h^2)));                       % scaling factor for potential field to noise component variance
 
psie = Fte .* psieo + sqrt(1 - Fte.^2) .* sge.*fL .* re;      % update eddy noise stream function
psix = Fte .* psixo + sqrt(1 - Fte.^2) .* sgx.*fL .* re;      % update particle eddy noise potential
psis = Fts .* psiso + sqrt(1 - Fts.^2) .* sgs.*fl .* rs;      % update particle settling noise potential


% filter noise amplitudes spatially to decorrelation length
% pad top/bot to avoid periodic bleed over on non-periodic boundaries
psie_padL0 = zeros(Nz+padL0,Nx); 
psix_padL0 = zeros(Nz+padL0,Nx);
psis_padl0 = zeros(Nz+padl0,Nx);

psie_padL0(1+padL0/2:end-padL0/2,:) = psie;
% psie_padL0(            1:padL0/2,:) = repmat(psie(1  ,:),padL0/2,1).*(1-open_cnv);
% psie_padL0(end-padL0/2+1:end    ,:) = repmat(psie(end,:),padL0/2,1).*(1-open_cnv);

psix_padL0(1+padL0/2:end-padL0/2,:) = psix;
% psix_padL0(            1:padL0/2,:) = repmat(psix(1  ,:),padL0/2,1).*(1-open_cnv);
% psix_padL0(end-padL0/2+1:end    ,:) = repmat(psix(end,:),padL0/2,1).*(1-open_cnv);

psis_padl0(1+padl0/2:end-padl0/2,:) = psis;
% psis_padl0(            1:padl0/2,:) = repmat(psis(1  ,:),padl0/2,1).*(1-open_sgr);
% psis_padl0(end-padl0/2+1:end    ,:) = repmat(psis(end,:),padl0/2,1).*(1-open_sgr);

psie_padL0 = real(ifft2(Gkpe_padL0 .* fft2(psie_padL0)));
psix_padL0 = real(ifft2(Gkpe_padL0 .* fft2(psix_padL0)));
psis_padl0 = real(ifft2(Gkps_padl0 .* fft2(psis_padl0)));

psie_flt = psie_padL0(1+padL0/2:end-padL0/2,:);
psix_flt = psix_padL0(1+padL0/2:end-padL0/2,:);
psis_flt = psis_padl0(1+padl0/2:end-padl0/2,:); 

% scale back to mean and std of raw noise amplitude fields
psie_flt = (psie_flt-mean(psie_flt(:))).*std(psie(:))./(std(psie_flt(:)) + eps) + mean(psie(:));
psix_flt = (psix_flt-mean(psix_flt(:))).*std(psix(:))./(std(psix_flt(:)) + eps) + mean(psix(:));
psis_flt = (psis_flt-mean(psis_flt(:))).*std(psis(:))./(std(psis_flt(:)) + eps) + mean(psis(:));

% taper noise potential fields towards closed boundaries and zero xtals
xtaper    = (1-exp(-chi./1e-4));
% psie_flt  = psie_flt  ;
psix_flt  = psix_flt   .* xtaper;
psis_flt  = psis_flt   .* xtaper;


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



