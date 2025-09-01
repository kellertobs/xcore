% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter==1
    % generate white noise
    re = randn(Nz+0, Nx+0);
    rx = randn(Nz+0, Nx+0);
    rs = randn(Nz+0, Nx+0);

    % filter white noise spatially to decorrelation length
    re = real(ifft2(Gkpe .* fft2(re)));
    rx = real(ifft2(Gkpe .* fft2(rx)));
    rs = real(ifft2(Gkps .* fft2(rs)));

    % rescale filtered noise to zero mean and unit variance
    re = (re - mean(re(:))) ./ std(re(:));
    rx = (rx - mean(rx(:))) ./ std(rx(:));
    rs = (rs - mean(rs(:))) ./ std(rs(:));
end


% noise decorrelation time
taue = L0h./(V +eps);
taus = l0h./(vx+eps);

% noise flux amplitudes
sge_raw = Xi * sqrt(     fReL.*ke./taue .* (L0h./(L0h+h)).^3);    % eddy mixture noise speed
sgx_raw = Xi * sqrt(chi.*fReL.*ke./taue .* (L0h./(L0h+h)).^3);    % eddy crystal noise speed
sgs_raw = Xi * sqrt(chi.*      ks./taus .* (l0h./(l0h+h)).^3);    % segregation noise speed


% filter noise amplitudes spatially to decorrelation length
% pad top/bot to avoid periodic bleed over on non-periodic boundaries
sge_padL0 = zeros(Nz+padL0,Nx); 
sgx_padL0 = zeros(Nz+padL0,Nx);
sgs_padl0 = zeros(Nz+padl0,Nx);

sge_padL0(1+padL0/2:end-padL0/2,:) = sge_raw;
sge_padL0(            1:padL0/2,:) = repmat(sge_raw(1  ,:),padL0/2,1);
sge_padL0(end-padL0/2+1:end    ,:) = repmat(sge_raw(end,:),padL0/2,1);

sgx_padL0(1+padL0/2:end-padL0/2,:) = sgx_raw;
sgx_padL0(            1:padL0/2,:) = repmat(sgx_raw(1  ,:),padL0/2,1);
sgx_padL0(end-padL0/2+1:end    ,:) = repmat(sgx_raw(end,:),padL0/2,1);

sgs_padl0(1+padl0/2:end-padl0/2,:) = sgs_raw;
sgs_padl0(            1:padl0/2,:) = repmat(sgs_raw(1  ,:),padl0/2,1);
sgs_padl0(end-padl0/2+1:end    ,:) = repmat(sgs_raw(end,:),padl0/2,1);

sge_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sge_padL0)));
sgx_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sgx_padL0)));
sgs_padl0 = real(ifft2(Gkps_padl0 .* fft2(sgs_padl0)));

sge = sge_padL0(1+padL0/2:end-padL0/2,:);
sgx = sgx_padL0(1+padL0/2:end-padL0/2,:);
sgs = sgs_padl0(1+padl0/2:end-padl0/2,:); 

% scale back to mean and std of raw noise amplitude fields
sge = (sge-mean(sge(:))).*std(sge_raw(:))./std(sge(:)) + mean(sge_raw(:));
sgx = (sgx-mean(sgx(:))).*std(sgx_raw(:))./std(sgx(:)) + mean(sgx_raw(:));
sgs = (sgs-mean(sgs(:))).*std(sgs_raw(:))./std(sgs(:)) + mean(sgs_raw(:));


% Ornstein–Uhlenbeck time update for evolving noise
Fte  = exp(-dt./taue);                                        % eddy noise time evolution factor
Fts  = exp(-dt./taus);                                        % settling noise time evolution factor

fL   = 2/sqrt(1 - exp(-h^2/(2*L0h^2)));                       % scaling factor for potential field to noise component variance
fl   = 2/sqrt(1 - exp(-h^2/(2*l0h^2)));                       % scaling factor for potential field to noise component variance
 
psie = Fte .* psieo + sqrt(1 - Fte.^2) .* sge.*fL .* re;      % update eddy noise stream function
psix = Fte .* psixo + sqrt(1 - Fte.^2) .* sgx.*fL .* rx;      % update particle eddy noise potential
psis = Fts .* psiso + sqrt(1 - Fts.^2) .* sgs.*fl .* rs;      % update particle settling noise potential


% take derivative of noise stream function and potentials to get flux components

% eddy noise components
psiec = (psie(icz(1:end-1),icx(1:end-1))+psie(icz(1:end-1),icx(2:end)) ...
      +  psie(icz(2:end  ),icx(1:end-1))+psie(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% particle eddy noise components
xixu = -ddx(psix(icz,icx),1);
xixw = -ddz(psix(icz,icx),1);

% particle settling noise components
xisu = -ddx(psis (icz,icx),1);
xisw = -ddz(psis (icz,icx),1);

% combine particle noise
xiwx = (xisw + xixw);
xiux = (xisu + xixu);

% get corresponding melt flux
xiwm = -xw(:,icx)./mw (:,icx).*xiwx;
xium = -xu(icz,:)./muu(icz,:).*xiux;



