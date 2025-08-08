% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter==1

    % Generate new white noise
    rsw  = randn(Nz+1, Nx+0);
    rsu  = randn(Nz+0, Nx+1);
    rexw = randn(Nz+1, Nx+0);
    rexu = randn(Nz+0, Nx+1);
    pse  = randn(Nz+1, Nx+1);

    rsw  = fft2(rsw);
    rsu  = fft2(rsu);
    rexw = fft2(rexw);
    rexu = fft2(rexu);
    pse  = fft2(pse);

    % Filter white noise spatially to decorrelation length
    rsw  = real(ifft2(Gkws .* rsw));
    rsu  = real(ifft2(Gkus .* rsu));
    rexw = real(ifft2(Gkwe .* rexw));
    rexu = real(ifft2(Gkue .* rexu));
    pse  = real(ifft2(Gkps .* pse));

    reu  =  ddz(pse,1);
    rew  = -ddx(pse,1);

    % rescale to unit RMS speed
    ss   = sqrt(mean(rsw (:).^2) + mean(rsu (:).^2));
    sex  = sqrt(mean(rexw(:).^2) + mean(rexu(:).^2));
    se   = sqrt(mean(rew (:).^2) + mean(reu (:).^2));
    rsw  = (rsw  - mean(rsw (:))) / ss;
    rsu  = (rsu  - mean(rsu (:))) / ss;
    rexw = (rexw - mean(rexw(:))) / sex;
    rexu = (rexu - mean(rexu(:))) / sex;
    rew  = (rew  - mean(rew (:))) / se;
    reu  = (reu  - mean(reu (:))) / se;

end

% noise decorrelation time
taue = elle/2./(V +eps);
taus = ells/2./(vx+eps);

% temporal evolution factor
Fe   = exp(-dt./taue);
Fs   = exp(-dt./taus);

Few = (Fe(icz(1:end-1),icx)+Fe(icz(2:end),icx))/2;
Feu = (Fe(icz,icx(1:end-1))+Fe(icz,icx(2:end)))/2;

Fsw = (Fs(icz(1:end-1),icx)+Fs(icz(2:end),icx))/2;
Fsu = (Fs(icz,icx(1:end-1))+Fs(icz,icx(2:end)))/2;

% random noise source amplitude
bndtapere = (1 - exp((-ZZ+h/2)/max(h,elle/2)) - exp(-(D-ZZ-h/2)/max(h,elle/2)).*(1-open));
bndtapers = (1 - exp((-ZZ+h/2)/max(h,ells/2)) - exp(-(D-ZZ-h/2)/max(h,ells/2)).*(1-open));

sgs   = Xi * sqrt(chi.*ks./taus .* (ells./(ells+h)).^3) .* bndtapers; % segregation noise speed
sge   = Xi * sqrt(     ke./taue .* (elle./(elle+h)).^3) .* bndtapere; % eddy mixture noise speed
sgex  = Xi * sqrt(chi.*ke./taue .* (elle./(elle+h)).^3) .* bndtapere; % eddy crystal noise speed

sgsw  = (sgs(icz(1:end-1),icx) + sgs(icz(2:end),icx))./2;
sgsu  = (sgs(icz,icx(1:end-1)) + sgs(icz,icx(2:end)))./2;

sgexw = (sgex(icz(1:end-1),icx) + sgex(icz(2:end),icx))./2;
sgexu = (sgex(icz,icx(1:end-1)) + sgex(icz,icx(2:end)))./2;

sgew  = (sge(icz(1:end-1),icx) + sge(icz(2:end),icx))./2;
sgeu  = (sge(icz,icx(1:end-1)) + sge(icz,icx(2:end)))./2;

% Ornstein–Uhlenbeck time update
xisw  =  Fsw .* xiswo  + sqrt(1 - Fsw.^2) .* sgsw  .* rsw(:,icx);
xisu  =  Fsu .* xisuo  + sqrt(1 - Fsu.^2) .* sgsu  .* rsu(icz,:);
xiexw =  Few .* xiexwo + sqrt(1 - Few.^2) .* sgexw .* rexw(:,icx);
xiexu =  Feu .* xiexuo + sqrt(1 - Feu.^2) .* sgexu .* rexu(icz,:);
xiew  =  Few .* xiewo  + sqrt(1 - Few.^2) .* sgew  .* rew(:,icx);
xieu  =  Feu .* xieuo  + sqrt(1 - Feu.^2) .* sgeu  .* reu(icz,:);

% update phase noise speeds
xiwx  = xisw + xiexw;
xiux  = xisu + xiexu;
xiwm  = -xw(:,icx)./mw (:,icx).*xiwx;
xium  = -xu(icz,:)./muu(icz,:).*xiux;