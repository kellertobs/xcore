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
    rs   = sqrt(((rsw(1:end-1,:)+rsw(2:end,:))/2).^2 ...
              + ((rsu(:,1:end-1)+rsu(:,2:end))/2).^2);
    re   = sqrt(((rew(1:end-1,:)+rew(2:end,:))/2).^2 ...
              + ((reu(:,1:end-1)+reu(:,2:end))/2).^2);
    rex  = sqrt(((rexw(1:end-1,:)+rexw(2:end,:))/2).^2 ...
              + ((rexu(:,1:end-1)+rexu(:,2:end))/2).^2);
    rsw  = (rsw  - mean(rsw (:))) / geomean(rs(:));
    rsu  = (rsu  - mean(rsu (:))) / geomean(rs(:));
    rexw = (rexw - mean(rexw(:))) / geomean(rex(:));
    rexu = (rexu - mean(rexu(:))) / geomean(rex(:));
    rew  = (rew  - mean(rew (:))) / geomean(re(:));
    reu  = (reu  - mean(reu (:))) / geomean(re(:));

end

% noise decorrelation time
taue = L0/2./(V +eps);
taus = l0/2./(vx+eps);

% temporal evolution factor
Fe   = exp(-dt./taue);
Fs   = exp(-dt./taus);

Few = (Fe(icz(1:end-1),icx)+Fe(icz(2:end),icx))/2;
Feu = (Fe(icz,icx(1:end-1))+Fe(icz,icx(2:end)))/2;

Fsw = (Fs(icz(1:end-1),icx)+Fs(icz(2:end),icx))/2;
Fsu = (Fs(icz,icx(1:end-1))+Fs(icz,icx(2:end)))/2;

% random noise source amplitude
sgs   = Xi * sqrt(chi.*      ks./taus .* (l0./(l0+h)).^3);  % segregation noise speed
sgex  = Xi * sqrt(chi.*fReL.*ke./taue .* (L0./(L0+h)).^3);  % eddy crystal noise speed
sge   = Xi * sqrt(     fReL.*ke./taue .* (L0./(L0+h)).^3);  % eddy mixture noise speed

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

% enforce closed boundary conditions
if ~open_sgr; xisw ([1 end],:) = 0; end
if ~open_cnv; xiew ([1 end],:) = 0; else; xiew (1,:) = 0; end
if ~open_cnv; xiexw([1 end],:) = 0; else; xiexw(1,:) = 0; end

% update phase noise speeds
xiwx  = xisw + xiexw;
xiux  = xisu + xiexu;
xiwm  = -xw(:,icx)./mw (:,icx).*xiwx;
xium  = -xu(icz,:)./muu(icz,:).*xiux;


