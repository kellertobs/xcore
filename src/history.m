%*****  RECORD HISTORY OF DIAGNOSTICS  ************************************

[dsumBdtoo, dsumBdto] = deal(dsumBdto, dsumBdt);
[dsumMdtoo, dsumMdto] = deal(dsumMdto, dsumMdt);
[dsumXdtoo, dsumXdto] = deal(dsumXdto, dsumXdt);

stp = max(1,step/nrh);

% record model time
HST.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
HST.sumB(stp  ) = sum(rho(:)*h*h*1)+eps;  % [kg]
HST.sumM(stp  ) = sum(  M(:)*h*h*1)+eps;  % [kg]
HST.sumX(stp  ) = sum(  X(:)*h*h*1)+eps;  % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumBdt =  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1);  % [kg/s]
dsumMdt =  sum(sum(-Gm*h*h*1)) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1);  % [kg/s]
dsumXdt =  sum(sum( Gx*h*h*1)) ...
        +  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1);  % [kg/s]

if step>=2; HST.DB(stp) = (a2*HST.DB(max(1,stp-1)) + a3*HST.DB(max(1,stp-2)) + (b1*dsumBdt + b2*dsumBdto + b3*dsumBdtoo)*dt)/a1; else; HST.DB(stp) = 0; end  % [kg]
if step>=2; HST.DM(stp) = (a2*HST.DM(max(1,stp-1)) + a3*HST.DM(max(1,stp-2)) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; HST.DM(stp) = 0; end  % [kg]
if step>=2; HST.DX(stp) = (a2*HST.DX(max(1,stp-1)) + a3*HST.DX(max(1,stp-2)) + (b1*dsumXdt + b2*dsumXdto + b3*dsumXdtoo)*dt)/a1; else; HST.DX(stp) = 0; end  % [kg]

% record conservation error of mass M, heat S, components C
HST.EB(stp  ) = (HST.sumB(stp) - HST.DB(stp) - HST.sumB(1))./HST.sumB(1);  % [kg/kg]
HST.EM(stp  ) = (HST.sumM(stp) - HST.DM(stp) - HST.sumM(1))./HST.sumB(1);  % [kg/kg]
HST.EX(stp  ) = (HST.sumX(stp) - HST.DX(stp) - HST.sumX(1))./HST.sumB(1);  % [kg/kg]

% record variable and coefficient diagnostics
HST.V(stp,1) = min(V(:));
HST.V(stp,2) = sqrt(mean(V(:).^2));
HST.V(stp,3) = max(V(:));

HST.vx(stp,1) = min(vx(:));
HST.vx(stp,2) = sqrt(mean(vx(:).^2));
HST.vx(stp,3) = max(vx(:));

HST.xie(stp,1) = min(xie(:));
HST.xie(stp,2) = sqrt(mean(xie(:).^2));
HST.xie(stp,3) = max(xie(:));

HST.xix(stp,1) = min(xix(:));
HST.xix(stp,2) = sqrt(mean(xix(:).^2));
HST.xix(stp,3) = max(xix(:));

HST.W(stp,1) = min(-W(:));
HST.W(stp,2) = mean(abs(W(:)));
HST.W(stp,3) = max(-W(:));

HST.U(stp,1) = min(U(:));
HST.U(stp,2) = mean(abs(U(:)));
HST.U(stp,3) = max(U(:));

HST.P(stp,1) = min(P(:));
HST.P(stp,2) = mean(P(:));
HST.P(stp,3) = max(P(:));

HST.Pt(stp,1) = min(Pt(:));
HST.Pt(stp,2) = mean(abs(Pt(:)));
HST.Pt(stp,3) = max(Pt(:));

HST.x(stp,1) = min(x(:));
HST.x(stp,2) = mean(x(:));
HST.x(stp,3) = max(x(:));

HST.m(stp,1) = min(m(:));
HST.m(stp,2) = mean(m(:));
HST.m(stp,3) = max(m(:));

HST.chi(stp,1) = min(chi(:));
HST.chi(stp,2) = mean(chi(:));
HST.chi(stp,3) = max(chi(:));

HST.mu(stp,1) = min(mu(:));
HST.mu(stp,2) = mean(mu(:));
HST.mu(stp,3) = max(mu(:));

HST.Gx(stp,1) = min(Gx(:));
HST.Gx(stp,2) = mean(Gx(:));
HST.Gx(stp,3) = max(Gx(:));

HST.rho(stp,1) = min(rho(:));
HST.rho(stp,2) = mean(rho(:));
HST.rho(stp,3) = max(rho(:));

HST.wx(stp,1) = min(abs(wx(:)));
HST.wx(stp,2) = mean(abs(wx(:)));
HST.wx(stp,3) = max(abs(wx(:)));

HST.wm(stp,1) = min(abs(wm(:)));
HST.wm(stp,2) = mean(abs(wm(:)));
HST.wm(stp,3) = max(abs(wm(:)));

HST.RaD(stp,1) = min(RaD(:));
HST.RaD(stp,2) = geomean(RaD(:));
HST.RaD(stp,3) = max(RaD(:));

HST.Rs(stp,1) = min(Rs(:));
HST.Rs(stp,2) = geomean(Rs(:));
HST.Rs(stp,3) = max(Rs(:));

HST.ReD(stp,1) = min(ReD(:));
HST.ReD(stp,2) = geomean(ReD(:));
HST.ReD(stp,3) = max(ReD(:));

HST.Red(stp,1) = min(Red(:));
HST.Red(stp,2) = geomean(Red(:));
HST.Red(stp,3) = max(Red(:));

HST.ks(stp,1) = min(ks(:));
HST.ks(stp,2) = geomean(ks(:));
HST.ks(stp,3) = max(ks(:));

HST.ke(stp,1) = min(ke(:));
HST.ke(stp,2) = geomean(ke(:));
HST.ke(stp,3) = max(ke(:));

HST.eta0(stp,1) = min(etamix(:));
HST.eta0(stp,2) = geomean(etamix(:));
HST.eta0(stp,3) = max(etamix(:));

HST.etax(stp,1) = min(etax(:));
HST.etax(stp,2) = geomean(etax(:));
HST.etax(stp,3) = max(etax(:));

HST.eta(stp,1) = min(eta(:));
HST.eta(stp,2) = geomean(eta(:));
HST.eta(stp,3) = max(eta(:));

HST.etae(stp,1) = min(etae(:));
HST.etae(stp,2) = geomean(etae(:));
HST.etae(stp,3) = max(etae(:));

HST.etas(stp,1) = min(etas(:));
HST.etas(stp,2) = geomean(etas(:));
HST.etas(stp,3) = max(etas(:));

HST.etat(stp,1) = min(etat(:));
HST.etat(stp,2) = geomean(etat(:));
HST.etat(stp,3) = max(etat(:));

% time-averaged diagnostics
stp0 = round((1./stp.^10 + 1./(30 + (stp-30)/2).^10).^-(1/10));
HST.x_tavg(stp,:)    = mean(HST.x(stp0:stp,:),1);

HST.V_tavg(stp,:)    = mean(HST.V(stp0:stp,:),1);
HST.vx_tavg(stp,:)   = mean(HST.vx(stp0:stp,:),1);

HST.xie_tavg(stp,:)  = mean(HST.xie(stp0:stp,:),1);
HST.xix_tavg(stp,:)  = mean(HST.xix(stp0:stp,:),1);

HST.Rs_tavg(stp,:)   = mean(HST.Rs(stp0:stp,:),1);
HST.RaD_tavg(stp,:)  = mean(HST.RaD(stp0:stp,:),1);

HST.ReD_tavg(stp,:)  = mean(HST.ReD(stp0:stp,:),1);
HST.Red_tavg(stp,:)  = mean(HST.Red(stp0:stp,:),1);

if stp>1
    HST.Dx_tavg(stp,:) = abs(mean(diff(HST.x(min(stp-1,stp0):stp,2))));
else
    HST.Dx_tavg(stp,:) = 1e-3;
end

if stp>1
    HST.DV_tavg(stp,:) = abs(mean(diff(HST.V(min(stp-1,stp0):stp,3))));
else
    HST.DV_tavg(stp,:) = 1e-3;
end
