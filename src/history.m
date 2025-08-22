%*****  RECORD HISTORY OF DIAGNOSTICS  ************************************

[dsumBdtoo, dsumBdto] = deal(dsumBdto, dsumBdt);
[dsumMdtoo, dsumMdto] = deal(dsumMdto, dsumMdt);
[dsumXdtoo, dsumXdto] = deal(dsumXdto, dsumXdt);

stp = max(1,step/nrh);

% record model time
HST.time(stp) = time;
HST.dt(stp)   = dt;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
HST.sumB(stp  ) = sum(rho(:)*h*h*1)+eps;  % [kg]
HST.sumM(stp  ) = sum(  M(:)*h*h*1)+eps;  % [kg]
HST.sumX(stp  ) = sum(  X(:)*h*h*1)+eps;  % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumBdt =  sum(qz_advn_X(1,2:end-1)*h*1) - sum(qz_advn_X(end,2:end-1)*h*1) ...
        +  sum(qz_advn_M(1,2:end-1)*h*1) - sum(qz_advn_M(end,2:end-1)*h*1) ...
        +  sum(qz_dffn_X(1,2:end-1)*h*1) - sum(qz_dffn_X(end,2:end-1)*h*1) ...
        +  sum(qz_dffn_M(1,2:end-1)*h*1) - sum(qz_dffn_M(end,2:end-1)*h*1);  % [kg/s]
dsumMdt =  sum(sum(-Gx*h*h*1)) ...
        +  sum(qz_advn_M(1,2:end-1)*h*1) - sum(qz_advn_M(end,2:end-1)*h*1) ...
        +  sum(qz_dffn_M(1,2:end-1)*h*1) - sum(qz_dffn_M(end,2:end-1)*h*1);  % [kg/s]
dsumXdt =  sum(sum( Gx*h*h*1)) ...
        +  sum(qz_advn_X(1,2:end-1)*h*1) - sum(qz_advn_X(end,2:end-1)*h*1) ...
        +  sum(qz_dffn_X(1,2:end-1)*h*1) - sum(qz_dffn_X(end,2:end-1)*h*1);  % [kg/s]

if step>=2; HST.DB(stp) = (a2*HST.DB(max(1,stp-1)) + a3*HST.DB(max(1,stp-2)) + (b1*dsumBdt + b2*dsumBdto + b3*dsumBdtoo)*dt)/a1; else; HST.DB(stp) = 0; end  % [kg]
if step>=2; HST.DM(stp) = (a2*HST.DM(max(1,stp-1)) + a3*HST.DM(max(1,stp-2)) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; HST.DM(stp) = 0; end  % [kg]
if step>=2; HST.DX(stp) = (a2*HST.DX(max(1,stp-1)) + a3*HST.DX(max(1,stp-2)) + (b1*dsumXdt + b2*dsumXdto + b3*dsumXdtoo)*dt)/a1; else; HST.DX(stp) = 0; end  % [kg]

% record conservation error of mass M, heat S, components C
HST.EB(stp  ) = (HST.sumB(stp) - HST.DB(stp) - HST.sumB(1))./HST.sumB(1);  % [kg/kg]
HST.EM(stp  ) = (HST.sumM(stp) - HST.DM(stp) - HST.sumM(1))./HST.sumB(1);  % [kg/kg]
HST.EX(stp  ) = (HST.sumX(stp) - HST.DX(stp) - HST.sumX(1))./HST.sumB(1);  % [kg/kg]

% record variable and coefficient diagnostics
HST.x(stp,1) = min(x(:));
HST.x(stp,2) = mean(x(:));
HST.x(stp,3) = max(x(:));
HST.x(stp,4) = std(x(:));

HST.V(stp,1) = min(V(:));
HST.V(stp,2) = geomean(V(:));
HST.V(stp,3) = max(V(:));
HST.V(stp,4) = rms(V(:));

HST.vx(stp,1) = min(vx(:));
HST.vx(stp,2) = geomean(vx(:));
HST.vx(stp,3) = max(vx(:));
HST.vx(stp,4) = rms(vx(:));

HST.vm(stp,1) = min(vm(:));
HST.vm(stp,2) = geomean(vm(:));
HST.vm(stp,3) = max(vm(:));
HST.vm(stp,4) = rms(vm(:));

HST.xie(stp,1) = min(xie(:));
HST.xie(stp,2) = geomean(xie(:));
HST.xie(stp,3) = max(xie(:));
HST.xie(stp,4) = rms(xie(:));

HST.xix(stp,1) = min(xix(:));
HST.xix(stp,2) = geomean(xix(:));
HST.xix(stp,3) = max(xix(:));
HST.xix(stp,4) = rms(xix(:));

HST.rho(stp,1) = min(rho(:));
HST.rho(stp,2) = geomean(rho(:));
HST.rho(stp,3) = max(rho(:));

HST.Ra(stp,1) = min(Ra(:));
HST.Ra(stp,2) = geomean(Ra(:));
HST.Ra(stp,3) = max(Ra(:));

HST.Rc(stp,1) = min(Rc(:));
HST.Rc(stp,2) = geomean(Rc(:));
HST.Rc(stp,3) = max(Rc(:));

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

HST.kx(stp,1) = min(kx(:));
HST.kx(stp,2) = geomean(kx(:));
HST.kx(stp,3) = max(kx(:));

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
stp0 = round((1./stp.^10 + 1./(30 + (stp-30)*0.75).^10).^-(1/10));
HST.x_tavg(stp,:)    = mean(HST.x(stp0:stp,:),1);

HST.V_tavg(stp,:)    = mean(HST.V(stp0:stp,:),1);
HST.vx_tavg(stp,:)   = mean(HST.vx(stp0:stp,:),1);

HST.xie_tavg(stp,:)  = mean(HST.xie(stp0:stp,:),1);
HST.xix_tavg(stp,:)  = mean(HST.xix(stp0:stp,:),1);

HST.Rc_tavg(stp,:)   = mean(HST.Rc(stp0:stp,:),1);
HST.Ra_tavg(stp,:)   = mean(HST.Ra(stp0:stp,:),1);

HST.ReD_tavg(stp,:)  = mean(HST.ReD(stp0:stp,:),1);
HST.Red_tavg(stp,:)  = mean(HST.Red(stp0:stp,:),1);
