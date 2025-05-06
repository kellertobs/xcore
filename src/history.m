% record run history

[dsumBdtoo, dsumBdto] = deal(dsumBdto, dsumBdt);
[dsumMdtoo, dsumMdto] = deal(dsumMdto, dsumMdt);
[dsumXdtoo, dsumXdto] = deal(dsumXdto, dsumXdt);

stp = max(1,step/nrh);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumB(stp  ) = sum(rho(:)*h*h*1)+eps;  % [kg]
hist.sumM(stp  ) = sum(  M(:)*h*h*1)+eps;  % [kg]
hist.sumX(stp  ) = sum(  X(:)*h*h*1)+eps;  % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumBdt =  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1);  % [kg/s]
dsumMdt =  sum(sum(Gm*h*h*1)) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1);  % [kg/s]
dsumXdt =  sum(sum(Gx*h*h*1)) ...
        +  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1);  % [kg/s]

if step>=2; hist.DB(stp) = (a2*hist.DB(max(1,stp-1)) + a3*hist.DB(max(1,stp-2)) + (b1*dsumBdt + b2*dsumBdto + b3*dsumBdtoo)*dt)/a1; else; hist.DB(stp) = 0; end  % [kg]
if step>=2; hist.DM(stp) = (a2*hist.DM(max(1,stp-1)) + a3*hist.DM(max(1,stp-2)) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; hist.DM(stp) = 0; end  % [kg]
if step>=2; hist.DX(stp) = (a2*hist.DX(max(1,stp-1)) + a3*hist.DX(max(1,stp-2)) + (b1*dsumXdt + b2*dsumXdto + b3*dsumXdtoo)*dt)/a1; else; hist.DX(stp) = 0; end  % [kg]

% record conservation error of mass M, heat S, components C
hist.EB(stp  ) = (hist.sumB(stp) - hist.DB(stp) - hist.sumB(1))./hist.sumB(1);  % [kg/kg]
hist.EM(stp  ) = (hist.sumM(stp) - hist.DM(stp) - hist.sumM(1))./hist.sumB(1);  % [kg/kg]
hist.EX(stp  ) = (hist.sumX(stp) - hist.DX(stp) - hist.sumX(1))./hist.sumB(1);  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(stp,1) = min(min(-W(:,2:end-1)));
hist.W(stp,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(stp,3) = max(max(-W(:,2:end-1)));

hist.U(stp,1) = min(min(U(2:end-1,:)));
hist.U(stp,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(stp,3) = max(max(U(2:end-1,:)));

hist.P(stp,1) = min(min(P(2:end-1,2:end-1)));
hist.P(stp,2) = mean(mean(P(2:end-1,2:end-1)));
hist.P(stp,3) = max(max(P(2:end-1,2:end-1)));

hist.Pt(stp,1) = min(min(Pt));
hist.Pt(stp,2) = mean(mean(abs(Pt)));
hist.Pt(stp,3) = max(max(Pt));

hist.x(stp,1) = min(min(x));
hist.x(stp,2) = mean(mean(x));
hist.x(stp,3) = max(max(x));

hist.m(stp,1) = min(min(m));
hist.m(stp,2) = mean(mean(m));
hist.m(stp,3) = max(max(m));

hist.chi(stp,1) = min(min(chi));
hist.chi(stp,2) = mean(mean(chi));
hist.chi(stp,3) = max(max(chi));

hist.mu(stp,1) = min(min(mu));
hist.mu(stp,2) = mean(mean(mu));
hist.mu(stp,3) = max(max(mu));

hist.Gm(stp,1) = min(min(Gm));
hist.Gm(stp,2) = mean(mean(Gm));
hist.Gm(stp,3) = max(max(Gm));

hist.Gx(stp,1) = min(min(Gx));
hist.Gx(stp,2) = mean(mean(Gx));
hist.Gx(stp,3) = max(max(Gx));

hist.dV(stp,1) = min(min(dV));
hist.dV(stp,2) = mean(mean(dV));
hist.dV(stp,3) = max(max(dV));

hist.rho(stp,1) = min(min(rho));
hist.rho(stp,2) = mean(mean(rho));
hist.rho(stp,3) = max(max(rho));

hist.eta(stp,1) = min(min(eta0));
hist.eta(stp,2) = geomean(geomean(eta0));
hist.eta(stp,3) = max(max(eta0));

hist.eta(stp,1) = min(min(eta));
hist.eta(stp,2) = geomean(geomean(eta));
hist.eta(stp,3) = max(max(eta));

hist.wx(stp,1) = min(min(abs(wx(:,2:end-1))));
hist.wx(stp,2) = mean(mean(abs(wx(:,2:end-1))));
hist.wx(stp,3) = max(max(abs(wx(:,2:end-1))));

hist.wx(stp,1) = min(min(abs(wm(:,2:end-1))));
hist.wx(stp,2) = mean(mean(abs(wm(:,2:end-1))));
hist.wx(stp,3) = max(max(abs(wm(:,2:end-1))));

hist.Ra(stp,1) = min(min(Ra));
hist.Ra(stp,2) = geomean(geomean(Ra));
hist.Ra(stp,3) = max(max(Ra));

hist.RaD(stp,1) = min(min(RaD));
hist.RaD(stp,2) = geomean(geomean(RaD));
hist.RaD(stp,3) = max(max(RaD));

hist.Rux(stp,1) = min(min(Rux));
hist.Rux(stp,2) = geomean(geomean(Rux));
hist.Rux(stp,3) = max(max(Rux));

hist.Re(stp,1) = min(min(Re));
hist.Re(stp,2) = geomean(geomean(Re));
hist.Re(stp,3) = max(max(Re));

hist.ReD(stp,1) = min(min(ReD));
hist.ReD(stp,2) = geomean(geomean(ReD));
hist.ReD(stp,3) = max(max(ReD));

hist.Rex(stp,1) = min(min(Rex));
hist.Rex(stp,2) = geomean(geomean(Rex));
hist.Rex(stp,3) = max(max(Rex));

hist.kwx(stp,1) = min(min(kwx));
hist.kwx(stp,2) = geomean(geomean(kwx));
hist.kwx(stp,3) = max(max(kwx));

hist.ke(stp,1) = min(min(ke));
hist.ke(stp,2) = geomean(geomean(ke));
hist.ke(stp,3) = max(max(ke));