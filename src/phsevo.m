%*****  PHASE FRACTION EVOLUTION  *****************************************

tic;

%***  update phase fraction densities

% phase advection rates and fluxes
[advn_X,qz_advn_X,qx_advn_X] = advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],{(  [xeq,x0]).*rhom0,'periodic'});
[advn_M,qz_advn_M,qx_advn_M] = advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],{(1-[xeq,x0]).*rhom0,'periodic'});
advn_rho = advn_X+advn_M;

% phase diffusion rates and fluxes
[dffn_X,qz_dffn_X,qx_dffn_X] = diffus(chi,X.*kx,h,[1,2],BCD);
[dffn_M,qz_dffn_M,qx_dffn_M] = diffus( mu,X.*kx,h,[1,2],BCD);

% boundary phase change rate
if open_sgr
    tau_x = h./(W0 + w0) + dt;
else
    tau_x = (L0h+l0h)^2/(1e-6 + fReL.*mean(ke(1,:)));
end
Gx       = R.*max(0,xeq-x).*rho./tau_x.*bndshape;

% total rates of change
dXdt     = - advn_X + dffn_X + Gx;

% residual of phase density evolution
res_X    = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);

% semi-implicit update of phase fraction densities
upd_X    = - alpha*res_X*dt/a1;
X        = X + upd_X;
X        = max(rho.*eps,min(rho.*(1-eps), X ));
M        = rho - X;

% update phase fractions
x        = X./rho; 
m        = M./rho;

% record timing
XEtime   = XEtime + toc;
