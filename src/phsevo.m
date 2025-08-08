%*****  PHASE FRACTION EVOLUTION  *****************************************

tic;

%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],{(  [xeq,x0]).*rhom0,'periodic'});
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],{(1-[xeq,x0]).*rhom0,'periodic'});
advn_rho = advn_X+advn_M;

% xtal diffusion rate
dffn_X   = diffus(chi,X.*kx,h,[1,2],BCD);

if ~bnchm
% boundary phase change rate
tau_x    = h./(W0 + w0) + dt;
Gx       = R.*max(0,xeq-x).*rho./tau_x.*bndshape;%[1;zeros(Nz-1,1)];
end

% total rates of change
dXdt     = advn_X + dffn_X + Gx;

% residual of phase density evolution
res_X    = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);

% semi-implicit update of phase fraction densities
upd_X    = - alpha*res_X*dt/a1;

X        = X + upd_X;
X        = max(rho.*eps,min(rho.*(1-eps), X ));
M        = rho - X;

%***  update phase fractions and component concentrations

% update phase fractions
x        = X./rho; 
m        = M./rho;

% record timing
TCtime   = TCtime + toc;
