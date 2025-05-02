%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_M;

% phase advection rates
dffn_X   = diffus(chi,rho.*kx,h,[1,2],BCD);
dffn_M   = diffus(mu ,rho.*kx,h,[1,2],BCD);

% phase change rates
% if iter==1
    tau_x = (h/2)./mean(abs(wx(2,:)));
    Gx    = (Gx + max(0,Da.*(xeq-x).*rho./tau_x.*topshape))/2;
    Gm    = -Gx;
% end

% total rates of change
dXdt   = advn_X + dffn_X + Gx;
dMdt   = advn_M + dffn_M + Gm;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
upd_X = max(-X/2, - alpha*res_X*dt/a1 + beta*upd_X );
upd_M = max(-M/2, - alpha*res_M*dt/a1 + beta*upd_M );

X     = X + upd_X;
M     = M + upd_M;

% get dynamically evolving mixture density 
RHO = X+M;

%***  update phase fractions and component concentrations

% update phase fractions
x = X./RHO; 
m = M./RHO;

% record timing
TCtime = TCtime + toc;% - eqtime;
