%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_M;

% phase change rates
%if iter==1; ishft = randi(Nx); end
Gx  =  Da.*(xeq-x).*rho./dt.*topshape;%.*(1+circshift(rp,ishft,2)/3);
Gm  =  -Gx;

% total rates of change
dXdt   = advn_X + Gx;
dMdt   = advn_M + Gm;

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
