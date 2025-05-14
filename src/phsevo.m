%*****  PHASE FRACTION EVOLUTION  *****************************************

tic;

%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_M;

% xtal diffusion rate
dffn_X   = diffus(chi,X.*kx,h,[1,2],BCD);

if ~bnchm

% smoothed random noise source
if iter==1
    rns   = randn(Nz,Nx);
    for i = 1:round(smth)  % apply same smoothing as for initial condition
        rns = rns + diffus(rns,1/8*ones(size(rns)),1,[1,2],{'periodic','periodic'});
    end
    rns  = (rns - mean(rns(:)))./std(rns(:)); % normalise to mean=0, var=1
end
var_rns  = 4/3*pi*d0^3.*chi.*kwx./(h+Delta_sgr)^3./(dt+Delta_sgr./vx);     % variance of random noise source
rns_X    = rho.*sqrt(var_rns)./h.*rns;                                     % approximate random flux divergence
rns_X    = rns_X - mean(rns_X(:));                                         % ensure global mass conservation

% boundary phase change rate
tau_x    = (h/2)./(norm(Wx+eps,'fro')./sqrt(length(Wx(:))));
Gx       = max(0,Da.*(xeq.*(1+rns./100)-x).*rho./tau_x.*topshape);
Gm       = -Gx;

end

% total rates of change
dXdt     = advn_X + dffn_X + Gx + rns_X;

% residual of phase density evolution
res_X    = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);

% semi-implicit update of phase fraction densities
upd_X    = - alpha*res_X*dt/a1;

X        = X + upd_X;
X        = max(rhox0.*eps^0.5,X);
M        = rho - X;

%***  update phase fractions and component concentrations

% update phase fractions
x        = X./rho; 
m        = M./rho;

% record timing
TCtime   = TCtime + toc;
