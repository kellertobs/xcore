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

% generate smooth random noise (once per timestep)
if iter==1
    rns   = randn(Nz,Nx);
    for i = 1:smth  % apply same smoothing as for initial condition
        rns = rns + diffus(rns,1/8*ones(size(rns)),1,[1,2],{'periodic','periodic'});
    end
    rns  = (rns - mean(rns(:)))./std(rns(:)); % normalise to mean=0, var=1
    isx  = randi(Nx,1); isz = randi(Nz,1);
    rns1 = circshift(circshift(rns,isx,2),isz,1);
    isx  = randi(Nx,1); isz = randi(Nz,1);
    rns2 = circshift(circshift(rns,isx,2),isz,1);
    isx  = randi(Nx,1); isz = randi(Nz,1);
    rns3 = circshift(circshift(rns,isx,2),isz,1);
    isx  = randi(Nx,1); isz = randi(Nz,1);
    rns4 = circshift(circshift(rns,isx,2),isz,1);
end

% random noise source for particle fluctuations
var_rnss = ks.*(Delta_sgr./(h+Delta_sgr)).^3./(dt+Delta_sgr./vx);          % variance of random noise source
rnsz     = (rns1(icz(1:end-1),:)+rns1(icz(2:end),:))/2;
rnsx     = (rns2(:,icx(1:end-1))+rns2(:,icx(2:end)))/2;
rns_Xs   = X.*sqrt(var_rnss).*(ddz(rnsz,h)+ddx(rnsx,h));                   % random flux divergence
rns_Xs   = rns_Xs - mean(rns_Xs(:));                                       % ensure global mass conservation

% random noise source for subgrid eddy fluctuations
var_rnse = ke.*(Delta_cnv./(h+Delta_cnv)).^3./(dt+Delta_cnv./V );          % variance of random noise source
rnsz     = (rns3(icz(1:end-1),:)+rns3(icz(2:end),:))/2;
rnsx     = (rns4(:,icx(1:end-1))+rns4(:,icx(2:end)))/2;
rns_Xe   = X.*sqrt(var_rnse).*(ddz(rnsz,h)+ddx(rnsx,h));                   % random flux divergence
rns_Xe   = rns_Xe - mean(rns_Xe(:));                                       % ensure global mass conservation

% boundary phase change rate
tau_x    = (h/2)./(W0 + w0) + dt;
Gx       = max(0,Da.*(xeq.*(1+rns./10)-x).*rho./tau_x.*topshape);
Gm       = -Gx;

end

% total rates of change
dXdt     = advn_X + dffn_X + Gx + rns_Xs + rns_Xe;

% residual of phase density evolution
res_X    = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);

% semi-implicit update of phase fraction densities
upd_X    = - alpha*res_X*dt/a1;

X        = X + upd_X;
X        = max(rho.*eps^0.5,X);
M        = rho - X;

%***  update phase fractions and component concentrations

% update phase fractions
x        = X./rho; 
m        = M./rho;

% record timing
TCtime   = TCtime + toc;
