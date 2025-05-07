%*****  FLUID MECHANICS SOLVER  *******************************************

tic;

if ~bnchm && step>0 && ~restart

%***  update mixture mass density
drhodt  = advn_rho;

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
upd_rho = - alpha*res_rho./b1./rho;
dV      = dV + upd_rho;  % correct volume source term by scaled residual

dVmean  = mean(dV,'all');

UBG     = - 0*dVmean./2 .* (L/2-XXu);
WBG     = - 2*dVmean./2 .* (0/2-ZZw);

end


%% assemble coefficients for matrix velocity diagonal and right-hand side

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

% assemble coefficients of z-stress divergence
    
% left boundary
ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% top boundary
ii  = MapW(1,2:end-1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii)) + WBG(1,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapW(end,2:end-1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii)) + WBG(end,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 =  etaco(2:end-1,1:end-1);   EtaC2 =  etaco(2:end-1,2:end);
EtaP1 =  eta  (1:end-1,:      );   EtaP2 =  eta  (2:end,:      );

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa  = a1.*rhow(2:end-1,:)./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = ddz(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right

% z-RHS vector
advn_mz = advect(rhow(2:end-1,:).*W(2:end-1,2:end-1),(U(2:end-2,:)+U(3:end-1,:))/2,(W(1:end-1,2:end-1)+W(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
rr  = + Drho(2:end-1,:) .* g0 ...
      + (a2.*rhoWo(2:end-1,:)+a3.*rhoWoo(2:end-1,:))/dt ...
      - advn_mz;
if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];

% assemble coefficients of x-stress divergence

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+top];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+bot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii    = MapU(2:end-1,:);
EtaC1 = etaco(1:end-1,:    );  EtaC2 = etaco(2:end,:    );
EtaP1 = eta  (:,icx(1:end-1));  EtaP2 = eta  (:,icx(2:end));

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,ifx(1:end-2)); jj2 = MapU(2:end-1,ifx(3:end)); jj3 = MapU(1:end-2,ifx(2:end-1)); jj4 = MapU(3:end,ifx(2:end-1));
aa  = (a1+gamma).*rhou./dt;

IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = ddx(rho(:,icx),h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); jj3 = MapW(2:end,1:end-1); jj4 = MapW(2:end,2:end);
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right

% x-RHS vector
advn_mx = advect(rhou.*U(2:end-1,:),(U(2:end-1,ifx(1:end-1))+U(2:end-1,ifx(2:end)))/2,(W(:,1:end-1)+W(:,2:end))/2,h,{ADVN,''},[1,2],BCA);
advn_mx(:,[1 end]) = repmat((advn_mx(:,1)+advn_mx(:,end))/2,1,2);
rr  = + (a2.*rhoUo+a3.*rhoUoo)/dt ...
    - advn_mx;
if bnchm
    rr = rr + src_U_mms(2:end-1,:);
end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];

% assemble coefficient matrix & right-hand side vector
KV  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV  = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator

if ~exist('GG','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    % coefficients for z-gradient
    ii  = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    % coefficients for x-gradient
    ii  = MapU(2:end-1,:);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    % assemble coefficient matrix
    GG  = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence operator

if ~exist('DD','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    %internal points
    ii  = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying velocities U, W
    %          left U          ||           right U       ||           top W           ||          bottom W
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below

    % assemble coefficient matrix
    DD  = sparse(IIL,JJL,AAL,NP,NW+NU);
end


%% assemble coefficients for matrix pressure diagonal and right-hand side

% if ~exist('KP','var') || bnchm || lambda1+lambda2>0
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    % boundary points
    ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
    jj1 = ii;
    jj2 = [MapP(2,:).'; MapP(end-1,:).'];
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    ii  = [MapP(2:end-1,1    ); MapP(2:end-1,end)]; % left & right
    jj1 = ii;
    jj2 = [MapP(2:end-1,end-1); MapP(2:end-1,2    )];
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
if ~exist('KP','var') || bnchm || lambda1+lambda2>0

    % internal points
    ii  = MapP(2:end-1,2:end-1);
    jj1 = MapP(1:end-2,2:end-1);
    jj2 = MapP(3:end-0,2:end-1);
    jj3 = MapP(2:end-1,1:end-2);
    jj4 = MapP(2:end-1,3:end-0);

    % coefficients multiplying matrix pressure P
    aa  = zeros(size(ii)) + lambda1*eps*h^2./eta;
    IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
    
    kP  = lambda2*h^2./eta;
    kP1 = (kP(icz(1:end-2),:).*kP(icz(2:end-1),:)).^0.5;   kP2 = (kP(icz(2:end-1),:).*kP(icz(3:end-0),:)).^0.5;
    kP3 = (kP(:,icx(1:end-2)).*kP(:,icx(2:end-1))).^0.5;   kP4 = (kP(:,icx(2:end-1)).*kP(:,icx(3:end-0))).^0.5;

    aa  = (kP1+kP2+kP3+kP4)/h^2;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL;-aa(:)     ];      % P on stencil centre
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; kP1(:)/h^2];      % P one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; kP2(:)/h^2];      % P one below
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; kP3(:)/h^2];      % P one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; kP4(:)/h^2];      % P one below
    
end

% assemble coefficient matrix
KP  = sparse(IIL,JJL,AAL,NP,NP);

% RHS
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

ii  = MapP(2:end-1,2:end-1);

rr  = dV;       % add volume source term
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);

if bnchm
    % fix P = P_mms in middle of domain
    nzp = round((Nz+2)/2);
    nxp = round((Nx+2)/2);
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = 1;
    DD(np0,:  ) = 0;
    RP(np0    ) = P_mms(nzp,nxp);
    
    % fix U = U_mms in middle of domain
    nzu = round((Nz+2)/2);
    nxu = round((Nx+2)/2);
    nu0 = MapU(nzu,nxu);
    KV(nu0,:  ) = 0;
    KV(nu0,nu0) = 1;
    GG(nu0,:  ) = 0;
    RV(nu0    ) = U_mms(nzu,nxu);
else
    % set P = 0 in fixed point
    nzp = 2;
    nxp = round((Nx+2)/2);
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = 1;
    DD(np0,:  ) = 0;
    RP(np0    ) = 0;
end


%% assemble and scale global coefficient matrix and right-hand side vector

LL  = [ KV   GG  ; ...
        DD   KP ];

RR  = [RV; RP];

etagh = ones(size(P));  etagh(2:end-1,2:end-1) = eta;
SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./(SCL + sqrt([zeros(NU+NW,1); h./etagh(:)])) ));

FF  = LL*SOL - RR;
FF  = SCL*FF;

LL  = SCL*LL*SCL;
RR  = SCL*RR;


%% Solve linear system of equations for vx, vz, P

SOL = SCL*(LL\RR);  % update solution

% map solution vector to 2D arrays
W = full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
U = full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
P = full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure
P = P - mean(mean(P(2:end-1,2:end-1)));             % reduce pressure by mean


%% Update phase segregation speeds
if ~bnchm

    % advection of xtal segregation speed
    advn_wx(2:end-1,2:end-1) = - advect(wx(2:end-1,2:end-1),(Ux(2:end-2,:)+Ux(3:end-1,:))/2,(Wx(1:end-1,2:end-1)+Wx(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
    advn_wx([1,end],:) = min(1,1-[top;bot]).*advn_wx([2,end-1],:);
    advn_wx(:,[1 end]) = advn_wx(:,[end-1 2]);

    % terminal xtal segregation speed
    wx0(2:end-1,2:end-1) = Drhox(2:end-1,:).*g0.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5; % solid segregation speed
    wx0([1,end],:) = min(1,1-[top;bot]).*wx0([2,end-1],:);
    wx0(:,[1 end]) = wx0(:,[end-1 2]);

    tw = max(dt,d0.^2.*rho./eta0);
    Gw = (wx0-wx)./((tw(icz(1:end-1),icx)+tw(icz(2:end),icx))./2);

    % total rate of change
    dwxdt  = advn_wx + Gw;

    % residual of xtal settling speed evolution
    res_wx = (a1*wx-a2*wxo-a3*wxoo)/dt - (b1*dwxdt + b2*dwxdto + b3*dwxdtoo);

    % semi-implicit update of xtal settling speed
    upd_wx = - alpha*res_wx*dt/a1 + beta*upd_wx;

    % update xtal settling speed
    wx     = wx + upd_wx;

    xw   = (x(icz(1:end-1),icx)+x(icz(2:end),icx))./2;
    mw   = (m(icz(1:end-1),icx)+m(icz(2:end),icx))./2;
    wm   = -xw./mw.*wx;

    % update phase velocities
    Wx  = W + wx;  % xtl z-velocity
    Ux  = U;       % xtl x-velocity
    Wm  = W + wm;  % mlt z-velocity
    Um  = U;       % mlt x-velocity

    
    %% update time step
    dtk = (h/2)^2/max(kx(:)); % diffusive time step size  
    dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:)]+eps));  % advective time step size
    dt  = (dt + min([1.1*dto,min(CFL*[dtk,dta]),dtmax]))/2; % time step size
end

FMtime = FMtime + toc;