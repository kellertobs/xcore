%*****  FLUID MECHANICS SOLVER  *******************************************


%% update mass flux divergence source term

if ~bnchm && step>0 && ~restart

%***  update mixture mass density
drhodt  = advn_rho;

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
upd_MFS = - alpha*res_rho./b1;
MFS     = MFS + upd_MFS;  % correct volume source term by scaled residual

dVmean = mean(MFS,'all');

UBG    = - 0*dVmean./2 .* (L/2-XXu);
WBG    = - 2*dVmean./2 .* (   -ZZw);

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
ii  = MapW(end,2:end-1); jj1 = ii; jj2 = MapW(end-1,2:end-1); jj3 = MapU(end-1,2:end); jj4 = MapU(end-1,1:end-1);
rho1 = rhow(end  ,:    );
rho2 = rhow(end-1,:    );
rho3 = rhou(end,2:end  );
rho4 = rhou(end,1:end-1);
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;+     rho1(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-open*rho2(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+open*rho3(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-open*rho4(:)/h];
aa  = open.*MFS(end,:) + closed.*rho1.*WBG(end,2:end-1)/h;
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


%% assemble coefficients for divergence of matrix mass flux (DM)

IIL = [];       % equation indeces into A
JJL = [];       % variable indeces into A
AAL = [];       % coefficients for A

%internal points
ii  = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
rho1 = rhou(:,1:end-1);
rho2 = rhou(:,2:end  );
rho3 = rhow(1:end-1,:);
rho4 = rhow(2:end  ,:);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -rho1(:)/h];  % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; +rho2(:)/h];  % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -rho3(:)/h];  % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; +rho4(:)/h];  % W one below

% Assemble coefficient matrix
DM  = sparse(IIL,JJL,AAL,NP,NW+NU);


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

    % coefficients multiplying matrix pressure P
    aa  = zeros(size(ii));
    
end

% assemble coefficient matrix
KP  = sparse(IIL,JJL,AAL,NP,NP);

% RHS
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

ii  = MapP(2:end-1,2:end-1);

rr  = MFS;       % add volume source term
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
    if open
        nzp = Nz+1;
        nxp = 1:Nx+2;
    else
        nzp = 2;
        nxp = round(Nx/2);
    end
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = speye(length(np0));
    if closed
        DD(np0,:  ) = 0;
        RP(np0    ) = 0;
    end
end


%% assemble and scale global coefficient matrix and right-hand side vector

LL  = [KV GG  ; ...
       DM KP ];

RR  = [RV; RP;];

scl = ones(size(P));  scl(2:end-1,2:end-1) = rho./eta;
SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./(SCL + sqrt([zeros(NU+NW,1); 1./scl(:)])) ));

FF  = SCL*(LL*SOL - RR);
LL  = SCL*LL*SCL;


%% Solve linear system of equations for vx, vz, P

if ~exist('pcol','var'); pcol = colamd(LL); end % get column permutation for sparsity pattern once per run
dLL         = decomposition(LL(:,pcol), 'lu');  % get LU-decomposition for consistent performance of LL \ RR
UPD(pcol,1) = dLL \ FF;                         % solve permuted decomposed system
UPD         = SCL*UPD;

% map solution vector to 2D arrays
upd_W = -full(reshape(UPD(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
upd_U = -full(reshape(UPD(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
upd_P = -full(reshape(UPD(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure
upd_P = upd_P - mean(mean(upd_P(2:end-1,2:end-1)));           % reduce pressure by mean

% update solution
W = W + upd_W;
U = U + upd_U;
P = P + upd_P;
SOL = [W(:);U(:);P(:)];


%% Update phase segregation speeds
if ~bnchm && step>=1

    % terminal xtal segregation speed for comparison
    wx0(2:end-1,2:end-1) = d0^2./etasw(2:end-1,:).*Drhox(2:end-1,:).*g0; % solid segregation speed
    wx0([1,end],:) = min(1,1-[top;bot]).*wx0([2,end-1],:);
    wx0(:,[1 end]) = wx0(:,[end-1 2]);

    % % advection of crystal momentum
    % % advn_Mx = - advect(Mx(2:end-1,:),(Ux(2:end-2,:)+Ux(3:end-1,:))/2,(Wx(1:end-1,2:end-1)+Wx(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
    % advn_Mx = - advect(Mx(2:end-1,:),(U(2:end-2,:)+U(3:end-1,:))/2,(W(1:end-1,2:end-1)+W(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
    % 
    % % crystal momentum transfer from mixture
    % Gvx     = - Cxw(2:end-1,:) .* wx(2:end-1,2:end-1);
    % 
    % % crystal momentum source
    % Qvx     = + chiw(2:end-1,:) .* Drhox(2:end-1,:) .* g0;
    % 
    % % get total rate of change
    % dMxdt(2:end-1,:) = advn_Mx + Gvx + Qvx;
    % 
    % % residual of xtal partial momentum deviation
    % res_Mx = (a1*Mx-a2*Mxo-a3*Mxoo)/dt - (b1*dMxdt + b2*dMxdto + b3*dMxdtoo);
    % 
    % % semi-implicit update of xtal momentum
    % upd_Mx = - alpha*res_Mx*dt/a1;
    % 
    % % update crystal momentum with dynamic under-relaxation
    % tau_p  = chiw.*rhow./Cxw;
    % relax  = dt ./ (dt + tau_p);
    % upd_Mx = (1-relax).*upd_Mx;
    % 
    % Mx     = Mx + upd_Mx;
    % 
    % % update crystal settling speed
    % wx(:,2:end-1) = Mx./Xw;
    % wx([1,end],:) = min(1,1-[top;bot]).*wx([2,end-1],:);
    % wx(:,[1 end]) = wx(:,[end-1 2]);

    wx  = wx0;
    wm  = -xw(:,icx)./mw (:,icx).*wx;

    % generate smooth random noise (once per timestep)
    if iter==1
        rweo  = rwe;
        rwe   = randn(Nz+1,Nx);
        for i = 1:ceil(smth)  % apply same smoothing as for initial condition
            ksmth = min(1,smth-i+1);
            rwe = rwe + diffus(rwe,ksmth/8*ones(size(rwe)),1,[1,2],{'periodic','periodic'});
        end
        rwe  = (rwe - mean(rwe(:)))./std(rwe(:)); % normalise to mean=0, var=1
        if step==1; rweo = rwe; end

        rueo = rue;
        rue  = randn(Nz,Nx+1);
        for i = 1:ceil(smth)  % apply same smoothing as for initial condition
            ksmth = min(1,smth-i+1);
            rue = rue + diffus(rue,ksmth/8*ones(size(rue)),1,[1,2],{'periodic','periodic'});
        end
        rue  = (rue - mean(rue(:)))./std(rue(:)); % normalise to mean=0, var=1
        if step==1; rueo = rue; end

        rwso = rws;
        ruso = rus;
        isx  = randi(Nx,1); isz = randi(Nz,1);
        rws  = circshift(circshift(rwe,isx,2),isz,1);
        isx  = randi(Nx,1); isz = randi(Nz,1);
        rus  = circshift(circshift(rue,isx,2),isz,1);
        if step==1; rwso = rws; end
        if step==1; rwso = rws; end

        tau  = sqrt(smth)/2*dt;
        rwe  = rweo + (rwe-rweo).*dt./tau;
        rue  = rueo + (rue-rueo).*dt./tau;
        rws  = rwso + (rws-rwso).*dt./tau;
        rus  = ruso + (rus-ruso).*dt./tau;
        rwe  = (rwe - mean(rwe(:)))./std(rwe(:)); % normalise to mean=0, var=1
        rue  = (rue - mean(rue(:)))./std(rue(:)); 
        rws  = (rws - mean(rws(:)))./std(rws(:));
        rus  = (rus - mean(rus(:)))./std(rus(:));
      
    end

    % random noise source variance
    xie  = Xi*sqrt(ke.*(Delta_cnv./(h+Delta_cnv)).^3./(dt+Delta_cnv./V )); % eddy noise speed
    xis  = Xi*sqrt(ks.*(Delta_sgr./(h+Delta_sgr)).^3./(dt+Delta_sgr./vx)); % segregation noise speed

    xisw = (xis(icz(1:end-1),icx) + xis(icz(2:end),icx))./2;
    xisu = (xis(icz,icx(1:end-1)) + xis(icz,icx(2:end)))./2;

    xiew = (xie(icz(1:end-1),icx) + xie(icz(2:end),icx))./2;
    xieu = (xie(icz,icx(1:end-1)) + xie(icz,icx(2:end)))./2;

    xiw = xisw .* rws(:,icx) + xiew .* rwe(:,icx); 
    xiu = xisu .* rus(icz,:) + xieu .* rue(icz,:);

    xiw = xiw - mean(xiw(:));
    xiu = xiu - mean(xiu(:)); 

    xiw([1,end],:) = xiw([2,end-1],:);
    xiw(:,[1 end]) = xiw(:,[end-1 2]);
    xiu([1 end],:) = xiu([2 end-1],:);
    xiu(:,[1 end]) = repmat(mean(xiu(:,[1 end]),2),1,2);

    xixw  = xiw;
    xixu  = xiu;
    ximw  = -xw(:,icx)./mw (:,icx).*xixw;
    ximu  = -xu(icz,:)./muu(icz,:).*xixu;

    % update phase velocities
    Wx  = W + wx + xixw;  % xtl z-velocity
    Ux  = U + 0  + xixu;  % xtl x-velocity
    Wm  = W + wm + ximw;  % mlt z-velocity
    Um  = U + 0  + ximu;  % mlt x-velocity

end

% record timing
FMtime = FMtime + toc;