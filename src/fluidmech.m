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

MFSmean = mean(MFS,'all');

MFBG    = MFSmean .* ZZw;

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
aa  = open.*MFS(end,:) + closed.*MFBG(end,2:end-1)/h;
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
        nzp = round(Nz/2);
        nxp = round(Nx/2);
    end
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = speye(length(np0));
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

        % Generate new white noise
        rsw = randn(Nz+1, Nx+0);
        rsu = randn(Nz+0, Nx+1);
        rew = randn(Nz+1, Nx+0);
        reu = randn(Nz+0, Nx+1);
        pse = randn(Nz+1, Nx+1);

        rsw = fft2(rsw);
        rsu = fft2(rsu);
        rew = fft2(rew);
        reu = fft2(reu);
        pse = fft2(pse);

        % Filter white noise spatially to decorrelation length
        rsw = real(ifft2(Gkws .* rsw));
        rsu = real(ifft2(Gkus .* rsu));
        rew = real(ifft2(Gkwe .* rew));
        reu = real(ifft2(Gkue .* reu));
        pse = real(ifft2(Gkps .* pse));

        % rescale to unit RMS speed
        ss  = sqrt(mean(rsw(:).^2) + mean(rsu(:).^2));
        se  = sqrt(mean(rew(:).^2) + mean(reu(:).^2));
        rsw = rsw / ss;
        rsu = rsu / ss;
        rew = rew / se;
        reu = reu / se;
        pse = pse / std(pse(:));

    end

    % noise decorrelation time
    taue = elle/2./(V +eps);
    taus = ells/2./(vx+eps);

    % temporal evolution factor
    Fe   = exp(-dt./taue);
    Fs   = exp(-dt./taus);

    Few = (Fe(icz(1:end-1),icx)+Fe(icz(2:end),icx))/2;
    Feu = (Fe(icz,icx(1:end-1))+Fe(icz,icx(2:end)))/2;

    Fec = (Fe(icz(1:end-1),icx(1:end-1))+Fe(icz(1:end-1),icx(2:end)))/4 ...
        + (Fe(icz(2:end  ),icx(1:end-1))+Fe(icz(2:end  ),icx(2:end)))/4;

    Fsw = (Fs(icz(1:end-1),icx)+Fs(icz(2:end),icx))/2;
    Fsu = (Fs(icz,icx(1:end-1))+Fs(icz,icx(2:end)))/2;

    % random noise source amplitude
    bndtaper = (1 - exp((-ZZ+h/2)/max(h,elle)) - closed.*exp(-(D-ZZ-h/2)/max(h,elle)));
    sgs   = Xi * sqrt(chi.*ks./taus .* (ells./(ells+h)).^3) .* bndtaper; % segregation noise speed
    sge   = Xi * sqrt(     ke./taue .* (elle./(elle+h)).^3) .* bndtaper; % eddy mixture noise speed
    sgex  = Xi * sqrt(chi.*ke./taue .* (elle./(elle+h)).^3) .* bndtaper; % eddy crystal noise speed

    sgsw  = (sgs(icz(1:end-1),icx) + sgs(icz(2:end),icx))./2;
    sgsu  = (sgs(icz,icx(1:end-1)) + sgs(icz,icx(2:end)))./2; 

    sgexw = (sgex(icz(1:end-1),icx) + sgex(icz(2:end),icx))./2;
    sgexu = (sgex(icz,icx(1:end-1)) + sgex(icz,icx(2:end)))./2;

    sgec  = (sge(icz(1:end-1),icx(1:end-1))+sge(icz(1:end-1),icx(2:end)))/4 ...
          + (sge(icz(2:end  ),icx(1:end-1))+sge(icz(2:end  ),icx(2:end)))/4;

    % Ornstein–Uhlenbeck time update
    xisw  =  Fsw .* xiswo  + sqrt(1 - Fsw.^2) .* sgsw  .* rsw(:,icx);
    xisu  =  Fsu .* xisuo  + sqrt(1 - Fsu.^2) .* sgsu  .* rsu(icz,:);
    xiexw =  Few .* xiexwo + sqrt(1 - Few.^2) .* sgexw .* rew(:,icx);
    xiexu =  Feu .* xiexuo + sqrt(1 - Feu.^2) .* sgexu .* reu(icz,:);  
    psie  =  Fec .* psieo  + sqrt(1 - Fec.^2) .* sgec  .* pse;
    xieu  =  ddz(psie,1); xieu = xieu(icz,:);
    xiew  = -ddx(psie,1); xiew = xiew(:,icx);

    % update phase noise speeds
    xiwx  = xisw + xiexw;
    xiux  = xisu + xiexu;
    xiwm  = -xw(:,icx)./mw (:,icx).*xiwx;
    xium  = -xu(icz,:)./muu(icz,:).*xiux;

    % update phase velocities
    Wx  = W + wx + xiwx + xiew;  % xtl z-velocity
    Ux  = U + 0  + xiux + xieu;  % xtl x-velocity
    Wm  = W + wm + xiwm + xiew;  % mlt z-velocity
    Um  = U + 0  + xium + xieu;  % mlt x-velocity

end

% record timing
FMtime = FMtime + toc;