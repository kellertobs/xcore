% scaling analysis of crystal-driven convection


%% laminar case (fReL = fReL = 0)

clear;
syms chi0 Dchi0 Drho0 rho0 g0 Ds0 d0 eta0 L0 l0

fReL = 0;
fRel = 0;

% turbulent eddy viscosity
clear W0; syms W0
ke0   = W0/Ds0/2*L0^2;
etae0 = eta0 + fReL*ke0*rho0;

% turbulent settling viscosity
clear w0; syms w0
ks0   = w0*l0;
etas0 = eta0 + fRel*ks0*rho0;

% Speed scale for crystal-driven convection
eq2W = W0 - Dchi0*Drho0*g0*Ds0^2/etae0 == 0;
W0   = solve(eq2W,W0); W0l = W0;

% Speed scale for crystal settling
eq2w = w0 - Drho0*g0*d0^2/etas0 == 0;
w0   = solve(eq2w,w0); w0l = w0;

% update diffusivities
ke0  = W0/Ds0/2*L0^2;
ks0  = w0*l0;
kx0  = ks0 + fReL*ke0;

% dimensionless numbers
Rc0  = W0/w0;                % Convection number
Ra0  = W0*Ds0/kx0;            % Rayleigh number
ReD0 = W0*Ds0/(etae0/rho0);   % Domain Reynolds number
Red0 = w0*d0/(etas0/rho0);   % Particle Reynolds number

fprintf(1,'\n\n*****  Scaling analysis for laminar case \n\n');
fprintf(1,'       w0    = %s \n'  ,simplify(w0l));
fprintf(1,'       W0    = %s \n\n',simplify(W0l));

fprintf(1,'       ks0   = %s \n'  ,simplify(ks0));
fprintf(1,'       ke0   = %s \n'  ,simplify(ke0));
fprintf(1,'       kx0   = %s \n\n',simplify(kx0));

fprintf(1,'       etas0 = %s \n'  ,simplify(etas0));
fprintf(1,'       etae0 = %s \n\n',simplify(etae0));

fprintf(1,'       Rc0   = %s \n'  ,simplify(Rc0));
fprintf(1,'       Ra0   = %s \n'  ,simplify(Ra0));
fprintf(1,'       ReD0  = %s \n'  ,simplify(ReD0));
fprintf(1,'       Red0  = %s \n\n',simplify(Red0));


%% turbulent case (fReL = fRel = 1, eta0 = 0)

syms chi0 Dchi0 Drho0 rho0 g0 D0 Ds0 d0 eta0 Ls0 l0 C0

fReL = 1;
fRel = 1;

% turbulent eddy viscosity
clear W0; syms W0
ke0   = W0*Ls0;
etae0 = 0 + fReL*ke0*rho0;

% turbulent settling viscosity
clear w0; syms w0
ks0   = w0*l0;
etas0 = 0 + fRel*ks0*rho0;

% Speed scale for inertial momentum build-up
eq2W = W0 - D0/(rho0*W0/(Dchi0*Drho0*g0)) == 0;
W0i  = solve(eq2W,W0);
W0i  = W0i(1);

% Speed scale for turbulent crystal-driven convection
eq2W = W0 - Dchi0*Drho0*g0*Ds0^2/etae0 == 0;
W0t  = solve(eq2W,W0);
W0t  = W0t(1);

% Speed scale for crystal settling
eq2w = w0 - Drho0*g0*d0^2/etas0 == 0;
w0t  = solve(eq2w,w0);
w0t  = w0t(1);


clear W0 w0 kx0 ks0 ke0 etas0 etae0; syms W0 w0;

% update diffusivities
ke0  = W0*Ls0;
ks0  = w0*l0;
kx0  = ks0 + fReL*ke0;

% update viscosities
etas0 = 0 + fRel*ks0*rho0;
etae0 = 0 + fReL*ke0*rho0;

% dimensionless numbers
Ri0  = W0i/W0t;
Rc0  = W0i/w0t;              % Convection number
Ra0  = W0*Ds0/kx0;           % Rayleigh number
ReD0 = W0*Ds0/(etae0/rho0);  % Domain Reynolds number
Red0 = w0*d0/(etas0/rho0);   % Particle Reynolds number

fprintf(1,'\n\n*****  Scaling analysis for turbulent case \n\n');
fprintf(1,'       w0t   = %s \n'  ,simplify(w0t));
fprintf(1,'       W0t   = %s \n'  ,simplify(W0t));
fprintf(1,'       W0i   = %s \n\n',simplify(W0i));

fprintf(1,'       ks0   = %s \n'  ,simplify(ks0));
fprintf(1,'       ke0   = %s \n'  ,simplify(ke0));
fprintf(1,'       kx0   = %s \n\n',simplify(kx0));

fprintf(1,'       etas0 = %s \n'  ,simplify(etas0));
fprintf(1,'       etae0 = %s \n\n',simplify(etae0));

fprintf(1,'       Rc0   = %s \n'  ,simplify(Rc0));
fprintf(1,'       Ra0   = %s \n'  ,simplify(Ra0));
fprintf(1,'       ReD0  = %s \n'  ,simplify(ReD0));
fprintf(1,'       Red0  = %s \n\n',simplify(Red0));


%% general case

syms chi0 Dchi0 Drho0 rho0 g0 Ds0 d0 eta0 Ls0 l0 fReL fRel etae

% turbulent eddy viscosity
clear W0; syms W0
ke0   = W0*Ls0;
etae0 = eta0 + fReL*ke0*rho0;

% turbulent settling viscosity
clear w0; syms w0
ks0   = w0*l0;
etas0 = eta0 + fRel*ks0*rho0;

% Speed scale for crystal-driven convection
eq2W = W0 - Dchi0*Drho0*g0*Ds0^2/etae0 == 0;
W0g  = solve(eq2W,W0);
W0g  = W0g(2);

% Speed scale for crystal settling
eq2w = w0 - Drho0*g0*d0^2/etas0 == 0;
w0g  = solve(eq2w,w0);
w0g  = w0g(2);

clear W0 w0 kx0 ks0 ke0 etas0 etae0; syms W0 w0;

% update diffusivities
ke0  = W0*Ls0;
ks0  = w0*l0;
kx0  = ks0 + fReL*ke0;

% update viscosities
etas0 = eta0 + fRel*ks0*rho0;
etae0 = eta0 + fReL*ke0*rho0;

% dimensionless numbers
Rc0  = W0/w0;                 % Convection number
Ra0  = W0*Ds0/kx0;            % Rayleigh number
ReD0 = W0*Ds0/(etae0/rho0);   % Domain Reynolds number
Red0 = w0*d0 /(etas0/rho0);   % Particle Reynolds number

fprintf(1,'\n\n*****  Scaling analysis for general case \n\n');
fprintf(1,'       w0    = %s \n'  ,simplify(w0g));
fprintf(1,'       W0    = %s \n\n',simplify(W0g));

fprintf(1,'       ks0   = %s \n'  ,simplify(ks0));
fprintf(1,'       ke0   = %s \n'  ,simplify(ke0));
fprintf(1,'       kx0   = %s \n\n',simplify(kx0));

fprintf(1,'       etas0 = %s \n'  ,simplify(etas0));
fprintf(1,'       etae0 = %s \n\n',simplify(etae0));

fprintf(1,'       Rc0   = %s \n'  ,simplify(Rc0));
fprintf(1,'       Ra0   = %s \n'  ,simplify(Ra0));
fprintf(1,'       ReD0  = %s \n'  ,simplify(ReD0));
fprintf(1,'       Red0  = %s \n\n',simplify(Red0));


%% numerical evaluation
clear; close all;

fig1 = figure(1); clf; set(gcf,'Units','centimeters','Position',[6,6,20,15]);
fig2 = figure(2); clf; set(gcf,'Units','centimeters','Position',[8,8,20,15]);
fig3 = figure(3); clf; set(gcf,'Units','centimeters','Position',[10,10,20,15]);
fig4 = figure(4); clf; set(gcf,'Units','centimeters','Position',[12,12,18,15]);
load ../src/colmap/vik.mat; colmap = vik;
if ~isfolder('../out/scaling'); mkdir('../out/scaling'); end

etait = 0;
for eta0=10.^linspace(-1,5,4)
    etait = etait+1;

    lnshd = (6-(log10(eta0)+1))/8;
    blk   = [0 0 0]; wht = [1 1 1];
    red   = [0.8 0 0.1]; blu = [0.1 0 0.8];
    prp   = [0.5 0.1 0.5]; grn = [0.2 0.6 0.3];

    % parameter range
    d0   = 10.^linspace(-4,-1,121);
    D0   = 10.^linspace(0,6,121);
    [D0,d0] = meshgrid(D0,d0);

    Drho0 = 500;
    chi0  = 0.01;
    Dchi0 = chi0/20;
    rho0  = 2700;
    g0    = 10;

    Ds0  = D0/20;
    L0   = D0/400;
    h0   = D0/200;
    l0   = d0*10;

    indD   = find(D0(1,:)==10);
    indd   = find(d0(:,1)==0.01);
    etaref = 1e+1;


    % Laminar case

    % speed scale for crystal settling
    w0l   =        Drho0.*g0.*d0 .^2./eta0;

    % speed scale for crystal-driven convection
    W0l   = Dchi0.*Drho0.*g0.*Ds0.^2./eta0;

    % particle diffusivity
    k0    = w0l.*l0;
    kx0   = k0;

    xix0  =  sqrt(chi0.*k0./(l0./2./w0l).*(l0./(l0+h0)).^3);

    % dimensionless numbers
    Rc0   = W0l./w0l;                 % Convection number
    Ra0   = W0l.*Ds0./k0;             % Rayleigh number
    ReD0  = W0l.*Ds0./(eta0./rho0);   % Domain Reynolds number
    Red0  = w0l.*d0 ./(eta0./rho0);   % Particle Reynolds number

    % Re-dependent ramp factors
    digits = 32;
    ReL0  = W0l.*L0 ./(eta0./rho0);
    Rel0  = w0l.*l0 ./(eta0./rho0);
    fReL  = vpa(1-exp(-ReL0));
    fRel  = vpa(1-exp(-Rel0));

    set(0,'CurrentFigure',fig1)
    subplot(2,2,1)
    p11 = loglog(d0(:,1),w0l(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,2)
    p12 = loglog(D0(1,:),W0l(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,3)
    p18 = loglog(d0(:,1),Ra0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p13 = loglog(d0(:,1),Rc0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p14 = loglog(d0(:,1),Red0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu);
    subplot(2,2,4)
    p15 = loglog(D0(1,:),Ra0(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p16 = loglog(D0(1,:),ReD0(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn);
    p17 = loglog(D0(1,:),Rc0(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on

    set(0,'CurrentFigure',fig2)
    subplot(2,2,1)
    p21 = loglog(d0(:,1),k0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    subplot(2,2,1)
    p22 = loglog(d0(:,1),kx0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,2)
    p23 = loglog(D0(1,:),kx0(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,3)
    p24 = loglog(d0(:,1),xix0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    subplot(2,2,3)
    p25 = loglog(d0(:,1),xix0(:,indD),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,4)
    p26 = loglog(D0(1,:),xix0(indd,:),'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on

    % Turbulent case

    % Navier-Stokes speed scale for crystal settling
    w0t    = sqrt(Drho0.*g0.*d0.^2./(l0.*rho0));
    ks0   = w0t.*l0;
    etas0 = 0 + rho0.*ks0;

    % Navier-Stokes speed scale for crystal-driven convection
    W0t   = sqrt(2.*Dchi0.*Drho0.*g0.*Ds0.^3./(L0.^2.*rho0));
    W0i   = sqrt(   Dchi0.*Drho0.*g0.*D0            ./rho0);
    ke0   = W0i./Ds0./2.*L0.^2;
    etae0 = 0 + rho0.*ke0;

    kx0   = ks0 + ke0;

    xie0  =  sqrt(      ke0./(L0./2./W0i).*(L0./(L0+h0)).^3);
    xiex0 =  sqrt(chi0.*ke0./(L0./2./W0i).*(L0./(L0+h0)).^3);
    xis0  =  sqrt(chi0.*ks0./(l0./2./w0t).*(l0./(l0+h0)).^3);
    xix0  =  xis0 + xiex0;

    % dimensionless numbers for non-turbulent case
    Rc0   = W0i./w0t;                 % Turbulent Convection number
    Ra0   = W0i.*Ds0./kx0;            % Turbulent Rayleigh number
    ReD0  = W0i.*Ds0./(etae0./rho0);  % Turbulent Domain Reynolds number
    Red0  = w0t.*d0 ./(etas0./rho0);  % Turbulent Particle Reynolds number

    set(0,'CurrentFigure',fig1)
    subplot(2,2,1)
    p31 = loglog(d0(:,1),w0t(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,2)
    p32 = loglog(D0(1,:),W0i(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,3)
    p38 = loglog(d0(:,1),Ra0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p33 = loglog(d0(:,1),Rc0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p34 = loglog(d0(:,1),Red0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu);
    subplot(2,2,4)
    p35 = loglog(D0(1,:),Ra0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p36 = loglog(D0(1,:),ReD0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn);
    p37 = loglog(D0(1,:),Rc0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on

    set(0,'CurrentFigure',fig2)
    subplot(2,2,1)
    p41 = loglog(d0(:,1),ks0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    p42 = loglog(d0(:,1),kx0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,2)
    p43 = loglog(D0(1,:),ke0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p44 = loglog(D0(1,:),kx0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,3)
    p45 = loglog(d0(:,1),xis0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    p46 = loglog(d0(:,1),xix0(:,indD),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,4)
    p47 = loglog(D0(1,:),xie0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p48 = loglog(D0(1,:),xix0(indd,:),'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on

    % General case

    % Navier-Stokes speed scale for crystal settling
    w0    = double((sqrt(4.*Drho0.*g0.*rho0.*fRel.*l0.*d0.^2 + eta0.^2) - eta0)./(2.*fRel.*l0.*rho0));
    ks0   = w0.*l0;
    etas0 = double(eta0 + fRel.*rho0.*ks0);

    % Navier-Stokes speed scale for crystal-driven convection
    Ri0   = W0i./W0t;
    W0    = double(Ds0.*(sqrt(Ri0.^-2.*2.*Dchi0.*Drho0.*g0.*rho0.*fReL.*L0.^2.*Ds0  + eta0.^2) - eta0)./(Ri0.^-2.*fReL.*L0.^2.*rho0));
    ke0   = W0./Ds0./2.*L0.^2;
    etae0 = double(eta0 + fReL.*rho0.*ke0);

    kx0   = double(ks0 + fReL.*ke0);

    xie0  =  double(sqrt(      fReL.*ke0./(L0./2./W0).*(L0./(L0+h0)).^3));
    xiex0 =  double(sqrt(chi0.*fReL.*ke0./(L0./2./W0).*(L0./(L0+h0)).^3));
    xis0  =  sqrt(chi0.*      ks0./(l0./2./w0).*(l0./(l0+h0)).^3);
    xix0  =  xis0 + xiex0;

    % dimensionless numbers for non-turbulent case
    Rc0   = W0./w0;                     % General Convection number
    Ra0   = W0.*Ds0./kx0;               % General Rayleigh number
    ReD0  = W0.*Ds0./(etae0./rho0);     % General Domain Reynolds number
    Red0  = w0.*d0 ./(etas0./rho0);     % General Particle Reynolds number

    set(0,'CurrentFigure',fig1)
    subplot(2,2,1)
    p51 = loglog(d0(:,1),w0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,2)
    p1(etait) = loglog(D0(1,:),W0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,3)
    p58 = loglog(d0(:,1),Ra0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p53 = loglog(d0(:,1),Rc0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p54 = loglog(d0(:,1),Red0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu);
    subplot(2,2,4)
    p55 = loglog(D0(1,:),Ra0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*prp); axis tight; box on; hold on
    p56 = loglog(D0(1,:),ReD0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn);
    p57 = loglog(D0(1,:),Rc0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    drawnow;

    set(0,'CurrentFigure',fig2)
    subplot(2,2,1)
    p61 = loglog(d0(:,1),ks0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    p62 = loglog(d0(:,1),kx0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,2)
    p2(etait) = loglog(D0(1,:),double(fReL(indd,:).*ke0(indd,:)),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p64       = loglog(D0(1,:),kx0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,3)
    p65 = loglog(d0(:,1),xis0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blu); axis tight; box on; hold on
    p66 = loglog(d0(:,1),xix0(:,indD),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    subplot(2,2,4)
    p67 = loglog(D0(1,:),xie0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*red); axis tight; box on; hold on
    p68 = loglog(D0(1,:),xix0(indd,:),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*grn); axis tight; box on; hold on
    drawnow;

    set(0,'CurrentFigure',fig3); clf;
    subplot(2,2,1);
    imagesc(log10(D0(1,:)),log10(d0(:,1)),log10(Rc0)); axis xy tight; box on; hold on; colorbar('Ticklabelinterpreter','latex','FontSize',11); colormap(colmap);
    clim([-1,1].*max(abs(log10(Rc0(:)))));
    subplot(2,2,2);
    imagesc(log10(D0(1,:)),log10(d0(:,1)),log10(Red0)); axis xy tight; box on; hold on; colorbar('Ticklabelinterpreter','latex','FontSize',11); colormap(colmap);
    clim([-1,1].*max(abs(log10(Red0(:)))));
    subplot(2,2,3);
    imagesc(log10(D0(1,:)),log10(d0(:,1)),log10(Ra0)); axis xy tight; box on; hold on; colorbar('Ticklabelinterpreter','latex','FontSize',11); colormap(colmap);
    clim([-1,1].*max(abs(log10(Ra0(:)))));
    subplot(2,2,4);
    imagesc(log10(D0(1,:)),log10(d0(:,1)),log10(ReD0)); axis xy tight; box on; hold on; colorbar('Ticklabelinterpreter','latex','FontSize',11); colormap(colmap);
    clim([-1,1].*max(abs(log10(ReD0(:)))));

    subplot(2,2,1)
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','k','LineStyle',':','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','k','LineStyle',':','LineWidth',1);
    for iD=0:1:2
        for id=-4:1:-2
            plot(iD,id,'ko','LineWidth',1,'MarkerFaceColor','w')
        end
    end
    % plot(0.5,-1.5,'ko','LineWidth',1,'MarkerFaceColor','w')
    plot([6 6 6],[-3 -2 -1],'ko','LineWidth',1,'MarkerFaceColor','w')
    plot(6,-2,'wo','LineWidth',1,'MarkerFaceColor','k')
    plot(log10(1e1),log10(1e-2),'wo','LineWidth',1,'MarkerFaceColor','k')

    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    % xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
    ylabel('log$_{10}$ Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
    title('log$_{10}$ Convection No. Rc','Interpreter','latex','FontSize',13);
    text(0.84,0.91,'\textbf{(a)}','Interpreter','latex','FontSize',13,'Units','normalized')

    subplot(2,2,2)
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','k','LineStyle',':','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','k','LineStyle',':','LineWidth',1);
    for iD=0:1:2
        for id=-4:1:-2
            plot(iD,id,'ko','LineWidth',1,'MarkerFaceColor','w')
        end
    end
    % plot(0.5,-1.5,'ko','LineWidth',1,'MarkerFaceColor','w')
    plot([6 6 6],[-3 -2 -1],'ko','LineWidth',1,'MarkerFaceColor','w')
    plot(6,-2,'wo','LineWidth',1,'MarkerFaceColor','k')
    plot(log10(1e1),log10(1e-2),'wo','LineWidth',1,'MarkerFaceColor','k')

    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    % xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
    % ylabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
    title('log$_{10}$ Settling Reynolds No. Re$_d$','Interpreter','latex','FontSize',13);
    text(0.84,0.91,'\textbf{(b)}','Interpreter','latex','FontSize',13,'Units','normalized')

    subplot(2,2,3)
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','k','LineStyle',':','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','k','LineStyle',':','LineWidth',1);
    for iD=0:1:2
        for id=-4:1:-2
            plot(iD,id,'ko','LineWidth',1,'MarkerFaceColor','w')
        end
    end
    % plot(0.5,-1.5,'ko','LineWidth',1,'MarkerFaceColor','w')
    plot([6 6 6],[-3 -2 -1],'ko','LineWidth',1,'MarkerFaceColor','w')
    plot(6,-2,'wo','LineWidth',1,'MarkerFaceColor','k')
    plot(log10(1e1),log10(1e-2),'wo','LineWidth',1,'MarkerFaceColor','k')

    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    xlabel('log$_{10}$ Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
    ylabel('log$_{10}$ Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
    title('log$_{10}$ Convect. Rayleigh No. Ra','Interpreter','latex','FontSize',13);
    text(0.84,0.91,'\textbf{(c)}','Interpreter','latex','FontSize',13,'Units','normalized')

    subplot(2,2,4)
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','k','LineStyle',':','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','w','LineStyle','-','LineWidth',1);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','k','LineStyle',':','LineWidth',1);
    for iD=0:1:2
        for id=-4:1:-2
            plot(iD,id,'ko','LineWidth',1,'MarkerFaceColor','w')
        end
    end
    % plot(0.5,-1.5,'ko','LineWidth',1,'MarkerFaceColor','w')
    plot([6 6 6],[-3 -2 -1],'ko','LineWidth',1,'MarkerFaceColor','w')
    plot(6,-2,'wo','LineWidth',1,'MarkerFaceColor','k')
    plot(log10(1e1),log10(1e-2),'wo','LineWidth',1,'MarkerFaceColor','k')

    set(gca,'TickLabelInterpreter','latex','FontSize',11)
    xlabel('log$_{10}$ Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
    % ylabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
    title('log$_{10}$ Convect. Reynolds No. Re$_D$','Interpreter','latex','FontSize',13);
    text(0.84,0.91,'\textbf{(d)}','Interpreter','latex','FontSize',13,'Units','normalized','Color','w')
    drawnow;

    % define regime diagram

    % X = log10([Rc0(:),Ra0(:),ReD0(:),Red0(:)]);
    % [PC_C, PC_A, PC_V] = pca(X,'Algorithm','svd','Centered','on','VariableWeights','variance');
    % [Ic,Fc]  = kmeans(PC_A,4,'Replicates',10,'Display','final');
    % 
    % [~,is]  = sort(sum(Fc,2),'ascend');
    % Ics = zeros(size(Ic));
    % for kc = 1:size(Fc,1)
    %     Ics(Ic==is(kc),1) = kc;
    % end

    Ics = zeros(size(D0));
    Ics(Ra0<=1 & ReL0 <= 1) = 1;
    Ics(Ra0> 1 & ReL0 <= 1) = 2;
    Ics(Ra0<=1 & ReL0 >  1) = 3;
    Ics(Ra0> 1 & ReL0 >  1) = 4;
    % Ics(Ra0<=1 & ReL0 <= 1 & Rel0 >  1) = 5;
    % Ics(Ra0> 1 & ReL0 <= 1 & Rel0 >  1) = 6;
    % Ics(Ra0<=1 & ReL0 >  1 & Rel0 >  1) = 7;
    % Ics(Ra0> 1 & ReL0 >  1 & Rel0 >  1) = 8;

    % Ics(Ra0<=1 & ReL0 <= 1 & Rel0 <= 1) = 1;
    % Ics(Ra0> 1 & ReL0 <= 1 & Rel0 <= 1) = 2;
    % Ics(Ra0<=1 & ReL0 >  1 & Rel0 <= 1) = 3;
    % Ics(Ra0> 1 & ReL0 >  1 & Rel0 <= 1) = 4;
    % Ics(Ra0<=1 & ReL0 <= 1 & Rel0 >  1) = 5;
    % Ics(Ra0> 1 & ReL0 <= 1 & Rel0 >  1) = 6;
    % Ics(Ra0<=1 & ReL0 >  1 & Rel0 >  1) = 7;
    % Ics(Ra0> 1 & ReL0 >  1 & Rel0 >  1) = 8;

    set(0,'CurrentFigure',fig4); clf;
    imagesc(log10(D0(1,:)),log10(d0(:,1)),Ics); axis xy tight; box on; hold on; colormap(colmap); clim([1 4])

    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','w','LineStyle','-','LineWidth',2);
    line(log10([1e0,1e6]),log10([1e-2,1e-2]),'Color','k','LineStyle',':','LineWidth',2);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','w','LineStyle','-','LineWidth',2);
    line(log10([1e1,1e1]),log10([1e-4,1e-1]),'Color','k','LineStyle',':','LineWidth',2);
    for iD=0:1:2
        for id=-4:1:-2
            plot(iD,id,'ko','LineWidth',2,'MarkerFaceColor','w','MarkerSize',10)
        end
    end
    % plot(0.5,-1.5,'ko','LineWidth',1,'MarkerFaceColor','w')
    plot([6 6 6],[-3 -2 -1],'ko','LineWidth',2,'MarkerFaceColor','w','MarkerSize',10)
    plot(6,-2,'wo','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10)
    plot(log10(1e1),log10(1e-2),'wo','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10)

    set(gca,'TickLabelInterpreter','latex','FontSize',13)
    xlabel('log$_{10}$ Layer depth $D_0$ [m]','Interpreter','latex','FontSize',15);
    ylabel('log$_{10}$ Crystal size $d_0$ [m]','Interpreter','latex','FontSize',15);
    title('Flow regimes','Interpreter','latex','FontSize',15);
    drawnow;

    name = ['../out/scaling/maps_eta',num2str(log10(eta0))];
    print(fig3,name,'-dpng','-r300','-image');
    name = ['../out/scaling/regimes_eta',num2str(log10(eta0))];
    print(fig4,name,'-dpng','-r300','-image');

end

figure(1)
subplot(2,2,1)
line([1e-2,1e-2],[1e-8,1e2],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e-4,1e-1]); ylim([1e-8,1e2]); yticks([1e-8 1e-6 1e-4 1e-2 1e0 1e2])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Settling speed $w_0$ [m/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(a)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p11,p31,p51],{'laminar','turbulent','general'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,2)
line([1e1,1e1],[1e-8,1e2],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e0,1e6]); ylim([1e-8,1e2]); yticks([1e-8,1e-6,1e-4,1e-2,1e0,1e2]);
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Convective speed $W_0$ [m/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(b)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p1(1),p1(2),p1(3),p1(4)],{'$\eta_0=10^{-1}$ Pas','$\eta_0=10^{1}$ Pas','$\eta_0=10^{3}$ Pas','$\eta_0=10^{5}$ Pas'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,3)
line([1e-2,1e-2],[1e-9,1e6],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e-4,1e-1]); ylim([1e-9,1e6]); yticks([1e-9,1e-6,1e-3,1e0,1e3,1e6]);
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Dimensionless numbers [1]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(c)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p53,p54],{'Rc','Re$_d$'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,4)
line([1e1,1e1],[1e-2,1e8],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e0,1e6]); ylim([1e-2,1e8]); yticks([1e-4 1e-2 1e0 1e2 1e4 1e6 1e8])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Dimensionless numbers [1]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(d)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p55,p56],{'Ra','Re$_D$'},'Interpreter','latex','FontSize',11,'Location','southeast');


figure(2)
subplot(2,2,1)
line([1e-2,1e-2],[1e-12,1e0],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e-4,1e-1]); ylim([1e-12,1e0]); yticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Diffusivity [m$^2$/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(a)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p21,p41,p61],{'laminar','turbulent','general'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,2)
line([1e1,1e1],[1e-12,1e4],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e0,1e6]); ylim([1e-12,1e4]); yticks([1e-12 1e-8 1e-4 1e0 1e4])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Diffusivity [m$^2$/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(b)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p2(1),p2(2),p2(3),p2(4)],{'$\eta_0=10^{-1}$ Pas','$\eta_0=10^{1}$ Pas','$\eta_0=10^{3}$ Pas','$\eta_0=10^{5}$ Pas'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,3)
line([1e-2,1e-2],[1e-12,1e0],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e-4,1e-1]); ylim([1e-12,1e0]); yticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal size $d_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Noise amplitude $\xi_{s,0}$ [m$^2$/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(c)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p65,p66],{'$\kappa_{s,0}$, $\xi_{s,0}$','$\kappa_{x,0}$, $\xi_{x,0}$'},'Interpreter','latex','FontSize',11,'Location','southeast');

subplot(2,2,4)
line([1e1,1e1],[1e-10,1e2],'Color','k','LineStyle',':','LineWidth',1);
xlim([1e0,1e6]); ylim([1e-10,1e2]); yticks([1e-10 1e-8 1e-6 1e-4 1e-2 1e0 1e2])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Layer depth $D_0$ [m]','Interpreter','latex','FontSize',13);
ylabel('Noise amplitude [m/s]','Interpreter','latex','FontSize',13);
text(0.02,0.91,'\textbf{(d)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p67,p68],{'$\kappa_{e,0}$, $\xi_{e,0}$','$\kappa_{x,0}$, $\xi_{x,0}$'},'Interpreter','latex','FontSize',11,'Location','southeast');

name = '../out/scaling/lines1';
print(fig1,name,'-dpng','-r300','-image');
name = '../out/scaling/lines2';
print(fig2,name,'-dpng','-r300','-image');

