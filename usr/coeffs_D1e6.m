clear; close all;

fig1 = figure(1); clf; set(gcf,'Units','centimeters','Position',[6,6,20,15]);
% fig2 = figure(2); clf; set(gcf,'Units','centimeters','Position',[8,8,20,15]);
% fig3 = figure(3); clf; set(gcf,'Units','centimeters','Position',[10,10,20,15]);
% fig4 = figure(4); clf; set(gcf,'Units','centimeters','Position',[12,12,18,15]);
load ../src/colmap/vik.mat; colmap = vik;
if ~isfolder('../out/scaling'); mkdir('../out/scaling'); end

par_default;

D     = 1e6;
D0    = D/10;
L0    = D/100;
d0    = 0.1;
l0    = 10*d0;

chi0  = linspace(0.001,0.7-0.001,699).';
mu0   = 1-chi0;

Dchi0 = chi0./10;

etait = 0;
for eta0=10.^linspace(-1,5,4)
    etait = etait+1;

    etam0 = eta0;
    rho0  = chi0.*rhox0 + mu0.*rhom0;
    Drho0 = rhox0-rho0;

    lnshd = (6-(log10(eta0)+1))/8;
    blk   = [0 0 0]; wht = [1 1 1];
    red   = [0.8 0 0.1]; blu = [0.1 0 0.8];
    prp   = [0.5 0.1 0.5]; grn = [0.2 0.6 0.3];

    % get coefficient contrasts
    kv = [etax0;etam0];
    Mv = [etax0;etam0].'./[etax0;etam0];

    % get permission weights
    ff = permute(cat(3,chi0,mu0 ),[3,1,2]);
    FF = permute(repmat(ff,1,1,1,2),[4,1,2,3]);
    Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
    Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

    % get momentum flux and transfer coefficients
    thtv = squeeze(prod(Mv.^Xf,2));
    etaf = (kv.*thtv);

    % get effective viscosity
    digits(64);
    etamix = vpa(squeeze(sum(ff.*etaf,1))).';

    % get laminar convection and settling speeds
    W0l =  Dchi0.*Drho0.*g0.*D0^2./etamix;  % laminar convection speed
    w0l =         Drho0.*g0.*d0^2./etamix;  % laminar settling speed

    w0t  =  sqrt(       Drho0.*g0.*d0^2./(l0  .*rho0));  % terminal turbulent settling speed
    W0t  =  sqrt(Dchi0.*Drho0.*g0.*D0^3./(L0^2.*rho0));  % terminal turbulent convective speed
    W0i  =  sqrt(Dchi0.*Drho0.*g0.*D   ./(      rho0));  % inertially limited convective speed
    Ri0  =  W0i./(W0t);

    W0 = W0l;
    w0 = w0l;
    tol    = 1e-5;
    res    = 1;
    while res>tol
        W0prv = W0;
        w0prv = w0;

        ReL0  = vpa(W0.*L0./(etamix./rho0));    % convective Reynolds No at L0, eta0
        Rel0  = vpa(w0.*l0./(etamix./rho0));    % settling Reynolds No at l0, eta0

        fReL0 = vpa((1-exp(-ReL0)));        % Re-dependent ramp factor
        fRel0 = vpa((1-exp(-Rel0)));        % Re-dependent ramp factor

        % general convective speed
        W0    = ((sqrt(4./Ri0.^2.*Dchi0.*Drho0.*g0.*rho0.*fReL0.*L0.^2.*D0 + etamix.^2) - etamix).*D0./(2.*fReL0.*L0.^2.*rho0./Ri0.^2));

        % general settling speed
        w0    = ((sqrt(4               .*Drho0.*g0.*rho0.*fRel0.*l0.*d0.^2 + etamix.^2) - etamix)    ./(2.*fRel0.*l0   .*rho0        ));

        res   = (norm(W0(~isnan(W0))-W0prv(~isnan(W0)))./norm(W0(~isnan(W0))) + norm(w0(~isnan(w0))-w0prv(~isnan(w0)))./norm(w0(~isnan(w0))));  % residual
    end

    % diffusivities
    ke0   =  W0./D0.*L0^2;
    ks0   =  w0.*l0;
    kx0   =  ks0 + fReL0.*ke0;

    etae0 =  fReL0.*ke0.*rho0;
    etat0 =  fRel0.*ks0.*rho0;

    set(0,'CurrentFigure',fig1)
    subplot(2,2,1)
    p11 = semilogy(chi0,etamix,'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p12 = semilogy(chi0,etae0,'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p13 = semilogy(chi0,etamix + etae0,'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,2)
    p14 = semilogy(chi0,d0^2./etamix,'--','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p15 = semilogy(chi0,d0^2./etat0,'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p16 = semilogy(chi0,d0^2./(etamix + etat0),'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,3)
    p17 = semilogy(chi0,ke0,'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p18 = semilogy(chi0,kx0,'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    subplot(2,2,4)
    p19 = semilogy(chi0,ks0,'-.','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
    p20 = semilogy(chi0,kx0,'-','LineWidth',1.5,'Color',lnshd.*wht + (1-lnshd).*blk); axis tight; box on; hold on
end


figure(1)
subplot(2,2,1)
xlim([0,0.7]); ylim([1e-2,1e12]); yticks([1e-2 1e2 1e6 1e10])
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal fraction [vol]','Interpreter','latex','FontSize',13);
ylabel('Viscosity [Pas]','Interpreter','latex','FontSize',13);
text(0.02,0.07,'\textbf{(a)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p11,p12,p13],{'$\bar{\eta}$','$\eta_e$','$\eta$'},'Interpreter','latex','FontSize',11,'Location','south');

subplot(2,2,2)
xlim([0,0.7]); ylim([1e-12,1e-0]); yticks([1e-12,1e-8,1e-4,1e0]);
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal fraction [vol]','Interpreter','latex','FontSize',13);
ylabel('Settling coeff. [m$^2$/Pas]','Interpreter','latex','FontSize',13);
text(0.02,0.07,'\textbf{(b)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p14,p15,p16],{'$d_0^2\,/\,\bar{\eta}$','$d_0^2\,/\,\eta_s$','$d_0^2\,/\,\eta$'},'Interpreter','latex','FontSize',11,'Location','south');

subplot(2,2,3)
xlim([0,0.7]); ylim([1e2,1e6]); yticks([1e2,1e3,1e4,1e5,1e6]);
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal fraction [vol]','Interpreter','latex','FontSize',13);
ylabel('Diffusivity [m$^2$/s]','Interpreter','latex','FontSize',13);
text(0.02,0.07,'\textbf{(c)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p17,p18],{'$\kappa_e$','$\kappa_x$'},'Interpreter','latex','FontSize',11,'Location','south');

subplot(2,2,4)
xlim([0,0.7]); ylim([1e-9,1e6]); yticks([1e-9,1e-6,1e-3,1e0,1e3,1e6]);
set(gca,'TickLabelInterpreter','latex','FontSize',11)
xlabel('Crystal fraction [vol]','Interpreter','latex','FontSize',13);
ylabel('Diffusivity [m$^2$/s]','Interpreter','latex','FontSize',13);
text(0.02,0.07,'\textbf{(d)}','Interpreter','latex','FontSize',13,'Units','normalized')
legend([p19,p20],{'$\kappa_s$','$\kappa_x$'},'Interpreter','latex','FontSize',11,'Location','south');
