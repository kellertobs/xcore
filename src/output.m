%*****  CREATE AND SAVE OUTPUT FIGURE AND DATA FRAMES  ********************

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',1.5};
XR = {'XDir','reverse'};
MS = {'MarkerSize',8};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

% set axis and border dimensions
if Nx>1
    axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
else
    axh = 6.0; axw = 6.0;
end
ahs = 0.60; avs = 0.80;
axb = 1.20; axt = 1.50;
axl = 1.50; axr = 0.60;

% initialize figures and axes
if ~exist('fh1','var'); fh1 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh1); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh1,UN{:},'Position',[1 1 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off','Resize','off');
ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(13) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(14) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(15) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(16) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

if ~exist('fh2','var'); fh2 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh2); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh2,UN{:},'Position',[3 3 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(23) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(24) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(25) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(26) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

if ~exist('fh3','var'); fh3 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh3); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh3,UN{:},'Position',[5 5 fw fh]);
set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(33) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(34) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

if ~exist('fh4','var'); fh4 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh4); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh4,UN{:},'Position',[7 7 fw fh]);
set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh4,'Color','w','InvertHardcopy','off','Resize','off');
ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(43) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(44) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

% plot velocity-pressure solution in Fig. 1
set(0,'CurrentFigure',fh1)
set(fh1,'CurrentAxes',ax(11));
imagesc(Xsc,Zsc,-W(:,2:end-1)/Wsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',Wun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',sun,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(12));
imagesc(Xsc,Zsc, V./Wsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',Wun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh1,'CurrentAxes',ax(13));
imagesc(Xsc,Zsc, P(2:end-1,2:end-1)/psc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [',pun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(14));
imagesc(Xsc,Zsc,rho/rsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [',dun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(15));
imagesc(Xsc,Zsc,log10(eta/(esc+eesc))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta$ [',eun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(16));
imagesc(Xsc,Zsc,MFS/MFSsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \bar{\rho} \mathbf{v}$ [',MFSun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot density, rheology, and segregation speeds in Fig. 2
set(0,'CurrentFigure',fh2)
set(fh2,'CurrentAxes',ax(21));
imagesc(Xsc,Zsc,x/xsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$x$ [',xun,']'],TX{:},FS{:});  ylabel(['Depth [',sun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
set(fh2,'CurrentAxes',ax(22));
imagesc(Xsc,Zsc,Gx./Gsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x$ [',Gun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh2,'CurrentAxes',ax(23));
imagesc(Xsc,Zsc,log10(kx/kxsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_x$ [',kun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh2,'CurrentAxes',ax(24));
imagesc(Xsc,Zsc,-wx(:,2:end-1)/wxsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_x$ [',wun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:});  xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh2,'CurrentAxes',ax(25));
imagesc(Xsc,Zsc,-wm(:,2:end-1)/wmsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_m$ [',wun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:}); 
set(fh2,'CurrentAxes',ax(26));
imagesc(Xsc,Zsc,log10(etas/(esc+essc))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta_s$ [',eun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh3)
set(fh3,'CurrentAxes',ax(31));
imagesc(Xsc,Zsc,log10(Ra/Rasc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Ra [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',sun,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(32));
imagesc(Xsc,Zsc,log10(Rc/Rcsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Rc [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.05,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh3,'CurrentAxes',ax(33));
imagesc(Xsc,Zsc,log10(ReD/ReDsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Re$_D$ [1]'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(34));
imagesc(Xsc,Zsc,log10(Red/Redsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Re$_d$ [1]'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh4)
set(fh4,'CurrentAxes',ax(41));
imagesc(Xsc,Zsc,log10(ks/kssc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_s$ [',kun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',sun,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(42));
imagesc(Xsc,Zsc,xix/xixsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_x$ [',Wun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.05,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh4,'CurrentAxes',ax(43));
imagesc(Xsc,Zsc,log10(ke/kesc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_e$ [',kun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(44));
imagesc(Xsc,Zsc,xie/xiesc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_e$ [',Wun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot model history
if ~exist('fh14','var'); fh14 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh14); clf;
end

subplot(3,1,1)
plot(HST.time/tsc,HST.x(:,2)/xsc,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/tsc,(HST.x(:,2)+HST.x(:,4))/xsc,'k--',LW{1},1);
plot(HST.time/tsc,HST.x(:,[1,3])/xsc,'k:',LW{:});
% plot(HST.time/tsc,HST.x_tavg(:,2)/xsc,'k-',LW{1},1,'Color',[1 1 1]*0.5);
% plot(HST.time/tsc,HST.x_tavg(:,[1,3])/xsc,'k:',LW{1},1,'Color',[1 1 1]*0.5);

plot(HST.time(end)/tsc,x0/xsc,'ks' ,LW{:},MS{:});
plot(HST.time(end)/tsc,xeq/xsc,'ko',LW{:},MS{:});
markerlabels([xeq/xsc,x0/xsc],{'$x_\mathrm{eq}$','$x_0$'},FS,'right');

set(gca,TL{:},FS{:},'Xticklabels',[])
legend('mean','mean + std','min/max',TX{:},FS{:},'Location','northwest')
title(['Crystallinity [',xun,']'],TX{:},FS{1},15);

subplot(3,1,2)
semilogy(HST.time/tsc,HST.V(:,2)/Wsc,'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.vx(:,2)/whsc,'k--',LW{:});
semilogy(HST.time/tsc,HST.xie(:,2)/xiesc,'-',LW{:},'Color',[1 1 1]/2);
semilogy(HST.time/tsc,HST.xix(:,2)/xixsc,'--',LW{:},'Color',[1 1 1]/2);

% semilogy(HST.time/tsc,HST.V_tavg(:,2)/Wsc,'-' ,LW{1},1,'Color',[1 1 1]/2);
% semilogy(HST.time/tsc,HST.vx_tavg(:,2)/whsc,'--',LW{1},1,'Color',[1 1 1]/2);
% semilogy(HST.time/tsc,HST.xie_tavg(:,2)/xiesc,'-',LW{1},1,'Color',[1 1 1]/1.5);
% semilogy(HST.time/tsc,HST.xix_tavg(:,2)/xixsc,'--',LW{1},1,'Color',[1 1 1]/1.5);

semilogy(HST.time(end)/tsc,W0/Wsc,'ko' ,LW{:},MS{:});
semilogy(HST.time(end)/tsc,w0/whsc,'kv',LW{:},MS{:});
semilogy(HST.time(end)/tsc,xie0/xiesc,'o',LW{:},'Color',[1 1 1]/2,MS{:});
semilogy(HST.time(end)/tsc,xix0/xixsc,'v',LW{:},'Color',[1 1 1]/2,MS{:});
markerlabels([W0/Wsc,w0/whsc,xie0/xiesc,xix0/xixsc],{'$W_0$','$w_{0}$','$\xi_{e,0}$','$\xi_{x,0}$'},FS,'right');

yticks = 10.^(-4:1:4);
yticklabels = {'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$','$10^{4}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
set(gca,TL{:},FS{:},'Xticklabel',[]);
legend('convection','settling','eddy noise','xtal noise',TX{:},FS{:},'Location','northwest')
title(['Flow speeds [',Wun,']'],TX{:},FS{1},15);

subplot(3,1,3)
semilogy(HST.time/tsc,HST.Ra(:,2)/Rasc,'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.ReD(:,2)/ReDsc,'k-.',LW{:});
semilogy(HST.time/tsc,HST.Red(:,2)/Redsc,'k:' ,LW{:});
semilogy(HST.time/tsc,HST.Rc (:,2)/Rcsc,'k--',LW{:});

% semilogy(HST.time/tsc,HST.Ra_tavg(:,2)/Rasc,'k-' ,LW{1},1,'Color',[1 1 1]*0.5);
% semilogy(HST.time/tsc,HST.ReD_tavg(:,2)/ReDsc,'k-.',LW{1},1,'Color',[1 1 1]*0.5);
% semilogy(HST.time/tsc,HST.Red_tavg(:,2)/Redsc,'k:' ,LW{1},1,'Color',[1 1 1]*0.5);
% semilogy(HST.time/tsc,HST.Rc_tavg (:,2)/Rcsc,'k--',LW{1},1,'Color',[1 1 1]*0.5);

semilogy(HST.time(end)/tsc,Ra0/Rasc,'ko' ,LW{:},MS{:});
semilogy(HST.time(end)/tsc,ReD0/ReDsc,'ks',LW{:},MS{:});
semilogy(HST.time(end)/tsc,Red0/Redsc,'kd',LW{:},MS{:});
semilogy(HST.time(end)/tsc,Rc0/Rcsc,'kv',LW{:},MS{:});
markerlabels([Ra0/Rasc,ReD0/ReDsc,Red0/Redsc,Rc0/Rcsc],{'Ra$_0$','Re$_{D,0}$','Re$_{d,0}$','Rc$_0$'},FS,'right');

yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
legend('Ra','Re$_D$','Re$_d$','Rc',TX{:},FS{:},'Location','northwest')
title(['Dimensionless numbers'],TX{:},FS{1},15);
xlabel(['Time [',tun,']'],TX{:},FS{1},15);

% plot model history
if ~exist('fh15','var'); fh15 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh15); clf;
end

subplot(2,1,1)
yyaxis left
semilogy(HST.time/tsc,HST.ke(:,2)/kesc,'-' ,LW{:},'Color',lncls(1,:)); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.ke(:,[1,3])/kesc,':',LW{:},'Color',lncls(1,:));

semilogy(HST.time(end)/tsc,ke0/kesc,'ko' ,LW{:},MS{:});
% markerlabels([ke0/kesc],{'k$_{e,0}$'},FS,'left');

set(gca,TL{:},FS{:},'YColor',lncls(1,:),'YScale','log');
yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');

yyaxis right
semilogy(HST.time/tsc,HST.etae(:,2)/eesc,'--' ,LW{:},'Color',lncls(2,:)); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.etae(:,[1,3])/eesc,':',LW{:},'Color',lncls(2,:));
set(gca,TL{:},FS{:},'YColor',lncls(2,:),'YScale','log');

semilogy(HST.time(end)/tsc,etae0/eesc,'ko' ,LW{:},MS{:});
markerlabels([etae0/eesc,etae0/eesc],{'k$_{e,0}$','$\eta_{e,0}$'},FS,'left');

yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
legend('mean','min/max',TX{:},FS{:},'Location','southeast')
title(['Eddy diffusivity [',kun,']  $|$  Eddy viscosity [',eun,']'],TX{:},FS{1},15);

subplot(2,1,2)
yyaxis left
semilogy(HST.time/tsc,HST.ks(:,2)/kssc,'-' ,LW{:},'Color',lncls(1,:)); hold on; box on;
semilogy(HST.time/tsc,HST.ks(:,[1,3])/kssc,':',LW{:},'Color',lncls(1,:));

semilogy(HST.time(end)/tsc,ks0/kssc,'ko' ,LW{:},MS{:});
% markerlabels([ks0/kssc],{'k$_{s,0}$'},FS,'left');

set(gca,TL{:},FS{:},'YColor',lncls(1,:),'YScale','log');

yyaxis right
semilogy(HST.time/tsc,HST.etat(:,2)/essc,'--' ,LW{:},'Color',lncls(2,:)); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.etat(:,[1,3])/essc,':',LW{:},'Color',lncls(2,:));

semilogy(HST.time(end)/tsc,etas0/essc,'ko' ,LW{:},MS{:});
markerlabels([etas0/essc,etas0/essc],{'k$_{s,0}$','$\eta_{s,0}$'},FS,'left');

set(gca,TL{:},FS{:},'YColor',lncls(2,:),'YScale','log');
legend('mean','min/max',TX{:},FS{:},'Location','southeast')
title(['Settling diffusivity [',kun,']  $|$  Settling viscosity [',eun,']'],TX{:},FS{1},15);
xlabel(['Time [',tun,']'],TX{:},FS{1},15);

if plot_cv
if ~exist('fh16','var'); fh16 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh16); clf;
end
plot(HST.time/tsc,HST.EB,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/tsc,HST.EM,'k-.',LW{:});
plot(HST.time/tsc,HST.EX,'k--',LW{:});
set(gca,TL{:},FS{:})
legend('xtal','melt','mixt',TX{:},FS{:},'Location','northwest')
ylabel('Rel. error [1]',TX{:},FS{1},15);
xlabel(['Time [',tun,']'],TX{:},FS{1},15);
title('Global Conservation Error',TX{:},FS{1},18)
end

drawnow

% save output to file
if save_op && ~restart
    name = [outdir,'/',runID,'/',runID,'_cnv_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_ndn_',num2str(floor(step/nop))];
    print(fh3,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_dfn_',num2str(floor(step/nop))];
    print(fh4,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_hnd'];
    print(fh14,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_hdf'];
    print(fh15,name,'-dpng','-r300','-image');

    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','dt','time','step','MFS','wx','wm','psie','psiex','psis','xisw','xisu','xiew','xieu','xiexw','xiexu');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','dt','time','step','MFS','wx','wm','psie','psiex','psis','xisw','xisu','xiew','xieu','xiexw','xiexu');
    name = [outdir,'/',runID,'/',runID,'_HST'];
    save(name,'HST');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
