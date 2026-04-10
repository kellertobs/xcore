%*****  CREATE AND SAVE OUTPUT FIGURE AND DATA FRAMES  ********************

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = colororder;
LW = {'LineWidth',1.5};
XR = {'XDir','reverse'};
MS = {'MarkerSize',6};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

if ndm_op
    tsc = t0;  tun = '1';
else
    % update time scale and units
    if time < 1e3
        tsc = 1;
        tun = 's';
    elseif time>= 1e3 && time < 1e3*hr
        tsc = hr;
        tun = 'hr';
    elseif time >= 1e3*hr && time < 1e2*yr
        tsc = yr;
        tun = 'yr';
    elseif time >= 1e2*yr
        tsc = 1e3*yr;
        tun = 'kyr';
    end
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
fw = axl + 2*axw + 1*ahs + axr;
set(fh2,UN{:},'Position',[3 3 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(23) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(24) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

if ~exist('fh3','var'); fh3 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh3); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh3,UN{:},'Position',[5 5 fw fh]);
set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(33) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(34) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(35) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(36) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

if ~exist('fh4','var'); fh4 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh4); clf;
end
colormap(colmap);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh4,UN{:},'Position',[7 7 fw fh]);
set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh4,'Color','w','InvertHardcopy','off','Resize','off');
ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(43) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(44) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(45) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(46) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

% plot velocity-pressure solution in Fig. 1
set(0,'CurrentFigure',fh1)
set(fh1,'CurrentAxes',ax(11));
imagesc(Xsc,Zsc,-W(:,2:end-1)./Wsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',Wun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',sun,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(12));
imagesc(Xsc,Zsc,U(2:end-1,:)./Wsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',Wun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh1,'CurrentAxes',ax(13));
imagesc(Xsc,Zsc, P(2:end-1,2:end-1)/psc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [',pun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(14));
imagesc(Xsc,Zsc,rho/rsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [',dun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(15));
imagesc(Xsc,Zsc,log10(eta/(esc+eesc))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta$ [',eun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(16));
imagesc(Xsc,Zsc,MFS/MFSsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \rho \mathbf{v}$ [',MFSun,']'],TX{:},FS{:}); xlabel(['Width [',sun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot density, rheology, and segregation speeds in Fig. 2
set(0,'CurrentFigure',fh2)
set(fh2,'CurrentAxes',ax(21));
if max(x(:))>10*chi0
    imagesc(Xsc,Zsc,log10(x/xsc)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $x$ [',xun,']'],TX{:},FS{:});  ylabel(['Depth [',sun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
else
    imagesc(Xsc,Zsc,x/xsc); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$x$ [',xun,']'],TX{:},FS{:});  ylabel(['Depth [',sun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
end
set(fh2,'CurrentAxes',ax(22));
imagesc(Xsc,Zsc,log10(etas/(esc+etsc))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta_s$ [',eun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.1,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh2,'CurrentAxes',ax(23));
imagesc(Xsc,Zsc,-wx(2:end-1,2:end-1)/wxsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_x$ [',wun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:});  xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh2,'CurrentAxes',ax(24));
imagesc(Xsc,Zsc,-wm(2:end-1,2:end-1)/wmsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_m$ [',wun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh3)
set(fh3,'CurrentAxes',ax(31));
imagesc(Xsc,Zsc,log10(ks/kssc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_s$ [',xun,']'],TX{:},FS{:});  ylabel(['Depth [',sun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
set(fh3,'CurrentAxes',ax(32));
imagesc(Xsc,Zsc,log10(kx/kxsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_x$ [',Gun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh3,'CurrentAxes',ax(33));
imagesc(Xsc,Zsc,log10(ke/kesc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_e$ [',kun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh3,'CurrentAxes',ax(34));
imagesc(Xsc,Zsc,xie/xiesc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_e$ [',xieun,']'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:});  xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(35));
imagesc(Xsc,Zsc,xix/xixsc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_x$ [',xieun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(36));
imagesc(Xsc,Zsc,xis/xissc); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_s$ [',xisun,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh4)
set(fh4,'CurrentAxes',ax(41));
imagesc(Xsc,Zsc,log10(ReD/ReDsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Re$_D$ [1]'],TX{:},FS{:});  ylabel(['Depth [',sun,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
set(fh4,'CurrentAxes',ax(42));
imagesc(Xsc,Zsc,log10(Ra/Rasc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Ra [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.05+Nx/Nz/10,['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh4,'CurrentAxes',ax(43));
imagesc(Xsc,Zsc,log10(Rc/Rcsc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Rc [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh4,'CurrentAxes',ax(44));
imagesc(Xsc,Zsc,log10(Noe/Noesc)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ No$_e$ [1]'],TX{:},FS{:}); ylabel(['Depth [',sun,']'],TX{:},FS{:});  xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(45));
imagesc(Xsc,Zsc,log10(Nox/Noxsc)); axis ij equal tight; box on; cb = colorbar;
clim([-4,max(-3, max(log10(Nox(:)/Noxsc)))])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ No$_x$ [1]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(46));
imagesc(Xsc,Zsc,log10(Nos/Nossc)); axis ij equal tight; box on; cb = colorbar;
clim([-4,max(-3, max(log10(Nos(:)/Nossc)))])
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ No$_s$ [1]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',sun,']'],TX{:},FS{:});

% plot 1D horizontal average profiles
if ~exist('fh13','var'); fh13 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh13); clf;
end

sgtitle(['time = ',num2str(time/tsc,3),' [',tun,']'],TX{:},'FontSize',14,'Color','k');
subplot(1,4,1)
plot(mean(x./xsc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(1,:)); axis ij tight; box on; hold on
plot(min(x./xsc,[],2),Zc./ssc,'k:','LineWidth',0.75,'Color',CL(1,:));
plot(max(x./xsc,[],2),Zc./ssc,'k:','LineWidth',0.75,'Color',CL(1,:));

set(gca,TL{:},FS{:})
ylabel(['Depth [',sun,']'],TX{:},FS{:});
legend('mean','min','max',TX{:},TS{:},'Location','east')
title(['Crystallinity [',xun,']'],TX{:},FS{1},14);

subplot(1,4,2)
plot(rms(vx./wxpsc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(1,:)); axis ij tight; box on; hold on
plot(rms(vm./wmpsc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(2,:));
plot(rms(V ./Wpsc ,2),Zc./ssc,'LineWidth',1.5,'Color',CL(3,:));
plot(min(vx./wxpsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
plot(max(vx./wxpsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
plot(min(vm./wmpsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
plot(max(vm./wmpsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
plot(min(V ./Wpsc ,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));
plot(max(V ./Wpsc ,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));

set(gca,TL{:},FS{:})
legend('$|\Delta \mathbf{v}_x|$','$|\Delta \mathbf{v}_m|$','$|\mathbf{v}|$', TX{:},TS{:},'Location','east')
title(['Velocity [',wpun,']'],TX{:},FS{1},14);

subplot(1,4,3)
semilogx(geomean(ke./kesc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(1,:)); axis ij tight; box on; hold on
semilogx(geomean(ks./kssc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(2,:));
semilogx(geomean(kx./kxsc,2),Zc./ssc,'LineWidth',1.5,'Color',CL(3,:));
semilogx(min(ke./kesc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
semilogx(max(ke./kesc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
semilogx(min(ks./kssc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
semilogx(max(ks./kssc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
semilogx(min(kx./kxsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));
semilogx(max(kx./kxsc,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));

set(gca,TL{:},FS{:})
legend('$\kappa_e$','$\kappa_s$','$\kappa_x$',TX{:},TS{:},'Location','east')
title(['Diffusivity [',kun,']'],TX{:},FS{1},14);

subplot(1,4,4)
semilogx(geomean(eta   ./(esc+eesc),2),Zc./ssc,'LineWidth',1.5,'Color',CL(1,:)); axis ij tight; box on; hold on
semilogx(geomean(etas  ./(esc+etsc),2),Zc./ssc,'LineWidth',1.5,'Color',CL(2,:));
semilogx(geomean(etamix./ esc      ,2),Zc./ssc,'LineWidth',1.5,'Color',CL(3,:));
semilogx(min(eta   ./(esc+eesc),[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
semilogx(max(eta   ./(esc+eesc),[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(1,:));
semilogx(min(etas  ./(esc+etsc),[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
semilogx(max(etas  ./(esc+etsc),[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(2,:));
semilogx(min(etamix./ esc      ,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));
semilogx(max(etamix./ esc      ,[],2),Zc./ssc,':','LineWidth',0.75,'Color',CL(3,:));

set(gca,TL{:},FS{:})
legend('$\eta$','$\eta_s$','$\bar{\eta}$',TX{:},TS{:},'Location','east')
title(['Viscosity [',eun,']'],TX{:},FS{1},14);

% plot model history
if ~exist('fh14','var'); fh14 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh14); clf;
end

subplot(3,1,1)
plot(HST.time/tsc,HST.x(:,2)/xsc,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/tsc,(HST.x(:,2)+HST.x(:,4))/xsc,'k--',LW{1},1);
plot(HST.time/tsc,HST.x(:,[1,3])/xsc,'k:',LW{:});

plot(HST.time(end)/tsc,Da/xsc,'ks' ,LW{:},MS{:});
markerlabels([Da/xsc],{'$x_0$'},FS,'right');

set(gca,TL{:},FS{:},'Xticklabels',[])
legend('mean','mean + std','min/max',TX{:},TS{:},'Location','northwest')
title(['Crystallinity [',xun,']'],TX{:},FS{1},14);

subplot(3,1,2)
semilogy(HST.time/tsc,HST.V(:,2)/Wsc,'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.vx(:,2)/whsc,'k--',LW{:});
semilogy(HST.time/tsc,HST.xie(:,2)/xiesc,'-',LW{:},'Color',[1 1 1]/2);
semilogy(HST.time/tsc,HST.xis(:,2)/xissc,'--',LW{:},'Color',[1 1 1]/2);
semilogy(HST.time/tsc,HST.xix(:,2)/xixsc,'-.',LW{:},'Color',[1 1 1]/2);

semilogy(HST.time(end)/tsc,W0/Wsc,'ko' ,LW{:},MS{:});
semilogy(HST.time(end)/tsc,w0/whsc,'kv',LW{:},MS{:});
semilogy(HST.time(end)/tsc,xie0/xiesc,'o',LW{:},'Color',[1 1 1]/2,MS{:});
semilogy(HST.time(end)/tsc,xis0/xissc,'v',LW{:},'Color',[1 1 1]/2,MS{:});
semilogy(HST.time(end)/tsc,xix0/xixsc,'d',LW{:},'Color',[1 1 1]/2,MS{:});
markerlabels([W0/Wsc,w0/whsc,xie0/xiesc,xis0/xissc,xix0/xixsc],{'$W_0$','$w_{0}$','$\xi_{e,0}$','$\xi_{s,0}$','$\xi_{x,0}$'},FS,'right');

yticks = 10.^(-8:2:4);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
set(gca,TL{:},FS{:},'Xticklabel',[]);
legend('convection','settling','eddy noise','xtal noise','xtal-eddy noise',TX{:},TS{:},'Location','northwest')
title(['Flow speeds [',Wun,']'],TX{:},FS{1},14);

subplot(3,1,3)
semilogy(HST.time/tsc,HST.Ra(:,2)/Rasc,'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/tsc,HST.ReD(:,2)/ReDsc,'k-.',LW{:});
semilogy(HST.time/tsc,HST.Red(:,2)/Redsc,'k:' ,LW{:});
semilogy(HST.time/tsc,HST.Rc (:,2)/Rcsc,'k--',LW{:});

semilogy(HST.time(end)/tsc,Ra0/Rasc,'ko' ,LW{:},MS{:});
semilogy(HST.time(end)/tsc,ReD0/ReDsc,'ks',LW{:},MS{:});
semilogy(HST.time(end)/tsc,Red0/Redsc,'kd',LW{:},MS{:});
semilogy(HST.time(end)/tsc,Rc0/Rcsc,'kv',LW{:},MS{:});
markerlabels([Ra0/Rasc,ReD0/ReDsc,Red0/Redsc,Rc0/Rcsc],{'Ra$_0$','Re$_{D,0}$','Re$_{d,0}$','Rc$_0$'},FS,'right');

yticks = 10.^(-4:1:2);
yticklabels = {'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
legend('Ra','Re$_D$','Re$_d$','Rc',TX{:},TS{:},'Location','northwest')
title(['Dimensionless numbers'],TX{:},FS{1},14);
xlabel(['Time [',tun,']'],TX{:},FS{1},14);


if plot_cv
if ~exist('fh16','var'); fh16 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh16); clf;
end
plot(HST.time/tsc,HST.EB,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/tsc,HST.EM,'k-.',LW{:});
plot(HST.time/tsc,HST.EX,'k--',LW{:});
set(gca,TL{:},FS{:})
legend('xtal','melt','mixt',TX{:},TS{:},'Location','northwest')
ylabel('Rel. error [1]',TX{:},FS{1},14);
xlabel(['Time [',tun,']'],TX{:},FS{1},14);
title('Global Conservation Error',TX{:},FS{1},14)
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
    name = [outdir,'/',runID,'/',runID,'_prf_',num2str(floor(step/nop))];
    print(fh13,name,'-dpng','-r300','-image');
    name = [outdir,'/',runID,'/',runID,'_hnd'];
    print(fh14,name,'-dpng','-r300','-image');

    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','Noe','Nos','Nox','dt','time','step','MFS','wx','wm','psie','psix','psis','xisw','xisu','xiew','xieu','xixw','xixu');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','x','m','chi','mu','X','M','dXdt','drhodt','Gx','rho','eta','etas','ke','ks','kx','Ra','ReD','Red','Rc','Noe','Nos','Nox','dt','time','step','MFS','wx','wm','psie','psix','psis','xisw','xisu','xiew','xieu','xixw','xixu');
    name = [outdir,'/',runID,'/',runID,'_HST'];
    save(name,'HST');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
