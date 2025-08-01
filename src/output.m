%*****  CREATE AND SAVE OUTPUT FIGURE AND DATA FRAMES  ********************

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',2};
XR = {'XDir','reverse'};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

% adjust scales and units for intuitive visualisation
if time < 1e3
    TimeScale = 1;
    TimeUnits = 's';
elseif time>= 1e3 && time < 1e3*hr
    TimeScale = hr;
    TimeUnits = 'hr';
elseif time >= 1e3*hr && time < 1e2*yr
    TimeScale = yr;
    TimeUnits = 'yr';
elseif time >= 1e2*yr
    TimeScale = 1e3*yr;
    TimeUnits = 'kyr';
end
if D < 1e3
    SpaceScale = 1;
    SpaceUnits = 'm';
elseif D >= 1e3
    SpaceScale = 1e3;
    SpaceUnits = 'km';
end
if max([V(:);vx(:)]) < 1000/yr
    SpeedScale = 1/yr;
    SpeedUnits = 'm/yr';
elseif max([V(:);vx(:)]) >= 1000/yr && max([V(:);vx(:)]) < 1000/hr
    SpeedScale = 1/hr;
    SpeedUnits = 'm/hr';
elseif max([V(:);vx(:)]) >= 1000/hr
    SpeedScale = 1;
    SpeedUnits = 'm/s';
end
if max(abs(P(:))) < 1e2
    PresScale = 1;
    PresUnits = 'Pa';
elseif max(abs(P(:))) >= 1e3 && max(abs(P(:))) < 1e6
    PresScale = 1e3;
    PresUnits = 'kPa';
elseif max(abs(P(:))) >= 1e6 && max(abs(P(:))) < 1e9
    PresScale = 1e6;
    PresUnits = 'MPa';
else
    PresScale = 1e9;
    PresUnits = 'GPa';
end
Xsc = Xc./SpaceScale;
Zsc = Zc./SpaceScale;
Zsf = Zf./SpaceScale;

% set axis and border dimensions
if Nx>1
    axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
else
    axh = 6.0; axw = 6.0;
end
ahs = 0.60; avs = 0.80;
axb = 1.20; axt = 1.50;
axl = 1.50; axr = 0.50;

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
imagesc(Xsc,Zsc,-W(:,2:end-1)/SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(12));
imagesc(Xsc,Zsc, U(2:end-1,:)./SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh1,'CurrentAxes',ax(13));
imagesc(Xsc,Zsc, P(2:end-1,2:end-1)/PresScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [',PresUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(14));
imagesc(Xsc,Zsc,rho); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
set(fh1,'CurrentAxes',ax(15));
imagesc(Xsc,Zsc,log10(eta)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta$ [Pas]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh1,'CurrentAxes',ax(16));
imagesc(Xsc,Zsc,MFS*TimeScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \bar{\rho} \mathbf{v}$ [kg/m$^3$/',TimeUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot density, rheology, and segregation speeds in Fig. 2
set(0,'CurrentFigure',fh2)
set(fh2,'CurrentAxes',ax(21));
imagesc(Xsc,Zsc,min(2*xeq,x)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$x$ [wt\%]'],TX{:},FS{:});  ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
set(fh2,'CurrentAxes',ax(22));
imagesc(Xsc,Zsc,Gx./rho*TimeScale*100.*(chi>eps^0.5)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/',TimeUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(0.5,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh2,'CurrentAxes',ax(23));
imagesc(Xsc,Zsc,log10(kx)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_x$ [m$^2$/s]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
set(fh2,'CurrentAxes',ax(24));
imagesc(Xsc,Zsc,-wx(2:end-1,2:end-1)/SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_x$ [',SpeedUnits,']'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});  xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
set(fh2,'CurrentAxes',ax(25));
imagesc(Xsc,Zsc,-wm(2:end-1,2:end-1)/SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta W_m$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); 
set(fh2,'CurrentAxes',ax(26));
imagesc(Xsc,Zsc,log10(etas)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\eta_s$ [Pas]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh3)
set(fh3,'CurrentAxes',ax(31));
imagesc(Xsc,Zsc,log10(RaD)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Ra [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(32));
imagesc(Xsc,Zsc,log10(Rs)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Rs [1]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.05,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh3,'CurrentAxes',ax(33));
imagesc(Xsc,Zsc,log10(ReD)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Re$_D$ [1]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
set(fh3,'CurrentAxes',ax(34));
imagesc(Xsc,Zsc,log10(Red)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Re$_d$ [1]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot phase and eddy diffusivities, Ra and Re numbers in Fig. 3
set(0,'CurrentFigure',fh4)
set(fh4,'CurrentAxes',ax(41));
imagesc(Xsc,Zsc,log10(ks)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_s$ [m$^2$/s]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(42));
imagesc(Xsc,Zsc,xix/SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_x$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
text(-0.05,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
set(fh4,'CurrentAxes',ax(43));
imagesc(Xsc,Zsc,log10(ke)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ $\kappa_e$ [m$^2$/s]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
set(fh4,'CurrentAxes',ax(44));
imagesc(Xsc,Zsc,xie/SpeedScale); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\xi_e$ [',SpeedUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

% plot model history
if ~exist('fh14','var'); fh14 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh14); clf;
end

subplot(3,1,1)
plot(HST.time/TimeScale,HST.x(:,2)*100,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/TimeScale,HST.x(:,[1,3])*100,'k:',LW{:});
plot(HST.time/TimeScale,HST.x_tavg(:,2)*100,'k-',LW{1},1,'Color',[1 1 1]*0.5);
plot(HST.time/TimeScale,HST.x_tavg(:,[1,3])*100,'k:',LW{1},1,'Color',[1 1 1]*0.5);

plot(HST.time(end)/TimeScale,x0*100,'ks' ,LW{:});
plot(HST.time(end)/TimeScale,xeq*100,'ko',LW{:});

set(gca,TL{:},FS{:},'Xticklabels',[])
legend('mean','min/max',TX{:},FS{:},'Location','northwest')
title('Crystallinity [wt\%]',TX{:},FS{1},15);

subplot(3,1,2)
semilogy(HST.time/TimeScale,HST.V(:,2)/SpeedScale,'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/TimeScale,HST.vx(:,2)/SpeedScale,'k--',LW{:});
semilogy(HST.time/TimeScale,HST.xie(:,2)/SpeedScale,'-',LW{:},'Color',[1 1 1]/2);
semilogy(HST.time/TimeScale,HST.xix(:,2)/SpeedScale,'--',LW{:},'Color',[1 1 1]/2);

semilogy(HST.time/TimeScale,HST.V_tavg(:,2)/SpeedScale,'-' ,LW{1},1,'Color',[1 1 1]/2);
semilogy(HST.time/TimeScale,HST.vx_tavg(:,2)/SpeedScale,'--',LW{1},1,'Color',[1 1 1]/2);
semilogy(HST.time/TimeScale,HST.xie_tavg(:,2)/SpeedScale,'-',LW{1},1,'Color',[1 1 1]/1.5);
semilogy(HST.time/TimeScale,HST.xix_tavg(:,2)/SpeedScale,'--',LW{1},1,'Color',[1 1 1]/1.5);

semilogy(HST.time(end)/TimeScale,W0/SpeedScale,'ko' ,LW{:});
semilogy(HST.time(end)/TimeScale,w0/SpeedScale,'kv',LW{:});
semilogy(HST.time(end)/TimeScale,xie0/SpeedScale,'o',LW{:},'Color',[1 1 1]/2);
semilogy(HST.time(end)/TimeScale,xix0/SpeedScale,'v',LW{:},'Color',[1 1 1]/2);

yticks = 10.^(-4:1:4);
yticklabels = {'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$','$10^{4}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
set(gca,TL{:},FS{:},'Xticklabel',[]);
legend('convection','settling','eddy noise flux','xtal noise flux',TX{:},FS{:},'Location','northwest')
title(['Flow speeds [',SpeedUnits,']'],TX{:},FS{1},15);

subplot(3,1,3)
semilogy(HST.time/TimeScale,HST.RaD(:,2),'k-' ,LW{:}); hold on; axis tight; box on;
semilogy(HST.time/TimeScale,HST.ReD(:,2),'k-.',LW{:});
semilogy(HST.time/TimeScale,HST.Red(:,2),'k:' ,LW{:});
semilogy(HST.time/TimeScale,HST.Rs (:,2),'k--',LW{:});

semilogy(HST.time/TimeScale,HST.RaD_tavg(:,2),'k-' ,LW{1},1,'Color',[1 1 1]*0.5);
semilogy(HST.time/TimeScale,HST.ReD_tavg(:,2),'k-.',LW{1},1,'Color',[1 1 1]*0.5);
semilogy(HST.time/TimeScale,HST.Red_tavg(:,2),'k:' ,LW{1},1,'Color',[1 1 1]*0.5);
semilogy(HST.time/TimeScale,HST.Rs_tavg (:,2),'k--',LW{1},1,'Color',[1 1 1]*0.5);

semilogy(HST.time(end)/TimeScale,Ra0,'ko' ,LW{:});
semilogy(HST.time(end)/TimeScale,ReD0,'ks',LW{:});
semilogy(HST.time(end)/TimeScale,Red0,'kd',LW{:});
semilogy(HST.time(end)/TimeScale,Rs0,'kv',LW{:});

yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
legend('Ra','Re$_D$','Re$_d$','Rs',TX{:},FS{:},'Location','northwest')
title(['Dimensionless numbers'],TX{:},FS{1},15);
xlabel(['Time [',TimeUnits,']'],TX{:},FS{1},15);

% plot model history
if ~exist('fh15','var'); fh15 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh15); clf;
end

subplot(2,1,1)
yyaxis left
semilogy(HST.time/TimeScale,HST.ke(:,2),'-' ,LW{:},'Color',lncls(1,:)); hold on; axis tight; box on;
semilogy(HST.time/TimeScale,HST.ke(:,[1,3]),':',LW{:},'Color',lncls(1,:));
set(gca,TL{:},FS{:},'YColor',lncls(1,:),'YScale','log');
yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');

yyaxis right
semilogy(HST.time/TimeScale,HST.etae(:,2),'--' ,LW{:},'Color',lncls(2,:)); hold on; axis tight; box on;
semilogy(HST.time/TimeScale,HST.etae(:,[1,3]),':',LW{:},'Color',lncls(2,:));
set(gca,TL{:},FS{:},'YColor',lncls(2,:),'YScale','log');
yticks = 10.^(-8:2:8);
yticklabels = {'$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$','$10^{8}$'};
set(gca,TL{:},FS{:},'Ytick',yticks,'Yticklabels',yticklabels,'YMinorTick','off');
legend('mean','min/max',TX{:},FS{:},'Location','southeast')
title('Eddy diffusivity [m$^2$/s]  $|$  Eddy viscosity [Pas]',TX{:},FS{1},15);

subplot(2,1,2)
yyaxis left
semilogy(HST.time/TimeScale,HST.ks(:,2),'-' ,LW{:},'Color',lncls(1,:)); hold on; box on;
semilogy(HST.time/TimeScale,HST.ks(:,[1,3]),':',LW{:},'Color',lncls(1,:));
set(gca,TL{:},FS{:},'YColor',lncls(1,:),'YScale','log');

yyaxis right
semilogy(HST.time/TimeScale,HST.etat(:,2),'--' ,LW{:},'Color',lncls(2,:)); hold on; axis tight; box on;
semilogy(HST.time/TimeScale,HST.etat(:,[1,3]),':',LW{:},'Color',lncls(2,:));
set(gca,TL{:},FS{:},'YColor',lncls(2,:),'YScale','log');
legend('mean','min/max',TX{:},FS{:},'Location','southeast')
title(['Segr. diffusivity [m$^2$/s]  $|$  Segr. viscosity [Pas]'],TX{:},FS{1},15);
xlabel(['Time [',TimeUnits,']'],TX{:},FS{1},15);

if plot_cv
if ~exist('fh16','var'); fh16 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh16); clf;
end
plot(HST.time/TimeScale,HST.EB,'k-' ,LW{:}); hold on; axis tight; box on;
plot(HST.time/TimeScale,HST.EM,'k-.',LW{:});
plot(HST.time/TimeScale,HST.EX,'k--',LW{:});
set(gca,TL{:},FS{:})
legend('xtal','melt','mixt',TX{:},FS{:},'Location','northwest')
ylabel('Rel. error [1]',TX{:},FS{1},15);
xlabel(['Time [',TimeUnits,']'],TX{:},FS{1},15);
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
    save(name,'U','W','P','Pt','x','m','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','etas','eII','tII','ke','ks','kx','RaD','ReD','Rs','Red','dt','time','step','MFS','wm','wx','Mx','dMxdt');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','x','m','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','etas','eII','tII','ke','ks','kx','RaD','ReD','Rs','Red','dt','time','step','MFS','wm','wx','Mx','dMxdt');
    name = [outdir,'/',runID,'/',runID,'_HST'];
    save(name,'HST');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
