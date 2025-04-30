
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
if time < 1e3*hr
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
if max(Vel(:)) < 1000/yr
    SpeedScale = 1/yr;
    SpeedUnits = 'm/yr';
elseif max(Vel(:)) >= 1000/yr && max(Vel(:)) < 1000/hr
    SpeedScale = 1/hr;
    SpeedUnits = 'm/hr';
elseif max(Vel(:)) >= 1000/hr
    SpeedScale = 1;
    SpeedUnits = 'm/s';
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

if Nx <= 1 && Nz <= 1  % create 0D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    subplot(3,1,1)
    plot(hist.time/tau_T,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/tau_T,hist.Tliq(:,2),CL{[1,3]},LW{:});
    plot(hist.time/tau_T,hist.Tsol(:,2),CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});

    subplot(3,1,2)
    plot(hist.time/tau_T,hist.mu (:,2)*100.*(hist.mu (:,2)>eps^0.5),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/tau_T,hist.chi(:,2)*100.*(hist.chi(:,2)>eps^0.5),CL{[1,4]},LW{:});
    plot(hist.time/tau_T,hist.phi(:,2)*100.*(hist.phi(:,2)>eps^0.5),CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});

    subplot(3,1,3)
    plot(hist.time/tau_T,    hist.Gm(:,2)./hist.rho(:,2)*hr*100.*(hist.mu (:,2)>eps^0.5),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/tau_T,    hist.Gx(:,2)./hist.rho(:,2)*hr*100.*(hist.chi(:,2)>eps^0.5),CL{[1,4]},LW{:});
    plot(hist.time/tau_T,10.*hist.Gf(:,2)./hist.rho(:,2)*hr*100.*(hist.phi(:,2)>eps^0.5),CL{[1,5]},LW{:});
    title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt \%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel(['Time [$\tau_T$]'],TX{:},FS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end
    subplot(3,1,1)
    plot(hist.time/tau_T,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on
    plot(hist.time/tau_T,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
    plot(hist.time/tau_T,hist.rho (:,2),'-',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(3,1,2)
    semilogy(hist.time/tau_T,hist.eta (:,2),'k-',LW{:}); axis xy tight; box on; hold on
    semilogy(hist.time/tau_T,hist.etam(:,2),'r-',LW{:});
    title('$\eta^m,\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(3,1,3)
    plot(hist.time/tau_T,hist.dV(:,2)*100*TimeScale,'k-',LW{:}); axis xy tight; box on;
    title(['$\dot{V}$ [\%/',TimeUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel(['Time [$\tau_T$]'],TX{:},FS{:});

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end
    subplot(3,1,1)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.c(:,2,:)).'); axis xy tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    Ttick = get(gca,'Xtick');
    Ptick = (hist.Pt(1) + (Ttick-hist.T(1)+273.15)*dPdT)/1e8;
    for i = 1:length(Ttick)
        TPtick{i} = [num2str(Ttick(i)),'; ',num2str(Ptick(i),3)]; 
    end
    title('Bulk CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,2)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cm(:,2,:)).'); axis xy tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    title('Melt CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,3)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cx(:,2,:)).'); axis xy tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    title('Xtal CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    xlabel(['T [$^\circ$C]; P [kbar]'],TX{:},FS{:});

    if ~exist('fh4','var'); fh4 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh4); clf;
    end
    subplot(3,1,1)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.c_oxd(:,2,:)/100).'); axis xy tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Bulk OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,2)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cm_oxd(:,2,:)/100).'); axis xy tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Melt OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,3)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cx_oxd(:,2,:)/100).'); axis xy tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Xtal OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    xlabel(['T [$^\circ$C]; P [kbar]'],TX{:},FS{:});

    if ~exist('fh5','var'); fh5 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh5); clf;
    end
    subplot(3,1,1)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.c_mem(:,2,:)/100).'); axis xy tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Bulk MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,2)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cm_mem(:,2,:)/100).'); axis xy tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Melt MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    subplot(3,1,3)
    imagesc(hist.T(:,2)-273.15,1:100,comp_plt(hist.cx_mem(:,2,:)/100).'); axis xy tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Xtal MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:},XR{:},'XTickLabel',TPtick);
    xlabel(['T [$^\circ$C]; P [kbar]'],TX{:},FS{:});

    if (fractxtl || fractmlt) && step>1
        if ~exist('fh6','var'); fh6 = figure(VIS{:});
        else; set(0, 'CurrentFigure', fh6); clf;
        end
        subplot(3,1,1)
        imagesc(cumsum(CML.d),1:100,comp_plt(CML.c        ).'); axis xy tight; colormap(gca,colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
        title('Cumulate CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(3,1,2)
        imagesc(cumsum(CML.d),1:100,comp_plt(CML.c_mem/100).'); axis xy tight; colormap(gca,colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
        title('Cumulate MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(3,1,3)
        imagesc(cumsum(CML.d),1:100,comp_plt(CML.c_oxd/100).'); axis xy tight; colormap(gca,colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
        title('Cumulate OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel(['Thickness [1]'],TX{:},FS{:});

        if ~exist('fh7','var'); fh7 = figure(VIS{:});
        else; set(0, 'CurrentFigure', fh7); clf;
        end
        subplot(3,1,1)
        colororder({'k','k'})
        yyaxis left
        plot(cumsum(CML.d),hist.T(:,2)-273.15,'k-',LW{:}); axis xy tight; box on;
        yyaxis right
        plot(cumsum(CML.d),hist.Pt(:,2)/1e8,'k-',LW{:}); axis xy tight; box on;
        title(['Cumulate T [$^\circ$C]  $|$  P [kbar]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(3,1,2)
        plot(cumsum(CML.d),CML.rho,'-',CL{[1,2]},LW{:}); axis xy tight; box on;
        title('Cumulate $\rho$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(3,1,3)
        semilogy(cumsum(CML.d),CML.eta,'-',CL{[1,2]},LW{:}); axis xy tight; box on;
        title('Cumulate $\eta$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel(['Thickness [1]'],TX{:},FS{:});

        if dPdT
            if ~exist('fh8','var'); fh8 = figure(VIS{:});
            else; set(0, 'CurrentFigure', fh8); clf;
            end
            for i=1:cal.ncmp-1
                plot(hist.Tm(:,i),hist.Pt(:,2)/1e8,'k','Color',[1 1 1].*i./cal.ncmp,LW{1},1.5); axis ij tight; box on; hold on
            end
            plot(hist.T(:,2)-273.15,hist.Pt(:,2)/1e8,CL{[1,2]},LW{:});
            plot(hist.Tliq(:,2),hist.Pt(:,2)/1e8,CL{[1,3]},LW{:});
            plot(hist.Tsol(:,2),hist.Pt(:,2)/1e8,CL{[1,4]},LW{:});
            set(gca,TL{:},TS{:});
            legend([cal.cmpStr(1:end-1),'$P,T$-path','$T_\mathrm{liq}$','$T_\mathrm{sol}$'],TX{:},FS{:})
            ylabel(['Pressure [kbar]'],TX{:},'FontSize',15);
            xlabel(['Temperature [$^\circ$C]'],TX{:},'FontSize',15);
        end
    end

elseif Nx <= 1  % create 1D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    plot(T-273.15,Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
    plot(cal.Tliq,Zsc.',CL{[1,3]},LW{:});
    plot(cal.Tsol,Zsc.',CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    plot(mu *100.*(mu >eps^0.5),Zsc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(chi*100.*(chi>eps^0.5),Zsc.',CL{[1,4]},LW{:});
    plot(phi*100.*(phi>eps^0.5),Zsc.',CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    plot(    Gm./rho*100*hr.*(mu >eps^0.5),Zsc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(    Gx./rho*100*hr.*(chi>eps^0.5),Zsc.',CL{[1,4]},LW{:});
    plot(10.*Gf./rho*100*hr.*(phi>eps^0.5),Zsc.',CL{[1,5]},LW{:});
    title(['$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    plot(-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1)/SpeedScale,Zsf.',CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
    plot(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)/SpeedScale,Zsf.',CL{[1,4]},LW{:});
    plot(-(mu ([1,1:end],:)+mu ([1:end,end],:))/2.*wm(:,2:end-1)/SpeedScale,Zsf.',CL{[1,3]},LW{:});
    plot(-                                         W (:,2:end-1)/SpeedScale,Zsf.',CL{[1,2]},LW{:});
    title(['$W$, $w_\Delta^f$, $w_\Delta^x$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end 
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    plot(rhox,Zsc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(rhom,Zsc.',CL{[1,3]},LW{:});
    plot(rho ,Zsc.',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    semilogx(etax,Zsc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    semilogx(etam,Zsc.',CL{[1,3]},LW{:});
    semilogx(eta ,Zsc.',CL{[1,2]},LW{:});
    title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    plot(dV*hr,Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title(['$\dot{V}$ [1/hr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    plot(P(2:end-1,2:end-1),Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
 
    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end 
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,3,1)
    imagesc(1:100,Zsc.',comp_plt(c)); axis ij tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    title('Bulk CMPs [wt\%]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,2)
    imagesc(1:100,Zsc.',comp_plt(cm)); axis ij tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    title('Melt CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,3)
    imagesc(1:100,Zsc.',comp_plt(cx)); axis ij tight; colormap(colmapcmp); clim([1,cal.ncmp]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.ncmp,'Ticklabels',cal.cmpStr);
    title('Xtal CMPs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    
    if ~exist('fh4','var'); fh4 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh4); clf;
    end
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,3,1)
    imagesc(1:100,Zsc.',comp_plt(c_oxd/100)); axis ij tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Bulk OXDs [wt\%]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,2)
    imagesc(1:100,Zsc.',comp_plt(cm_oxd/100)); axis ij tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Melt OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,3)
    imagesc(1:100,Zsc.',comp_plt(cx_oxd/100)); axis ij tight; colormap(colmapoxd); clim([1,cal.noxd]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.noxd,'Ticklabels',cal.oxdStr);
    title('Xtal OXDs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});

    if ~exist('fh5','var'); fh5 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh5); clf;
    end
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,3,1)
    imagesc(1:100,Zsc.',comp_plt(c_mem/100)); axis ij tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Bulk MEMs [wt\%]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,2)
    imagesc(1:100,Zsc.',comp_plt(cm_mem/100)); axis ij tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Melt MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,3,3)
    imagesc(1:100,Zsc.',comp_plt(cx_mem/100)); axis ij tight; colormap(colmapmem); clim([1,cal.nmem]); cb = colorbar; set(cb,TL{:},TS{:},'Ticks',1:cal.nmem,'Ticklabels',cal.memStr);
    title('Xtal MEMs [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});

else % create 2D plots

    % initialize figures and axes
    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end 
    colormap(colmap);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh1,UN{:},'Position',[1 1 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off','Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

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
    imagesc(Xsc,Zsc,-W(:      ,2:end-1)./SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh1,'CurrentAxes',ax(12));
    imagesc(Xsc,Zsc, U(2:end-1,:      )./SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh1,'CurrentAxes',ax(13));
    imagesc(Xsc,Zsc, P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [Pa]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh1,'CurrentAxes',ax(14));
    imagesc(Xsc,Zsc,Div_V*TimeScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \mathbf{v}$ [1/',TimeUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot phase fractions and reaction rates in Fig. 3
    set(0,'CurrentFigure',fh3)
    set(fh3,'CurrentAxes',ax(31));
    imagesc(Xsc,Zsc,chi.*100.*(chi>eps^0.5) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh3,'CurrentAxes',ax(32));
    imagesc(Xsc,Zsc,mu .*100.*(mu >eps^0.5)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\mu$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh3,'CurrentAxes',ax(33));
    imagesc(Xsc,Zsc,Gx./rho*hr*100.*(chi>eps^0.5)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(34));
    imagesc(Xsc,Zsc,Gm./rho*hr*100.*(mu >eps^0.5)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_m/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot density, rheology, and segregation speeds in Fig. 4
    set(0,'CurrentFigure',fh4)
    set(fh4,'CurrentAxes',ax(41));
    imagesc(Xsc,Zsc,rho-mean(rho,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta \bar{\rho}_h$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh4,'CurrentAxes',ax(42));
    imagesc(Xsc,Zsc,log10(eta)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh4,'CurrentAxes',ax(43));
    imagesc(Xsc,Zsc,-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1).*1e3/SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m',SpeedUnits,']'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(44));
    imagesc(Xsc,Zsc,-(mu ([1,1:end],:)+mu ([1:end,end],:))/2.*wm(:,2:end-1).*1e3/SpeedScale.*(any(mu>eps^0.5))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^m$ [m',SpeedUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

end

% plot model history
if plot_cv
    if ~exist('fh14','var'); fh13 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh13); clf;
    end 
    subplot(4,1,1);
    plot(hist.time/TimeScale,hist.DS./hist.sumS,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    plot(hist.time/TimeScale,hist.DB./hist.sumB,'k-',LW{:}); hold on; axis tight; box on; hold on
    plot(hist.time/TimeScale,hist.DM./hist.sumB,'k--',LW{:}); hold on; axis tight; box on;
    plot(hist.time/TimeScale,hist.DX./hist.sumB,'k-.',LW{:}); hold on; axis tight; box on;
    plot(hist.time/TimeScale,hist.DF./hist.sumB,'k:',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $\bar{\rho},F^i$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    plot(hist.time/TimeScale,hist.DC./hist.sumC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $C_j$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    plot(hist.time/TimeScale,hist.DT./hist.sumT,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $\Theta_k$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});

    if ~exist('fh15','var'); fh14 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh14); clf;
    end 
    subplot(4,1,1);
    plot(hist.time/TimeScale,hist.ES,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    plot(hist.time/TimeScale,hist.EB,'k-',LW{:}); hold on; axis tight; box on; hold on
    plot(hist.time/TimeScale,hist.EM,'k--',LW{:}); hold on; axis tight; box on;
    plot(hist.time/TimeScale,hist.EX,'k-.',LW{:}); hold on; axis tight; box on;
    plot(hist.time/TimeScale,hist.EF,'k:',LW{:}); hold on; axis tight; box on;
    ylabel('error $\bar{\rho},F^i$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    plot(hist.time/TimeScale,hist.EC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $C_j$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    plot(hist.time/TimeScale,hist.ET,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $\Theta_k$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});
end

drawnow

% save output to file
if save_op && ~restart
    if Nx <= 1  % print 0D/1D figures
        name = [outdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_oxd_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_mem_',num2str(floor(step/nop))];
        print(fh5,name,'-dpng','-r300','-image');
        if (fractxtl || fractmlt) && step>1 && Nz==1
            name = [outdir,'/',runID,'/',runID,'_cml_',num2str(floor(step/nop))];
            print(fh6,name,'-dpng','-r300','-image');
            name = [outdir,'/',runID,'/',runID,'_clp_',num2str(floor(step/nop))];
            print(fh7,name,'-dpng','-r300','-image');
            if dPdT
                name = [outdir,'/',runID,'/',runID,'_PT_',num2str(floor(step/nop))];
                print(fh8,name,'-dpng','-r300','-image');
            end
        end
    else      % print 2D figures
        name = [outdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
    end
    
    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','x','m','xq','mq','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','eII','tII','dt','time','step','dV','wx','wm');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','x','m','xq','mq','chi','mu','X','M','dXdt','dMdt','drhodt','Gx','Gm','rho','eta','eII','tII','dt','time','step','dV','wx','wm');
    name = [outdir,'/',runID,'/',runID,'_hist'];
    save(name,'hist');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
    