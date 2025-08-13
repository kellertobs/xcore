% print diagnostics
fprintf(1,'\n         total time to solution = %3.3f sec\n\n',toc(TTtime));
fprintf(1,'         fluid-mechanics solve  = %1.3e sec\n'  ,FMtime/(iter-1));
fprintf(1,'         phase evolution solve  = %1.3e sec\n'  ,TCtime/(iter-1));
fprintf(1,'         coefficients update    = %1.3e sec\n\n',UDtime/(iter-1));

fprintf(1,'         min x   = %1.6f;  mean x   = %1.6f;  max x   = %1.6f;   [%s]\n'   ,min(x(:))/xsc,mean(x(:))/xsc,max(x(:))/xsc,xun(xun~='\'&xun~=' '));
fprintf(1,'         min m   = %1.5f;  mean m   = %1.5f;  max m   = %1.5f;   [%s]\n\n' ,min(m(:))/xsc,mean(m(:))/xsc,max(m(:))/xsc,xun(xun~='\'&xun~=' '));

fprintf(1,'         min ks  = %1.2e;  mean ks  = %1.2e;  max ks  = %1.2e;   [%s]\n',min(ks(:))/kssc,geomean(ks(:))/kssc,max(ks(:))/kssc,kun(kun~='$'&kun~='^'));
fprintf(1,'         min ke  = %1.2e;  mean ke  = %1.2e;  max ke  = %1.2e;   [%s]\n\n',min(ke(:))/kesc,geomean(ke(:))/kesc,max(ke(:))/kesc,kun(kun~='$'&kun~='^'));
if ndm_op
fprintf(1,'         min rho = %1.6f;  mean rho = %1.6f;  max rho = %1.6f;   [%s]\n'  ,min(rho(:))/rsc,geomean(rho(:))/rsc,max(rho(:))/rsc,run(run~='$'&run~='^'));
else
fprintf(1,'         min rho = %4.3f;  mean rho = %4.3f;  max rho = %4.3f;   [%s]\n'  ,min(rho(:))/rsc,geomean(rho(:))/rsc,max(rho(:))/rsc,run(run~='$'&run~='^'));
end 
fprintf(1,'         min eta = %1.2e;  mean eta = %1.2e;  max eta = %1.2e;   [%s]\n\n',min(eta(:))/esc,geomean(eta(:))/esc,max(eta(:))/esc,eun);

fprintf(1,'         min V   = %1.2e;  mean V   = %1.2e;  max V   = %1.2e;   [%s]\n'  ,min(V(:))/Wsc,geomean(V(:))/Wsc,max(V(:))/Wsc,Wun);
fprintf(1,'         min vx  = %1.2e;  mean vx  = %1.2e;  max vx  = %1.2e;   [%s]\n'  ,min(vx(:))/wxsc,geomean(vx(:))/wxsc,max(vx(:))/wxsc,wun);
fprintf(1,'         min vm  = %1.2e;  mean vm  = %1.2e;  max vm  = %1.2e;   [%s]\n\n'  ,min(vm(:))/wmsc,geomean(vm(:))/wmsc,max(vm(:))/wmsc,wun);

fprintf(1,'         min xix = %1.2e;  mean xix = %1.2e;  max xix = %1.2e;   [%s]\n'  ,min(xix(:))/xixsc,geomean(xix(:))/xixsc,max(xix(:))/xixsc,Wun);
fprintf(1,'         min xie = %1.2e;  mean vm  = %1.2e;  max xie = %1.2e;   [%s]\n\n'  ,min(xie(:))/xiesc,geomean(xie(:))/xiesc,max(xie(:))/xiesc,Wun);

fprintf(1,'         min Rc  = %1.2e;  mean Rc  = %1.2e;  max Rc  = %1.2e;   [1]\n'  ,min(Rc(:))/Rcsc,geomean(Rc(:))/Rcsc,max(Rc(:))/Rcsc);
fprintf(1,'         min Ra  = %1.2e;  mean Ra  = %1.2e;  max Ra  = %1.2e;   [1]\n'  ,min(Ra(:))/Rasc,geomean(Ra(:))/Rasc,max(Ra(:))/Rasc);
fprintf(1,'         min ReD = %1.2e;  mean ReD = %1.2e;  max ReD = %1.2e;   [1]\n'  ,min(ReD(:))/ReDsc,geomean(ReD(:))/ReDsc,max(ReD(:))/ReDsc);
fprintf(1,'         min Red = %1.2e;  mean Red = %1.2e;  max Red = %1.2e;   [1]\n\n\n'  ,min(Red(:))/Redsc,geomean(Red(:))/Redsc,max(Red(:))/Redsc);
