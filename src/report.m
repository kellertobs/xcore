% get residual of thermochemical equations
res_mass = norm(upd_X  (:))./(norm(rho(:))+eps) ...
         + norm(upd_MFS(:))./(norm(MFS(:))+eps);
res_mmnt = norm(upd_W  (:))./(norm(W  (:))+eps) ...
         + norm(upd_U  (:))./(norm(U  (:))+eps) ...
         + norm(upd_P  (:))./(norm(P  (:))+eps);

resnorm = res_mass + res_mmnt;

if iter==1 || resnorm>resnorm0; resnorm0 = resnorm + 1e-32; end  % reset reference residual

% check for solver divergence or failing
if isnan(resnorm); error('!!! Solver failed with NaN: end run !!!'); end

% report iterations
if     iter >=  0  && iter <  10
    fprintf(1,'    ---  iter =   %d;  abs = %1.2e;  rel = %1.2e;  mass = %1.2e;  mmnt = %1.2e \n',iter,resnorm,resnorm/resnorm0,res_mass,res_mmnt);
elseif iter >= 10  && iter < 100
    fprintf(1,'    ---  iter =  %d;  abs = %1.2e;  rel = %1.2e;  mass = %1.2e;  mmnt = %1.2e \n',iter,resnorm,resnorm/resnorm0,res_mass,res_mmnt);
elseif iter >= 100 && iter < 1000
    fprintf(1,'    ---  iter = %d;  abs = %1.2e;  rel = %1.2e;  mass = %1.2e;  mmnt = %1.2e \n',iter,resnorm,resnorm/resnorm0,res_mass,res_mmnt);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if iter==1; clf; else; hold on; end
    plot(iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
    drawnow;
end
