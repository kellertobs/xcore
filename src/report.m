% get residual of thermochemical equations
res_mass = norm(upd_X (:))./(norm(rho  (:))+1) ...
         + norm(upd_dV(:))./(norm(Div_V(:))+1);
res_mmnt = norm(upd_W (:))./(norm(V    (:))+1) ...
         + norm(upd_U (:))./(norm(V    (:))+1) ...
         + norm(upd_P (:))./(norm(Pt   (:))+1);

resnorm = res_mass + res_mmnt/1e3;

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
