%*****  initialise model run  *********************************************
init;

%*****  physical time stepping loop  **************************************
while time <= tend && step <= Nt
    
    %***  time step info
    timing;

    %***  store previous solution
    store;
    
    %***  reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;
    if frst; alpha = alpha/2; beta = beta/2; end

    %***  non-linear iteration loop
    while resnorm/resnorm0 >= rtol/(1 + frst*10) && resnorm >= atol/(1 + frst*10) && iter <= maxit*(1 + frst)
        
        %***  solve phase evolution equations
        phsevo;

        %***  solve fluid-mechanics equations
        fluidmech;

        %***  update non-linear parameters and auxiliary variables
        update;

        %***  report convergence
        report;

        iter = iter+1;  % increment iteration count

    end % end non-linear iterations

    %***  update correlation length for convective/turbulent regularisation
    corrl;

    %***  record model history
    if ~mod(step,nrh); history; end

    %***  print model diagnostics
    diagnose;

    %***  plot model results
    if ~mod(step,nop); output; end

    %***  increment time/step
    time = time+dt;
    step = step+1;
    if frst; alpha = alpha*2; beta = beta*2; frst=0; end
    
end % end time stepping

%***  save final state of model
output;

diary off
