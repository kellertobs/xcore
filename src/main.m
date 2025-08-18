%*****  XCORE MAIN MODEL ROUTINE  *****************************************

%*****  initialise model run  *********************************************

init;

%*****  physical time stepping loop  **************************************
while time <= tend && time/t0 <= t0end && step <= Nt && HST.x(end,2) <= xend
   
    
    %***  time step info
    timing;

    %***  store previous solution
    store;
    
    %***  reset residuals and iteration count
    resnorm  = 1;
    resnorm0 = resnorm;
    iter     = 1;

    %***  non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit
        
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

    %***  record model history
    if ~mod(step,nrh); history; end

    %***  print model diagnostics
    diagnose;

    %***  plot model results
    if ~mod(step,nop); output; end

    %***  increment time/step
    time = time+dt;
    step = step+1;
    
end % end time stepping

%***  save final state of model
output;

diary off
