% crustal magmatism parameter variations  (27 simulations, PA)
D     = [1e+0; 1e+1 (ref); 1e+2;]; % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref)]; % [m]
etam0 = [1e-1; 1e+1 (ref); 1e+3]; % [Pas]

% additional parameter variations on reference run (10 simulations, TK)
N          = 200;    % [100; 300];
CFL        = 0.5;    % [0.25; 1]
elle       = h/2;    % [1/4; 2]*h
ells       = 10*d0;  % [5; 20]*d0
Da         = 1.0;    % [0.1]
Xi         = 
xeq        = 0.01;   % [0.1]

% magma ocean parameter variations  (4 simulations, TK)
D     =  1e6;          % [m]
d0    = [1e-2 (ref); 1e-1];  % [m]
etam0 = [1e-1 (ref); 1e+1];  % [Pas]

% additional parameter variations on magma ocean ref run (2 simulations
gamma      = 1e-3;   % [0,1e-2];
elle       = h/2;    % [0,1e-2];

% total # simulations: 27++6+12 = 51

% naming convention for RunIDs
D = 1e1; d = 1e-2; etam0 = 1e+1;  % => D1_d2_e1 reference run
D = 1e6; d = 1e-1; etam0 = 1e-1;  % => D6_d1_em1 magma ocean run
N = 100;                          % => D1_d2_e1_N100 resolution test of ref run
                                  % etc.
