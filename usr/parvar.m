% crustal magmatism parameter variations  (28 simulations, PA)
D     = [1e-1; 1e+1 (ref); 1e+3;];  % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref);];  % [m]
etam0 = [1e-1; 1e+1 (ref); 1e+3];   % [Pas]
d0    = 3e-2; D = 3e0; etam0 = 1e1; % (settling regime D30_d32_e1)

% additional parameter variations on reference run (11 simulations, TK)
N          = 200;    % [100; 300];
CFL        = 0.5;    % [0.25; 1]
elle       = h/2;    % [1/4; 2]*h
ells       = d0*10;  % [5; 20]*d0
Da         = 1.0;    % [0.1]
Xi         = 1.0;    % [0.1]
xeq        = 0.01;   % [0.1]

% magma ocean parameter variations  (4 simulations, TK)
D     =  1e6;          % [m]
d0    = [1e-2 (ref); 1e-1];  % [m]
etam0 = [1e-1; 1e+1 (ref)];  % [Pas]

% additional parameter variations on magma ocean ref run (5 simulations, TK)
gamma      = 1e-3;   % [0,1e-2];
elle       = h/2;    % [h/4,2*h,4*h];

% total # simulations: 28+11+4+5 = 48

% naming convention for RunIDs
D = 1e1; d = 1e-2; etam0 = 1e+1;  % => D1_d2_e1 reference run
D = 1e6; d = 1e-1; etam0 = 1e-1;  % => D6_d1_em1 magma ocean run
N = 100;                          % => D1_d2_e1_N100 additional tests for ref run
                                  % etc.
