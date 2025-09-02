% crustal magmatism parameter variations  (18 simulations, PA)
D     = [1e0; 1e1 (ref); 1e2;];  % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref);];  % [m]
etam0 = [1e+1; (ref); 1e+3];   % [Pas]

% additional parameter variations on reference run (13 simulations, TK)
N          = 200;    % [100; 300];
CFL        = 0.5;    % [0.25; 1]
L0         = h/2;    % [1/4; 2]*h
l0         = d0*10;  % [5; 20]*d0
R          = 0.5;    % [1,0.1]
Xi         = 0.5;    % [1,0.1]
xeq        = 0.01;   % [0.1]

% magma ocean parameter variations  (6 simulations, TK)
D     =  1e6;          % [m]
d0    = [1e-3; 1e-2 (ref); 1e-1];  % [m]
etam0 = [1e-1; 1e+1 (ref)];  % [Pas]

% additional parameter variations on magma ocean ref run (7 simulations, TK)
gamma      = 1e-3;   % [0,1e-2];
L0         = h/2;    % [h/4,2*h,4*h];
Xi         = 0.5;    % [[1,0.1]]

% total # simulations: 18+13+6+7 = 44

% naming convention for RunIDs
D = 1e1; d = 1e-2; etam0 = 1e+1;  % => D1_d2_e1 reference run
D = 1e6; d = 1e-1; etam0 = 1e-1;  % => D6_d1_em1 magma ocean run
N = 100;                          % => D1_d2_e1_N100 additional tests for ref run
                                  % etc.
