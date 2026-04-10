% crustal magmatism parameter variations  (27 simulations, PA)
D     = [1e0; 1e1 (ref); 1e2;];  % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref);];  % [m]
etam0 = [1e+1 (ref); 1e+2; 1e+3];   % [Pas]

% additional parameter variations on reference run (14 simulations, TK)
N     = 200;    % [100; 300];
CFL   = 0.5;    % [0.25; 1]
L0    = D/100;  % [1/200; 1/50]*D
l0    = d0*10;  % [5; 20]*d0
Da    = 0.01;   % [0.1]
Xi    = 0.5;    % [1,0.1,0]
rtol  = 1e-4;   % [1e-3,1e-5]

% additional run inertial settling regime (1 simulation, TK)
D     = 10;
d0    = 0.1;
etam0 = 0.1;

% additional runs for cumulate stack to t0end = 200 (5 simulations, TK)
D = 1e+1; d0 = 1e-2; etam0 = 1e+1;  % reference, regime transition
D = 1e+0; d0 = 1e-2; etam0 = 1e+1;  % laminar settling regime
D = 1e+1; d0 = 1e-1; etam0 = 1e-1;  % inertial settling regime
D = 1e+1; d0 = 1e-4; etam0 = 1e+3;  % laminar convection regime
D = 1e+6; d0 = 1e-1; etam0 = 1e-1;  % inertial convection regime

% magma ocean parameter variations  (9 simulations, TK)
D     =  1e6;          % [m]
d0    = [1e-3; 1e-2; 1e-1 (ref)];  % [m]
etam0 = [1e-1 (ref); 1e+0; 1e+1];  % [Pas]

% additional parameter variations on magma ocean ref run (7 simulations, TK)
gamma = 1e-3;   % [0,1e-2];
L0    = D/100;  % [D/200,D/50,D/10];
Xi    = 0.5;    % [1,0.1]

% total # simulations: 27+14+1+5+9+7 = 63

% naming convention for RunIDs
D = 1e1; d0 = 1e-2; etam0 = 1e+1;  % => D1_dm2_e1 reference run
D = 1e6; d0 = 1e-1; etam0 = 1e-1;  % => D6_dm1_em1 magma ocean run
N = 100;                           % => D1_dm2_e1_N100 additional tests run
