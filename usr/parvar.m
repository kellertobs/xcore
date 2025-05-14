% reference run (similar magnitude of kwx and ke, W and wx, Ra~Re~Ru~1)
D     = 1e1;
d0    = 1e-2;
etam0 = 1e1;

% crustal magmatism parameter variations  (36 simulations)
D     = [1e1 (ref); 1e2; 1e3;]; % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref);];  % [m]
etam0 = [1e0; 1e1 (ref); 1e2; 1e3]; % [Pas]

% magma ocean parameter variations  (4 simulations)
D     = [1e6]; % [m]
d0    = [1e-2; 1e-1];  % [m]
etam0 = [1e-1; 1e0]; % [Pas]

% additional parameter variations (16 simulations)
N          = 200;    % [100; 300];
CFL        = 0.5;    % [0.25; 1]
Delta_cnv  = h;      % [1/3,3]*h
Delta_sgr  = 10*d0;  % [10/3; 30]*d0
Rec        = 1;      % [0.1; 10]
Scx        = 1;      % [10]
Da         = 0.5;    % [1,0.1]
gamma      = 1e-3;   % [0,1e-4,1e-2]
closed_bot = 0;      % [1]

% total # simulations: 36+4+16 = 56

% naming convention for RunIDs
D = 1e1; d = 1e-2; etam0 = 1e1;  % => D1_d2_e1 reference run
D = 1e6; d = 1e-1; etam0 = 1e-1; % => D6_d1_em1 magma ocean run
N = 100;                         % => D1_d2_e1_N100 resolution test of ref run
% etc.



