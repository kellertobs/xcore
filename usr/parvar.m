% reference run (similar magnitude of kwx and ke, W and wx, Ra~Re~Ru~1)
D     = 1e+1;
d0    = 1e-2;
etam0 = 1e+1;

% crustal magmatism parameter variations  (27 simulations, PA)
D     = [1e+1 (ref); 1e+2; 1e+3]; % [m]
d0    = [1e-4; 1e-3; 1e-2 (ref)]; % [m]
etam0 = [1e+1 (ref); 1e+2; 1e+3]; % [Pas]

% settling-dominated regime demo (4 simulations, PA)
D     = [5];          % [m]
d     = [1e-2; 1(3?)e-3]; % [m]
etam0 = [1e+0; 1e+1]  % [Pas]
% Additional settling: (2 simulations, PA)
D = [5]; % [m]
d = [1e-4; 3e-2] % [m]
etam0 = [1e+0]; % [Pas]

% magma ocean parameter variations  (6 simulations, TK)
D     =  1e6;          % [m]
d0    = [1e-2 (ref); 1e-1];  % [m]
etam0 = [1e-1 (ref); 1e+0];  % [Pas]
Rex   = 1;             % [0.1; 10]

% additional parameter variations (12 simulations, TK)
N          = 200;    % [100; 300];
CFL        = 0.5;    % [0.25; 1]
% Delta_cnv  = h;      % [1/3; 3]*h  % needed with the new implementation ?
% Delta_sgr  = 20*d0;  % [10; 40]*d0 % needed with the new implementation ?
Rec        = 1;      % [0.1; 10]
Scx        = 1;      % [10]
Da         = 0.5;    % [1,0.25]
gamma      = 1e-3;   % [1e-4,1e-2]; D = 1e3;
% closed_bot = 0;      % [1] % dont think we need this in the end!
xend       = xeq/2 = 0.05;

% total # simulations: 27+6+6+12 = 51

% naming convention for RunIDs
D = 1e1; d = 1e-2; etam0 = 1e+1;  % => D1_d2_e1 reference run
D = 1e6; d = 1e-1; etam0 = 1e-1;  % => D6_d1_em1 magma ocean run
N = 100;                          % => D1_d2_e1_N100 resolution test of ref run
% etc.
