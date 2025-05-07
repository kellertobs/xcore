% main parameter variations  (64 simulations)
D     = [1e0; 1e2; 1e4; 1e6]; % [m]
d0    = [1e-4; 1e-3; 1e-2; 1e-1];  % [m]
etam0 = [1e-1; 1e1; 1e3; 1e5]; % [Pas]

% reference run
D     = 1e2;
d0    = 1e-2;
etam0 = 1e1;

% additional parameter variations (12 simulations)
N         = 200;    % [100; 300]; [200x400];
CFL       = 0.5;    % [0.25; 1]
Delta_cnv = h;      % [1/3,3]*h
Delta_sgr = 10*d0;  % [10/3; 30]*d0
Rec       = 1;      % [0.1; 10]
Scx       = 1;      % [10]

% potential additional variations (4 simulations)
x0        = 0.01;   % [1e-5; 0.1]
xeq       = 0.1;    % [      0.2]
Da        = 0.1;    % [1; 0.01]


