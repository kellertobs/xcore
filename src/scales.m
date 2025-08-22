%*****  calculate and print characteristic scales  ************************

% length scales
D0      =  D/20;
d0      =  d0;
L0      =  L0;  L0h = (L0+h)/2;
l0      =  l0;  l0h = (l0+h)/2;
h0      =  h;

% material parameter scales
rho0    =  rhom0;
Drho0   =  rhox0-rhom0;
Dchi0   =  xeq/10;
chi0    =  xeq;
eta0    =  etam0;

% speed scales
W0      =  Dchi0*Drho0*g0*D0^2/eta0;  % laminar convection speed
w0      =        Drho0*g0*d0^2/eta0;  % laminar settling speed

ReL0    =  W0*L0/(eta0/rho0);         % laminar convective Reynolds No at L0
Rel0    =  w0*l0/(eta0/rho0);         % laminar settling Reynolds No at L0

fReL    =  1-exp(-ReL0);              % Re-dependent ramp factor
fRel    =  1-exp(-Rel0);              % Re-dependent ramp factor

% general convective speed
if fReL>1e-4; W0  =  D0/(1/2*fReL*L0^2*rho0) * (sqrt(Dchi0*Drho0*fReL*g0*rho0*L0^2*D0 + eta0^2) - eta0); end

% general settling speed
if fRel>1e-4; w0  =  1 /(2  *fRel*l0  *rho0) * (sqrt(4*    Drho0*g0*rho0*fRel*l0*d0^2 + eta0^2) - eta0); end

% diffusivities
eII0    =  W0/D0/4;
ke0     =  eII0*L0^2;
ks0     =  w0*l0;
kx0     =  ks0 + fReL*ke0;

% times
tW0     =  D/W0;
tw0     =  D/w0;
tk0     =  D^2/kx0;
tin0    =  rho0*W0/(Dchi0*Drho0*g0);
t0      =  tin0/2 + min([tW0, tw0, tk0]);
dt0     =  min([(h0/2)^2/kx0 , (h0/2)/(W0+w0)]);

% noise flux amplitudes
xie0    =  Xi*sqrt(     fReL*ke0/(L0h/W0)*(L0h./(L0h+h0))^3);
xiex0   =  Xi*sqrt(chi0*fReL*ke0/(L0h/W0)*(L0h./(L0h+h0))^3);
xis0    =  Xi*sqrt(chi0*     ks0/(l0h/w0)*(l0h./(l0h+h0))^3);
xix0    =  xis0 + xiex0;

% phase change rate
tau0    =  h0./(W0 + w0) + dt0;
G0      =  R*chi0*rho0/tau0;

% viscosities, stress/pressure
etae0   =  fReL*ke0*rho0;
etas0   =  fRel*ks0*rho0;
p0      =  (eta0+etae0)*eII0;

% general dimensionless numbers
Da0     =  G0/(rho0/t0);                % Dahmköhler number
Ne0     =  xie0/W0;                     % Eddy noise number
Ns0     =  xix0/w0;                     % Particle noise number
Rc0     =  W0/w0;                       % Convection number
Ra0     =  W0*D0/kx0;                   % Rayleigh number
ReD0    =  W0*D0/((eta0+etae0)/rho0);   % Convection Reynolds number
Red0    =  w0*d0/((eta0+etas0)/rho0);   % Particle Reynolds number

% print scaling analysis to standard output
fprintf(1,'\n  Scaled domain depth D0    = %1.0e [m]',D0);
fprintf(1,'\n  Crystal size        d0    = %1.0e [m]',d0);
fprintf(1,'\n  Eddy  corrl. length L0    = %1.0e [m]',L0);
fprintf(1,'\n  Segr. corrl. length l0    = %1.0e [m]\n',l0);

fprintf(1,'\n  Density             rho0  = %1.0f  [kg/m3]',rho0);
fprintf(1,'\n  Density contrast    Drho0 = %1.0f   [kg/m3]',Drho0);
fprintf(1,'\n  Cristal. contrast   Dchi0 = %1.3f [wt]',Dchi0);
fprintf(1,'\n  Viscosity           eta0  = %1.0e [Pas]\n',eta0);

fprintf(1,'\n  Convection  speed   W0    = %1.2e [m/s]',W0);
fprintf(1,'\n  Segregation speed   w0    = %1.2e [m/s]\n',w0);

fprintf(1,'\n  Eddy  diffusivity   ke0   = %1.1e [m2/s]',ke0);
fprintf(1,'\n  Segr. diffusivity   ks0   = %1.1e [m2/s]',ks0);
fprintf(1,'\n  Eddy  viscosity     etae  = %1.1e [Pas]',etae0);
fprintf(1,'\n  Segr. viscosity     etas0 = %1.1e [Pas]\n',etas0);

fprintf(1,'\n  Eddy noise rate     xie0  = %1.2e [m/s]',xie0);
fprintf(1,'\n  Segr. noise rate    xix0  = %1.2e [m/s]\n',xix0);

fprintf(1,'\n  Reaction rate       G0    = %1.2e [kg/m3/s]\n',G0);

fprintf(1,'\n  Inertial    time    tin0  = %1.2e [s]',tin0);
fprintf(1,'\n  Convection  time    tW0   = %1.2e [s]',tW0);
fprintf(1,'\n  Segregation time    tw0   = %1.2e [s]',tw0);
fprintf(1,'\n  Diffusion   time    tk0   = %1.2e [s]\n',tk0);

fprintf(1,'\n  Dahmköhler No       Da0   = %1.2e [1]',Da0);
fprintf(1,'\n  Eddy Noise No       Ne0   = %1.2e [1]',Ne0);
fprintf(1,'\n  Settl. Noise No     Ns0   = %1.2e [1]\n',Ns0);

fprintf(1,'\n  Convection No       Rc0   = %1.2e [1]',Rc0);
fprintf(1,'\n  Rayleigh No         Ra0   = %1.2e [1]',Ra0);
fprintf(1,'\n  Domain  Reynolds No ReD0  = %1.2e [1]',ReD0);
fprintf(1,'\n  Crystal Reynolds No Red0  = %1.2e [1]\n\n\n',Red0);


% adjust scales and units for visualisation
if ndm_op
    tsc = t0;  tun = '1';
    Wsc = W0;  Wun = '1';
    wxsc = w0;  wun = '1';
    wmsc = w0*(chi0/(1-chi0));
    whsc = w0; 
    psc = p0;  pun = '1';
    kssc = ks0;  kun = '1';
    kesc = ke0; 
    kxsc = kx0;
    xiesc = xie0;  xieun = '1';
    xixsc = xix0;  xixun = '1';
    esc   = eta0;  eun   = '1';
    eesc  = etae0; 
    essc  = etas0;
    rsc   = rho0;  dun   = '1';
    MFSsc = rho0/t0; MFSun = '1';
    xsc   = chi0;  xun = '1';
    Gsc   = rho0/t0;  Gun = '1';
    ssc   = D0;  sun = '1';
    Rasc  = Ra0;
    ReDsc = ReD0;
    Redsc = Red0;
    Rcsc  = Rc0;
else
    kssc = 1;  kun = 'm$^2$/s';
    kesc = 1;  
    kxsc = 1;  
    esc   = 0;  eun   = 'Pas';
    eesc  = 1;
    essc  = 1;
    rsc   = 1;  dun   = 'kg/m$^3$';
    MFSsc = 1; MFSun = 'kg/m$^3$/s';
    xsc   = 1/100;  xun = 'wt \%';
    Gsc   = 1;  Gun = 'kg/m$^3$/s';
    Rasc  = 1;
    ReDsc = 1;
    Redsc = 1;
    Rcsc  = 1;
if t0 < 1e3
    tsc = 1;
    tun = 's';
elseif t0>= 1e3 && t0 < 1e3*hr
    tsc = hr;
    tun = 'hr';
elseif t0 >= 1e3*hr && t0 < 1e2*yr
    tsc = yr;
    tun = 'yr';
elseif t0 >= 1e2*yr
    tsc = 1e3*yr;
    tun = 'kyr';
end
if D0 < 1e3
    ssc = 1;
    sun = 'm';
elseif D0 >= 1e3
    ssc = 1e3;
    sun = 'km';
end
if W0 < 1000/yr
    Wsc = 1/yr;
    Wun = 'm/yr';
elseif W0 >= 1000/yr && W0 < 1000/hr
    Wsc = 1/hr;
    Wun = 'm/hr';
elseif W0 >= 1000/hr
    Wsc = 1;
    Wun = 'm/s';
end
if w0 < 1000/yr
    wxsc = 1/yr;
    wun  = 'm/yr';
elseif w0 >= 1000/yr && w0 < 1000/hr
    wxsc = 1/hr;
    wun  = 'm/hr';
elseif w0 >= 1000/hr
    wxsc = 1;
    wun  = 'm/s';
end
wmsc = wxsc;
if p0 < 1e2
    psc = 1;
    pun = 'Pa';
elseif p0 >= 1e3 && p0 < 1e6
    psc = 1e3;
    pun = 'kPa';
elseif p0 >= 1e6 && p0 < 1e9
    psc = 1e6;
    pun = 'MPa';
else
    psc = 1e9;
    pun = 'GPa';
end
xiesc = Wsc;  xieun = Wun;
xixsc = Wsc;  xixun = Wun;
whsc  = Wsc;  
end
Xsc = Xc./ssc;
Zsc = Zc./ssc;
Zsf = Zf./ssc;
