%*****  calculate and print characteristic scales  ************************

% length scales
h0    =  h;
D0    =  D/10;
d0    =  d0;
L0    =  L0;  L0h = (L0+h0)/2;
l0    =  l0;  l0h = (l0+h0)/2;
bnd_w =  l0/2 + D/100;

% material parameter scales
rho0  =  rhom0;
Drho0 =  rhox0-rhom0;
if Da
     chi0 = Da;
    Dchi0 = Da/10;
else
     chi0 = x0;
    Dchi0 = x0/10;
end
eta0 =  etam0;

W0l  =  Dchi0*Drho0*g0*D0^2/eta0;  % laminar convection speed
w0l  =        Drho0*g0*d0^2/eta0;  % laminar settling speed

w0t  =  sqrt(      Drho0*g0*d0^2/(l0  *rho0));  % terminal turbulent settling speed
W0t  =  sqrt(Dchi0*Drho0*g0*D0^3/(L0^2*rho0));  % terminal turbulent convective speed
W0i  =  sqrt(Dchi0*Drho0*g0*D   /(     rho0));  % inertially limited convective speed
Ri0  =  W0i/W0t;

W0 = W0l;
w0 = w0l;
digits(24);
tol    = 1e-9;
res    = 1;
while res>tol
    W0prv = W0;
    w0prv = w0;

    ReL0  = W0*L0/(eta0/rho0);        % convective Reynolds No at L0, eta0
    Rel0  = w0*l0/(eta0/rho0);        % settling Reynolds No at l0, eta0

    fReL0 = vpa(1-exp(-ReL0));        % Re-dependent ramp factor
    fRel0 = vpa(1-exp(-Rel0));        % Re-dependent ramp factor

    % general convective speed
    W0    = double((sqrt(4./Ri0.^2.*Dchi0.*Drho0.*g0.*rho0.*fReL0.*L0.^2.*D0 + eta0.^2) - eta0).*D0./(2.*fReL0.*L0.^2.*rho0./Ri0.^2));

    % general settling speed
    w0    = double((sqrt(4               .*Drho0.*g0.*rho0.*fRel0.*l0.*d0.^2 + eta0.^2) - eta0)    ./(2.*fRel0.*l0   .*rho0        ));

    res   = abs(W0-W0prv)./W0 + abs(w0-w0prv)./w0;  % residual
end
fRel0 = double(fRel0);
fReL0 = double(fReL0);

% diffusivities
eII0    =  W0/D0;
ke0     =  eII0*L0^2;
ks0     =  w0*l0;
kx0     =  ks0 + fReL0*ke0;

% times
tW0     =  D/W0;
tw0     =  D/w0;
tk0     =  D^2/kx0;
ti0     =  D/W0i;
txi0    =  rhox0*d0^2/18/eta0;
t0      =  (1/(ti0+tW0) + 1/tw0 + 1/tk0).^-1;
dt0     =  min([(h0/2)^2/kx0 , (h0/2)/(W0+w0)]);

% noise flux amplitudes
taue0   =  L0/2/W0;
taus0   =  l0/2/w0;
St0     =  txi0/taue0;
xie0    =  Xi*sqrt(     fReL0*ke0/taue0);
xix0    =  Xi*sqrt(chi0*fReL0*ke0/taue0*St0/(1+St0^2));
xis0    =  Xi*sqrt(chi0*      ks0/taus0);

% phase change rate
G0      =  Da*rho0/t0*D/D0;

% viscosities, stress/pressure
etae0   =  fReL0*ke0*rho0;
etat0   =  fRel0*ks0*rho0;
p0      =  (eta0+etae0)*eII0;

% general dimensionless numbers
Noe0    =  xie0/W0;                     % Mixture-Eddy Noise number
Nox0    =  xix0/w0;                     % Particle-Eddy Noise number
Nos0    =  xis0/w0;                     % Particle-Settling Noise number
Rc0     =  W0/w0;                       % Convection number
Ra0     =  W0*D0/kx0;                   % Rayleigh number
ReD0    =  W0*D0/((eta0+etae0)/rho0);   % Convection Reynolds number
Red0    =  w0*d0/((eta0+etat0)/rho0);   % Particle Reynolds number

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
fprintf(1,'\n  Segr. viscosity     etas0 = %1.1e [Pas]\n',etat0);

fprintf(1,'\n  Mixture-Eddy noise  xie0  = %1.2e [m/s]',xie0);
fprintf(1,'\n  Particle-Eddy noise xix0  = %1.2e [m/s]',xix0);
fprintf(1,'\n  Settling noise      xis0  = %1.2e [m/s]\n',xis0);

fprintf(1,'\n  Reaction rate       G0    = %1.2e [kg/m3/s]\n',G0);

fprintf(1,'\n  Inertial    time    ti0   = %1.2e [s]',ti0);
fprintf(1,'\n  Convection  time    tW0   = %1.2e [s]',tW0);
fprintf(1,'\n  Segregation time    tw0   = %1.2e [s]',tw0);
fprintf(1,'\n  Diffusion   time    tk0   = %1.2e [s]\n',tk0);

fprintf(1,'\n  Dahmköhler No       Da0   = %1.2e [1]',Da);
fprintf(1,'\n  Mixt.-Eddy Noise No Noe0  = %1.2e [1]',Noe0);
fprintf(1,'\n  Part.-Eddy Noise No Nox0  = %1.2e [1]',Nox0);
fprintf(1,'\n  Settling Noise No   Nos0  = %1.2e [1]\n',Nos0);

fprintf(1,'\n  Convection No       Rc0   = %1.2e [1]',Rc0);
fprintf(1,'\n  Rayleigh No         Ra0   = %1.2e [1]',Ra0);
fprintf(1,'\n  Domain  Reynolds No ReD0  = %1.2e [1]',ReD0);
fprintf(1,'\n  Crystal Reynolds No Red0  = %1.2e [1]\n\n\n',Red0);


% adjust scales and units for visualisation
if ndm_op
    Wsc   = W0;  Wpsc = Wsc;  Wun = '1';
    wxsc  = w0;  wxpsc = wxsc; wun = '1'; wpun = '1';
    wmsc  = w0*(chi0/(1-chi0)); wmpsc = wmsc;
    whsc  = w0; 
    psc   = p0;  pun = '1';
    kssc  = ks0;  kun = '1';
    kesc  = ke0; 
    kxsc  = kx0;
    xiesc = xie0;  xieun = '1';
    xissc = xis0;  xisun = '1';
    xixsc = xix0;  xixun = '1';
    esc   = eta0;  eun   = '1';
    eesc  = etae0; 
    etsc  = etat0;
    rsc   = rho0;  dun   = '1';
    MFSsc = rho0/t0; MFSun = '1';
    xsc   = chi0;  xun = '1';
    Gsc   = G0;  Gun = '1';
    ssc   = D;  sun = '1';
    Rasc  = Ra0;
    ReDsc = ReD0;
    Redsc = Red0;
    Rcsc  = Rc0;
    Noesc = Noe0;
    Noxsc = Nox0;
    Nossc = Nos0;
else
    kssc  = 1;  kun = 'm$^2$/s';
    kesc  = 1;  
    kxsc  = 1;  
    esc   = 1;  eun   = 'Pas';
    eesc  = 0;
    etsc  = 0;
    rsc   = 1;  dun   = 'kg/m$^3$';
    MFSsc = 1; MFSun = 'kg/m$^3$/s';
    xsc   = 1/100;  xun = 'wt \%';
    Gsc   = 1;  Gun = 'kg/m$^3$/s';
    Rasc  = 1;
    ReDsc = 1;
    Redsc = 1;
    Rcsc  = 1;
    Noesc = 1;
    Noxsc = 1;
    Nossc = 1;

if D < 1e3
    ssc = 1;
    sun = 'm';
elseif D >= 1e3 && D < 1e6
    ssc = 1e3;
    sun = 'km';
else
    ssc = 1e6;
    sun = 'Mm';
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

if max([W0,w0]) < 1000/yr
    Wpsc = 1/yr; wxpsc = Wpsc; wmpsc = Wpsc;
    wpun = 'm/yr';
elseif max([W0,w0]) >= 1000/yr && max([W0,w0]) < 1000/hr
    Wpsc = 1/hr; wxpsc = Wpsc; wmpsc = Wpsc;
    wpun = 'm/hr';
elseif max([W0,w0]) >= 1000/hr
    Wpsc = 1; wxpsc = Wpsc; wmpsc = Wpsc;
    wpun = 'm/s';
end

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
xissc = wxsc;  xisun = wun;
xixsc = Wsc;  xixun = Wun;
whsc  = Wsc;  
end
Xsc = Xc./ssc;
Zsc = Zc./ssc;
Zsf = Zf./ssc;
