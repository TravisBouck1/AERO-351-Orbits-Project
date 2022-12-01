%% LEO to Crazy

% RUN JUSTINS CODE FIRST --> GroupProjectScript_JS.m

clc 
% Common Parameters
mu = 398600;        % muEarth
radiusEarth = 6378; % km

%% Convert varables and names from changes upstream
leo1.rp = leo1.peri + radiusEarth;
leo1.ra = leo1.apo + radiusEarth;
leo1.inc = deg2rad(leo1.inc);
leo1.RAAN = deg2rad(leo1.RAAN);
leo1.w = deg2rad(leo1.w);
leo1.Me = deg2rad(leo1.Me);

%% Time for Five LEO Periods

Time.LEOFivePeriods = 5*leo1.T;

%% Crazy Trajectory Parameters

Crazy.rp = 298 + radiusEarth;
Crazy.ra = 33699 + radiusEarth;
Crazy.ecc = (Crazy.ra - Crazy.rp)/(Crazy.ra+Crazy.rp);
Crazy.inc = deg2rad(0.6194);
Crazy.w = deg2rad(134.3946);
Crazy.a = 0.5*(Crazy.ra + Crazy.rp);
Crazy.T = (2*pi*Crazy.a^(3/2))/sqrt(mu);
Crazy.theta = deg2rad(110.976822868941); % Given by Justin
Crazy.h = sqrt(Crazy.rp*mu*(1+Crazy.ecc));
Crazy.RAAN = deg2rad(9.5320);


%% Circularize orbit at perigee of Crazy

% Time calc
Time.Apo2PeriForHohmann1 = leo1.T/2;

% Hohmann Transfer
rp1 = leo1.rp;
ra1 = leo1.ra;
rp3 = Crazy.rp;
ra3 = rp3;

[DeltaV.Hohmann1, Time.Hohmann1,Vp1_1,Vp2_1,deltaV2_1,h1_1,h2_1,h3_1] = Hohmann(leo1.rp,leo1.ra,Crazy.rp,Crazy.rp,mu);

%% Speed Change and inc Change to get into final orbit

% Time to line up
w = Crazy.w - leo1.w;
theta = w + pi; % Need travel from perigee to apogee and then a little more
h_circ = sqrt(2*mu)*sqrt((rp3*ra3)/(rp3 + ra3));
T_circ = (2*pi*rp3^(3/2))/sqrt(mu);
Time.Time2ApseLineAlligned  = ((h_circ^3)/(mu^2))*theta;


% were not really doing a Hohmann but the code solves for stuff we need
[~,~,Vp1,Vp2,deltaV2,h1_2,h2_2,h3_2] = Hohmann(Crazy.rp,Crazy.rp,Crazy.rp,Crazy.ra,mu);

dihedralAngle = leo1.inc - Crazy.inc;
DeltaV.MassiveOne = sqrt(Vp1^2 + Vp2^2 - 2*Vp1*Vp2*cos(dihedralAngle)); % inc and burn to get on crazy orbit

%% Phasing

divisions = 1;

[DeltaV.Phasing,Time.Phasing,hPhasing,TPhasing,eccPhasing,raPhasing] = phasing(Crazy.rp,Crazy.ra,mu,Crazy.theta,divisions);

%% Five Periods

Time.fivePeriods = 5*Crazy.T;

%% Final R and V vector ( just r and v of crazy at perigee)

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, Crazy.h, 0,Crazy.ecc);
vperifocal = v_vect_perifocal(mu, Crazy.h, 0,Crazy.ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(Crazy.RAAN,Crazy.inc,Crazy.w);

% Q transpose * r_peri = r_eci
r_final = (Qperi2geo' * rperifocal)';
v_final =(Qperi2geo' * vperifocal)';

disp("--------------------------")
disp("Final r vect: " + r_final(1) + " "+ r_final(2) + " "+ r_final(3))
disp("Final v vect: "+ v_final(1) + " "+ v_final(2) + " "+ v_final(3))
disp("--------------------------")

%% Final time and deltaV

Time.Total = Time.LEOFivePeriods + Time.Apo2PeriForHohmann1 + Time.Hohmann1 + ...
    Time.Time2ApseLineAlligned + Time.Phasing + Time.fivePeriods;

DeltaV.Total = DeltaV.Hohmann1 + DeltaV.MassiveOne + DeltaV.Phasing;

%% Report total DeltaV and Time since mission start

disp("--------------------------")
disp("deltaV: "+ DeltaV.Total + " km/s")
disp("Time: "+ Time.Total/3600 + " h" )
disp("--------------------------")

%% Hohmann Graph
h = h2_1;
ecc = (leo1.rp-Crazy.rp)/(leo1.rp+Crazy.rp);
theta = 0;
a = 0.5*( Crazy.rp+leo1.rp);
T = (2*pi*a^(3/2))/(sqrt(mu));

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h, theta,ecc);
vperifocal = v_vect_perifocal(mu, h, theta,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(leo1.RAAN,leo1.inc,leo1.w);

% Q transpose * r_peri = r_eci
r_final = (Qperi2geo' * rperifocal);
v_final =(Qperi2geo' * vperifocal);

% Burn
timespan = [0 T/2];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time_burn,Circ_state_burn] = ode45(@coast_ODE, timespan, init, options,mu);

% Circle
h = h3_1;
ecc = 0;
theta = 0;
a = Crazy.rp;
T = (2*pi*a^(3/2))/(sqrt(mu));

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h, theta,ecc);
vperifocal = v_vect_perifocal(mu, h, theta,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(leo1.RAAN,leo1.inc,leo1.w);

% Q transpose * r_peri = r_eci
r_final = (Qperi2geo' * rperifocal);
v_final =(Qperi2geo' * vperifocal);


timespan = [0 T];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time,Circ_state] = ode45(@coast_ODE, timespan, init, options,mu);

figure
h1 = gca;
earth_sphere(h1)
hold on

% LEO Trajectory
plot3(leo1.state_t0(:,1),leo1.state_t0(:,2),leo1.state_t0(:,3),'r','LineWidth',2)

% Circular Trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'b','LineWidth',2)

% Circularise Orbit- First Hohmann Burn
plot3(Circ_state_burn(:,1),Circ_state_burn(:,2),Circ_state_burn(:,3),'k','LineWidth',2)

title("1.1 Hohmann Transfer")


%% Inc and Speed Change Graph

% We are not going to see anything cause we change orbits instantly
figure
h1 = gca;
earth_sphere(h1)
hold on

% Circular Trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'b','LineWidth',2)

% Object II
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'r','LineWidth',2)

title("1.2 Inc and Speed Chnage")


%% Phasing Graph

h = hPhasing;
TA = 0;
ecc = eccPhasing;
RAAN = Crazy.RAAN;
inc = Crazy.inc;
w = Crazy.w;

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h,TA,ecc);
vperifocal = v_vect_perifocal(mu, h,TA,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(RAAN,inc,w);

% Q transpose * r_peri = r_eci
r_epoch = Qperi2geo' * rperifocal;
v_epoch = Qperi2geo' * vperifocal;

% convert to col vects
r_final = r_epoch;
v_final = v_epoch;

timespan = [0 TPhasing];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Phasing_time,Phasing_state] = ode45(@coast_ODE, timespan, init, options,mu);


figure
h1 = gca;
earth_sphere(h1)
hold on

% Crazy Trajectory
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)

% Phasing Trajectory
plot3(Phasing_state(:,1),Phasing_state(:,2),Phasing_state(:,3),'r','LineWidth',2)

title("1.3 Phasing Transfer")







