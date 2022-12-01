%% Testing it out 2


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


%% Circularize orbit 

% Time calc
Time.Apo2PeriForHohmann1 = leo1.T/2;

% Hohmann Transfer
rp1 = leo1.rp;
ra1 = leo1.ra;
rp3 = Crazy.rp;
ra3 = rp3;

[DeltaV.Hohmann1, Time.Hohmann1,Vp1_1,Vp2_1,deltaV2_1,h1_1,h2_1,h3_1] = Hohmann(ra1,rp1,rp3,ra3,mu);









h =  sqrt(2*mu)*sqrt((leo1.rp*Crazy.rp)/(leo1.rp + Crazy.rp));
ecc = (leo1.rp-Crazy.rp)/(leo1.rp+Crazy.rp);
theta = pi;
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
timespan = [T T/2];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time_burn,Circ_state_burn] = ode45(@coast_ODE, timespan, init, options,mu);

Circ_state_burn(:,1) = Circ_state_burn(:,1) + (-3080.17 + 3287.72);
Circ_state_burn(:,2) = Circ_state_burn(:,2) + (5549.06 -5956.57);
Circ_state_burn(:,3) = Circ_state_burn(:,3) + (2071.1 - 2217.71);


% Circle
h = h3_1;
ecc = 0;
theta = pi;
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
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Circularise Orbit- First Hohmann Burn
plot3(Circ_state_burn(:,1),Circ_state_burn(:,2),Circ_state_burn(:,3),'k','LineWidth',2)

% Hohmann Burn 1
plot3(Circ_state_burn(1,1),Circ_state_burn(1,2),Circ_state_burn(1,3),'*','LineWidth',5)

% Hohmann Burn 2
plot3(Circ_state_burn(end,1),Circ_state_burn(end,2),Circ_state_burn(end,3),'*','LineWidth',5)


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Hohmann Transfer');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','LEO1 orbit','Circular Transfer Trajectory','Hohmann Transfer Trajectory','Burn #2','Burn #1','best')




figure
hold on

% LEO Trajectory
plot(leo1.state_t0(:,1),leo1.state_t0(:,2),'r','LineWidth',2)

% Circular Trajectory
plot(Circ_state(:,1),Circ_state(:,2),'c','LineWidth',2)

% Circularise Orbit- First Hohmann Burn
plot(Circ_state_burn(:,1),Circ_state_burn(:,2),'k','LineWidth',2)

% Hohmann Burn 1
plot(Circ_state_burn(1,1),Circ_state_burn(1,2),'*','LineWidth',5)

% Hohmann Burn 2
plot(Circ_state_burn(end,1),Circ_state_burn(end,2),'*','LineWidth',5)


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Hohmann Transfer X-Y Plane');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('LEO1 orbit','Circular Transfer Trajectory','Hohmann Transfer Trajectory','Burn #2','Burn #1','best')











%% Find 3 points on circle orbit and crazy orbit

% Circle
h = h3_1;
ecc = 0;
theta = pi;
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



% three points for circle
r11 = Circ_state(1,1:3);
r12 = Circ_state(20,1:3);
r13 = Circ_state(140,1:3);

P = r12 - r11;
Q = r13 - r11;

n1 = cross(P,Q);

const = (n1(1)*r11(1) + n1(2)*r11(2) + n1(3)*r11(3)); 

% plane equation
P1 = [n1,const];

% three points for crazy
r21 = crazy.state_t0(1,1:3);
r22 = crazy.state_t0(20,1:3);
r23 = crazy.state_t0(140,1:3);

P = r22 - r21;
Q = r23 - r21;

n2 = cross(P,Q);

const = (n2(1)*r21(1) + n2(2)*r21(2) + n2(3)*r21(3)); 

% plane equation
P2 = [n2,const];

% Intersect line

m = cross(n1,n2)';

x2 = (P2(4) - (P1(4)*P2(1))/P1(1))/(P2(2) - (P2(1)*P1(2))/P1(1));

x1 = (-P1(2)*x2 + P1(4))/P1(1);

x3 = 0;

b = [x1;x2;x3];



x = linspace(1e-11,1e-10);
y = m*x + b;


y = y';









timespan = [0 T]; % 0 to T/3.175
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time,Circ_state] = ode45(@coast_ODE, timespan, init, options,mu);

figure
h1 = gca;
earth_sphere(h1)
hold on

% Circular Trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Plane Intersection point
plot3(Circ_state_burn(1,1),Circ_state_burn(1,2),Circ_state_burn(1,3),'*','LineWidth',5)

plot3(-3.9080e3,-5.4126e3,-7.9984,'*','LineWidth',5)


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Space Craft Position and Plane Intersection Point');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Circular Transfer Trajectory','Space Craft','Plane Intersection Point','best')





%% Inc change at intersection point


Dang = leo1.inc - crazy.inc;
DeltaVinc = abs(2*norm([5.8964,-4.2534,-2.6167])*sin(Dang/2));



% Circle
h = h3_1;
ecc = 0;
theta = pi;
a = Crazy.rp;
T = (2*pi*a^(3/2))/(sqrt(mu));

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h, theta,ecc);
vperifocal = v_vect_perifocal(mu, h, theta,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(leo1.RAAN,Crazy.inc,leo1.w);

% Q transpose * r_peri = r_eci
r_final = (Qperi2geo' * rperifocal);
v_final =(Qperi2geo' * vperifocal);


timespan = [0 T];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time,Circ_state_2] = ode45(@coast_ODE, timespan, init, options,mu);

% Circle
h = h3_1;
ecc = 0;
theta = pi;
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

[Circ_time,Circ_state_1] = ode45(@coast_ODE, timespan, init, options,mu);

figure
h1 = gca;
earth_sphere(h1)
hold on

% 
plot3(Circ_state_1(:,1),Circ_state_1(:,2),Circ_state_1(:,3),'LineWidth',2)


% Circular Trajectory
plot3(Circ_state_2(:,1),Circ_state_2(:,2),Circ_state_2(:,3),'LineWidth',2)

% Plane Intersection point

plot3(-3.9080e3,-5.4126e3,-7.9984,'*','LineWidth',5)


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Inc Change');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Circular Transfer Trajectory','Circular Transfer Trajectory with Changed Inc','Space Craft at Plane Intersection Point','best')

%% Wait until space craft comes to perigee of crazy



for i = 1:length(crazy.state_t0)
    temp(i) = norm([crazy.state_t0(i,1),crazy.state_t0(i,2),crazy.state_t0(i,3)]);
end

[minTemp,index] = min(temp);

timespan = [0 T/8.2];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time,Circ_state] = ode45(@coast_ODE, timespan, init, options,mu);



figure
h1 = gca;
earth_sphere(h1)
hold on

% Circular Trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Plane Intersection point
plot3(-3.9080e3,-5.4126e3,-7.9984,'*','Color','g','LineWidth',1)

% perigee of object 2
plot3(crazy.state_t0(index,1),crazy.state_t0(index,2),crazy.state_t0(index,3),'*','Color','r','LineWidth',1)



% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Wait Until Perigee Position of Object II');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Circular Transfer Trajectory with Changed Inc','Space Craft at Plane Intersection Point','Perigee Position of Object II','best')



%%


temp = 8.4428 + 0.8443 +  0.7950  + 0.4749 + 1.4405 + 8.0556 + 49.4024

%% Phasing Graph


divisions = 10;

[DeltaV.Phasing,Time.Phasing,hPhasing,TPhasing,eccPhasing,raPhasing] = phasing(Crazy.rp,Crazy.ra,mu,Crazy.theta,divisions);


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

for i = 1:length(crazy.state_t0)
    temp(i) = norm([crazy.state_t0(i,1),crazy.state_t0(i,2),crazy.state_t0(i,3)]);
end

[minTemp,index] = min(temp);


duration = 11.9975 ; % in HOURS (specified by user)

durationsec = duration *3600; % in seconds now

% *IMPORTANT: Don't forget that Me (and thus E) changes with time and can
% be used to calculate theta, True Anomaly

tstart = 0;
tstop = crazy.t0 + durationsec;
tspan = [tstart tstop];
ecc = crazy.ecc;
n = crazy.n;
Me0 = crazy.Me; % from TLEs
state = [crazy.r_epoch';crazy.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees,E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

crazy.state_afterphase = newstate(end,1:6);
crazy.time_aferphase = newtime(end);


figure
h1 = gca;
earth_sphere(h1)
hold on

% Crazy Trajectory
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)

% Phasing Trajectory
plot3(Phasing_state(:,1),Phasing_state(:,2),Phasing_state(:,3),'r','LineWidth',2)

% S/C
plot3(crazy.state_t0(index,1),crazy.state_t0(index,2),crazy.state_t0(index,3),'*','LineWidth',5)

% Object 2
plot3(crazy.state_afterphase(end,1),crazy.state_afterphase(end,2),crazy.state_afterphase(end,3),'*','LineWidth',5) % current location




for i = 1:length(Phasing_state)
    temp(i) = norm([Phasing_state(i,1),Phasing_state(i,2),Phasing_state(i,3)]);
end

minTemp_phasing = min(temp);
% Doesnt run into Earth









title("Phasing Transfer")

% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Phasing Maneuver');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Object II Trajectory','Phasing Trajectroy','Space Craft Location','Object II','best')


% 2D Phasing
figure
hold on

% Crazy Trajectory
plot(crazy.state_t0(:,1),crazy.state_t0(:,2),'b','LineWidth',2)

% Phasing Trajectory
plot(Phasing_state(:,1),Phasing_state(:,2),'r','LineWidth',2)

% S/C
plot(crazy.state_t0(index,1),crazy.state_t0(index,2),'*','LineWidth',5)

% Object 2
plot(crazy.state_afterphase(end,1),crazy.state_afterphase(end,2),'*','LineWidth',5) % current location


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Phasing Maneuver X-Y Plane');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Object II Trajectory','Phasing Trajectroy','Space Craft Location','Object II','best')



%% Crazy and Circle Trajectory with point at Crazy perigee for burn 

% were not really doing a Hohmann but the code solves for stuff we need
[deltaV,~,Vp1,Vp2,deltaV2,h1_2,h2_2,h3_2] = Hohmann(Crazy.rp,Crazy.rp,Crazy.rp,Crazy.ra,mu);

% Graph it 

figure
h1 = gca;
earth_sphere(h1)
hold on

% Crazy Trajectory
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)

% Circular trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Space Craft
plot3(crazy.state_t0(index,1),crazy.state_t0(index,2),crazy.state_t0(index,3),'*','LineWidth',5)









% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('Burn to line up with Object II');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Object II Trajectory','Circular Transfer Trajectory with Changed Inc','Space Craft Location','best')





