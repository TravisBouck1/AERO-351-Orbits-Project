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
rp3 = 411.7981 + radiusEarth;
ra3 = rp3;

[DeltaV.Hohmann1, Time.Hohmann1,Vp1_1,Vp2_1,deltaV2_1,h1_1,h2_1,h3_1] = Hohmann(ra1,rp1,rp3,ra3,mu);



%% Hohmann Graph
Crazy.rp = 411.7981 + radiusEarth;
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

Circ_state_burn(:,1) = Circ_state_burn(:,1) + (-3080.17 + 3287.72) + (3248.98-3314.7);
Circ_state_burn(:,2) = Circ_state_burn(:,2) + (5549.06 -5956.57) + (-5980.12 + 6066.79);
Circ_state_burn(:,3) = Circ_state_burn(:,3) + (2071.1 - 2217.71) + (-2211.47 + 2248.99);


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
plot3(Circ_state_burn(end,1),Circ_state_burn(end,2),Circ_state_burn(end,3),'*','LineWidth',5, 'Color', '#EDB120')

% Hohmann Burn 2
plot3(Circ_state_burn(1,1),Circ_state_burn(1,2),Circ_state_burn(1,3),'*','LineWidth',5, 'Color', '#D95319')




% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('2.2 Hohmann Transfer');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','LEO1 orbit','Final Circular Trajectory','Hohmann Transfer Trajectory','Burn #1','Burn #2','best')


Crazy.rp = 6676;



%%

figure
h1 = gca;
earth_sphere(h1)
hold on

% Circularise Orbit- 
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Space Craft Location 1
plot3(Circ_state_burn(1,1),Circ_state_burn(1,2),Circ_state_burn(1,3),'*','LineWidth',5)

% Space Craft Location 2
plot3(-4144.8,-5377.91,-8.65679,'*','LineWidth',5)

% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('2.3 Waiting For Intersection Point');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Final Circular Trajectory','Space Craft Location','Intersection Point','best')

%%




% Circle
h = h3_1;
ecc = 0;
theta = pi;
a = 411.7981+radiusEarth;
T = (2*pi*a^(3/2))/(sqrt(mu));

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h, theta,ecc);
vperifocal = v_vect_perifocal(mu, h, theta,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(leo1.RAAN,leo1.inc,leo1.w);

% Q transpose * r_peri = r_eci
r_final = (Qperi2geo' * rperifocal);
v_final =(Qperi2geo' * vperifocal);


timespan = [0 T/6.05];
options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
init = [r_final; v_final];

[Circ_time_asdf,Circ_state_asdf] = ode45(@coast_ODE, timespan, init, options,mu);








figure
hold on

% Circularise Orbit- 
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Space Craft Location 2
plot3(-4144.8,-5377.91,-8.65679,'*','LineWidth',1)

% Matching point
plot3(Circ_state_asdf(end,1),Circ_state_asdf(end,2),Circ_state_asdf(end,3),'*','LineWidth',1)

legend("circ","set point","match point")

%% Inc change

% Circle Speed
r = 411.7981+radiusEarth;
h = sqrt(mu*r);
Vaz1 = h/r;
Vr1 = 0;

% Crazy speed
h = sqrt(Crazy.rp*mu*(1+Crazy.ecc));
theta = acos((((h^2)/(r*mu))-1)/Crazy.ecc);

Vaz2 = h/r;
Vr2 = (mu/h)*Crazy.ecc*sin(theta);

dihedralAngle = leo1.inc - Crazy.inc;
MassiveOne = sqrt((Vr2-Vr1)^2 + Vaz1^2 + Vaz2^2 - 2*Vaz1*Vaz2*cos(dihedralAngle)); % inc and burn to get on crazy orbit

disp("big DeltaV: " + MassiveOne + " km/s")



figure
h1 = gca;
earth_sphere(h1)
hold on

% Circular Trajectory
plot3(Circ_state(:,1),Circ_state(:,2),Circ_state(:,3),'c','LineWidth',2)

% Object
plot3(-4144.8,-5377.91,-8.65679,'*','LineWidth',1)

% Object II
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)


% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('2.4 Inc and Speed Change');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Final Circular Trajectory','Space Craft','Object II Trajectory','best')

%%





figure
h1 = gca;
earth_sphere(h1)
hold on

% Crazy Trajectory
plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)

% Phasing Trajectory
plot3(Phasing_state(:,1),Phasing_state(:,2),Phasing_state(:,3),'r','LineWidth',2)

%
for i = 1:length(crazy.state_t0)
    temp(i) = norm([crazy.state_t0(i,1),crazy.state_t0(i,2),crazy.state_t0(i,3)]);
end

[minTemp,index] = min(temp);


% Space Craft Trajectory
plot3(crazy.state_t0(index,1),crazy.state_t0(index,2),crazy.state_t0(index,3),'*','Color','r','LineWidth',5)


% Object 22
plot3(crazy.state_t0(index+ 37,1),crazy.state_t0(index+ 37,2),crazy.state_t0(index+ 37,3),'*','Color','r','LineWidth',5)




title("2.5 Phasing Transfer")

% Graph pretty
ylim padded
xlim padded
zlim padded
xLab = xlabel('{\itx [km]}');
yLab = ylabel('{\ity [km]}');
zLab = zlabel('{\itz [km]}');
plotTitle = title('2.5 Phasing Maneuver');
set(plotTitle,'FontSize',14,'FontWeight','bold')
set(gca,'FontName','Palatino Linotype')
set([xLab, yLab],'FontName','Palatino Linotype')
set(gca,'FontSize', 8)
set([xLab, yLab, zLab],'FontSize', 14)
grid on
legend('Earth','Object II Trajectory','Phasing Trajectroy','Space Craft Location','Object II','best')










