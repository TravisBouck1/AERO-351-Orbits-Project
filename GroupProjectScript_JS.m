% Self, Justin
% Fall 2022
% Orbits AERO351
% Group Project

clear all; close all; clc; 

disp("AERO351 Orbital Debris Clean-Up Project; Fall 2022")
disp("Travis Bouck, Justin Self, Anmol Sharma, Shaniya Singh")
disp("---------------------------------------------")
disp("---------------------------------------------")
disp("---------------------------------------------")


%% Finding LEO, MEO, and GEO objects. 
% The strategy: find LEO and MEO RIGHT on the border of the zones. 
% LEO = 2000 km altitude; or 8378 km orbital radius
% MEO = 2000 km --> 36,000 km altitude; or 8378 km --> 42378 km
% GEO = 36,000 km --> or 42378 km and higher

%% OBJECT I: "PSLV R/B" (LEO OBJECT)
% Catalog ID: 37842

rearth = 6378; % km
mu = 398600; % km3/s2
hours = 1 / 3600; % converts from sec to hour
h2seconds = 3600; % converts to hour to sec when multiplied by hours

% LEO1 TLEs

disp("----------------")
disp("----------------")
disp("Object I: LEO -- PSLV R/B")
disp("----------------")

leo1.epoch = "20 November 2022, 06:09:52"; % UTC
leo1.ecc = 0.0059996;
leo1.inc = 19.7942; % deg
leo1.peri = 778; % km (THIS IS ALTITUDE, not orbital radius)
leo1.apo = 864; % km
leo1.RAAN = 293.8926; % deg
leo1.w = 126.0205; % deg
leo1.revperday = 14.21211707; % rev/day
leo1.Me = 234.5776; % deg

leo1.UTC = [2022 11 20 6 09 52];

% Call function to find r,v, T, and a at epoch
[leo1.r_epoch,leo1.v_epoch,leo1.T,leo1.a,leo1.TA] = rvepoch_fromTLEs(leo1.UTC,leo1.ecc,leo1.inc,leo1.peri,leo1.apo,leo1.RAAN,leo1.w,leo1.Me);

disp("LEO1 period is: " + leo1.T * hours + " hours")
disp("HEART CHECK: This is comforting. Circular LEO is ~90 minutes; we are only 19 minutes longer!")

disp("LEO1 true anomaly at time of TLE retrieval is: " + rad2deg(leo1.TA) + " degrees.")

% HEART CHECKS ON OUR WORK
disp("HEART CHECKS SO FAR:")
disp("LEO1 position at epoch is: " + norm(leo1.r_epoch) + " km")
disp("LEO1 altitude at epoch is: " + (norm(leo1.r_epoch) - rearth) + " km")
disp("LEO1 velocity at epoch is: " + norm(leo1.v_epoch) + " km/s")
disp("-----")
disp("These all seem reasonable for a LEO object.")
disp("------------------")



% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%% Finding MISSION START TIME, t0
% Design: mission will start when LEO1 is at apogee, since we will burn out after 5 orbits 
% at apogee.

% mean motion, n. 
leo1.n = 2*pi / leo1.T; % seconds

% find time since perigee
leo1.timesinceperigee = deg2rad(leo1.Me) / leo1.n; % seconds
disp("LEO1 time since perigee is: " + leo1.timesinceperigee*hours + " hours")
% did you convert Me to radians? *****

% Time FORWARD to Apogee
leo1.timeTOapogee = (3/2)*leo1.T - leo1.timesinceperigee; % seconds
disp("LEO1 time forward to apogee is: " + leo1.timeTOapogee * hours + " hours")

leo1.timetot0_sec = leo1.timeTOapogee + 28*leo1.T; % time from TLE retrieval to t0
leo1.timetot0_hours = leo1.timetot0_sec/3600; % in hours

% report mission t0: START MISSION AT LEO1 APOGEE.
disp("Mission will start exactly " + leo1.timetot0_hours + " hours after TLE retrieval.")

missiont0 = "November 22, 2022, 06:52:36"; % <-----------------------------
missiont0_UTC = [06 52 36];

[missionstart.jd, missionstart.UT_decimal, missionstart.j0] = juliandate_justin(2022,11,22,missiont0_UTC);

disp("Mission start date (UTC) is: " + missiont0)
disp("Mission start date (JD) is: " + missionstart.jd + " days")

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
%% Find position and velocity of LEO1 at t0.
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------

% Coast to t0 from starting epoch
tspan = [0 leo1.timetot0_sec]; % seconds
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [leo1.r_epoch, leo1.v_epoch]'; % leo1 at TLE epoch

% call ode here
[leo1.time_TLE_to_t0, leo1.state_t0] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % this should give us r, v at t0

% LOCATION AND VELOCITY OF LEO AT t0 MISSION START TIME
% t0 position
leo1.r_t0 = leo1.state_t0(end,1:3)';
leo1.r_t0MAG = norm(leo1.r_t0);

% t0 velocity
leo1.v_t0 = leo1.state_t0(end,4:6)';
leo1.v_t0MAG = norm(leo1.v_t0);

disp("Location of LEO1 at t0 is: " + leo1.r_t0MAG + " km")
disp("Altitude of LEO1 at t0 is: " + (leo1.r_t0MAG - rearth) + " km")
disp("Altitude of LEO at apogee is: " + leo1.apo + " km")
disp("-----")
disp("Good, we hit apogee for mission start as planned.")
disp("---------------")
disp("Velocity of LEO1 at t0 is: " + leo1.v_t0MAG + " km/s")
disp("HEART CHECK: These are still in LEO!")
disp("---------")

%% OBJECT II: "crazy" - BLOCK DM-SL (free choice object)
% Catalog ID: 38750 

% FIND AND PLOT CURRENT LOCATION.

disp("----------------")
disp("----------------")
disp("Object II: MEO/LEO -- BLOCK DM-SL (free choice object) ")
disp("----------------")

% Epoch info from Heavens-Above
% define function inputs from TLEs

% Epoch (UTC):	21 November 2022 09:40:33
% Eccentricity:	0.7143988
% inclination:	0.6194°
% perigee height:	298 km
% apogee height:	33699 km
% right ascension of ascending node:	9.5320°
% argument of perigee:	134.3946°
% revolutions per day:	2.42895513
% mean anomaly at epoch:	296.8969°
% orbit number at epoch:	8867

crazy.UTC = [2022 11 21 09 40 33]; % yyyy mm dd hh mm ss
crazy.ecc = 0.7143988;
crazy.inc = 0.6194; % deg
crazy.peri = 298; % km
crazy.apo = 33699; % km
crazy.RAAN = 9.5320; % deg
crazy.w = 134.3946; % deg
crazy.Me = 296.8969; % deg

% Call function
[crazy.r_epoch,crazy.v_epoch,crazy.T,crazy.a,crazy.TA] = rvepoch_fromTLEs(crazy.UTC,crazy.ecc,crazy.inc,crazy.peri,crazy.apo,crazy.RAAN,crazy.w,crazy.Me);

% HEART CHECKS ON OUR WORK
disp("----------------------------------")
disp("HEART CHECKS ON Object II at time of epoch:")
disp("Object II period is: " + crazy.T * hours + " hours (good for MEO)")
disp("Object II true anomaly at time of TLE retrieval is: " + rad2deg(crazy.TA) + " degrees.")
disp("Object II position at epoch is: " + norm(crazy.r_epoch) + " km")
disp("Object II altitude at epoch is: " + (norm(crazy.r_epoch) - rearth) + " km")
disp("At time of epoch, velocity is: " + norm(crazy.v_epoch) + " km/s (good; MEOish)")

%% Find r and v vectors for CRAZY object at t0. 
% Coast to t0 from starting epoch

% Epoch (UTC):	            21 November 2022 09:40:33       FROM HERE
% Mission start t0 (UTC):   November 22, 2022, 06:52:36     TO HERE

% Using "timeanddate".com, we find that:

crazy.timetot0_sec = 76323; % sec

tspan = [0 crazy.timetot0_sec]; % seconds
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [crazy.r_epoch, crazy.v_epoch]'; % leo1 at TLE epoch

crazy.n = 2*pi / crazy.T; % seconds

% find time since perigee
crazy.timesinceperigee = deg2rad(crazy.Me) / crazy.n; % seconds
disp("Object II (Crazy) time since perigee is: " + crazy.timesinceperigee*hours + " hours")

% call ode here
[crazy.time_TLE_to_t0, crazy.state_t0] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % this should give us r, v at t0

crazy.t0 = crazy.time_TLE_to_t0(end);

% LOCATION AND VELOCITY OF Object II AT t0 MISSION START TIME
% t0 position
crazy.r_t0 = crazy.state_t0(end,1:3)';
crazy.r_t0MAG = norm(crazy.r_t0);

% t0 velocity
crazy.v_t0 = crazy.state_t0(end,4:6)';
crazy.v_t0MAG = norm(crazy.v_t0);

disp("Location of Object II at t0 is: " + crazy.r_t0MAG + " km")
disp("Velocity of Object II at t0 is: " + crazy.v_t0MAG + " km/s")
disp("HEART CHECK: Location is good; MEO, as expected.")
disp("Velocity is good; MEO speeds!")
disp("---------")

%% OBJECT III: MEO - BREEZE-M R/B 
% Catalog ID: 37934

% FIND AND PLOT CURRENT LOCATION.
disp("----------------")
disp("----------------")
disp("Object II: MEO -- BREEZE-M R/B ")
disp("----------------")

% Epoch info from Heavens-Above
% define function inputs from TLEs

meo1.UTC = [2022 11 20 15 00 57];
meo1.ecc = 0.3883671;
meo1.inc = 1.0827;
meo1.peri = 11647;
meo1.apo = 34538;
meo1.RAAN = 346.4374;
meo1.w = 293.5371;
meo1.Me = 202.1083;

% Call function
[meo1.r_epoch,meo1.v_epoch,meo1.T,meo1.a,meo1.TA] = rvepoch_fromTLEs(meo1.UTC,meo1.ecc,meo1.inc,meo1.peri,meo1.apo,meo1.RAAN,meo1.w,meo1.Me);

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Calculate [r,v] for object at mission start time, t0. 
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Time BACKWARD TO PERIGEE.

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% mean motion, n, from TLEs.
meo1.n = 2*pi / meo1.T; % seconds

% find time since perigee
meo1.timesinceperigee = deg2rad(meo1.Me) / meo1.n; % seconds
% DID YOU CONVERT TO RADIANS?

disp("MEO1 time since perigee is: " + meo1.timesinceperigee*hours + " hours")


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Find new r and v vectors at mission start. 
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Coast to t0 from starting epoch

% Epoch (UTC):	            20 November 2022 15:00:57       FROM HERE
% Mission start t0 (UTC):   November 22, 2022, 06:52:36     TO HERE

meo1.timetot0_sec = 143499; % sec

% Coast to t0 from starting epoch
tspan = [0 meo1.timetot0_sec]; % seconds
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [meo1.r_epoch, meo1.v_epoch]'; % leo1 at TLE epoch

% call ode here
[meo1.time_TLE_to_t0, meo1.state_t0] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % this should give us r, v at t0

meo1.t0 = meo1.time_TLE_to_t0(end);

% LOCATION AND VELOCITY OF MEO AT t0 MISSION START TIME
% t0 position
meo1.r_t0 = meo1.state_t0(end,1:3)';
meo1.r_t0MAG = norm(meo1.r_t0);

% t0 velocity
meo1.v_t0 = meo1.state_t0(end,4:6)';
meo1.v_t0MAG = norm(meo1.v_t0);

disp("Location of MEO1 at t0 is: " + meo1.r_t0MAG + " km")
disp("Velocity of MEO1 at t0 is: " + meo1.v_t0MAG + " km/s")
disp("---------")


%% OBJECT IV: GEO - AGENA D R/B
% Catalog ID: 4511

% FIND AND PLOT CURRENT LOCATION.
disp("----------------")
disp("----------------")
disp("Object IV: GEO -- AGENA D R/B ")
disp("----------------")

% Epoch info from Heavens-Above
% define function inputs from TLEs

% Epoch (UTC):	21 November 2022 01:01:02
% Eccentricity:	0.0030049
% inclination:	1.6751°
% perigee height:	35618 km
% apogee height:	35872 km
% right ascension of ascending node:	277.7484°
% argument of perigee:	92.4831°
% revolutions per day:	1.00418284
% mean anomaly at epoch:	267.9209°
% orbit number at epoch:	1245

geo1.UTC = [2022 11 21 01 01 02];
geo1.ecc = 0.0030049;
geo1.inc = 1.6751;
geo1.peri = 35618;
geo1.apo = 35872;
geo1.RAAN = 277.7484;
geo1.w = 92.4831;
geo1.Me = 267.9209;

% Call function
[geo1.r_epoch,geo1.v_epoch,geo1.T,geo1.a,geo1.TA] = rvepoch_fromTLEs(geo1.UTC,geo1.ecc,geo1.inc,geo1.peri,geo1.apo,geo1.RAAN,geo1.w,geo1.Me);

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Calculate [r,v] for object at mission start time, t0. 

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% Time BACKWARD TO PERIGEE.

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% mean motion, n from TLEs
geo1.n = 2*pi / geo1.T; % seconds

% find time since perigee
geo1.timesinceperigee = deg2rad(geo1.Me) / geo1.n; % seconds
% DID YOU CONVERT TO RADIANS?
disp("GEO1 time since perigee is: " + meo1.timesinceperigee*hours + " hours")

% Coast to t0 from starting epoch

% Epoch (UTC):	            21 November 2022 01:01:02       FROM HERE
% Mission start t0 (UTC):   November 22, 2022, 06:52:36     TO HERE

geo1.timetot0_sec = 107494; % sec

% Coast to t0 from starting epoch
tspan = [0 geo1.timetot0_sec]; % seconds
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [geo1.r_epoch, geo1.v_epoch]'; % leo1 at TLE epoch

% call ode here
[geo1.time_TLE_to_t0, geo1.state_t0] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % this should give us r, v at t0

geo1.t0 = geo1.time_TLE_to_t0(end);

% LOCATION AND VELOCITY OF GEO AT t0 MISSION START TIME
% t0 position
geo1.r_t0 = geo1.state_t0(end,1:3)';
geo1.r_t0MAG = norm(geo1.r_t0);

% t0 velocity
geo1.v_t0 = geo1.state_t0(end,4:6)';
geo1.v_t0MAG = norm(geo1.v_t0);

disp("---------")
disp("Location of GEO1 at t0 is: " + geo1.r_t0MAG + " km")
disp("Velocity of GEO1 at t0 is: " + geo1.v_t0MAG + " km/s")
disp("GEO period is: " + geo1.T/3600 + " hours.")

disp("HEART CHECKS FOR GEO at t0:")
disp("Period is right on. Altitude is good. Velocity is right on. Great! ")


%% Plots of all objects at Mission Time (t) = t0.

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % Object I: LEO1 ORBIT AT MISSION START TIME, t0
       plot3(leo1.state_t0(:,1),leo1.state_t0(:,2),leo1.state_t0(:,3),'r','LineWidth',2)
       plot3(leo1.state_t0(end,1),leo1.state_t0(end,2),leo1.state_t0(end,3),'*','LineWidth',5)
      
   % Object II: CRAZY ORBIT AT MISSION START TIME, t0
       plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2)
       plot3(crazy.state_t0(end,1),crazy.state_t0(end,2),crazy.state_t0(end,3),'*','LineWidth',5)
       
   % Object III: MEO1 ORBIT AT MISSION START TIME, t0
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2)
       plot3(meo1.state_t0(end,1),meo1.state_t0(end,2),meo1.state_t0(end,3),'*','LineWidth',5)

   % Object IV: GEO1 ORBIT AT MISSION START TIME, t0
       plot3(geo1.state_t0(:,1),geo1.state_t0(:,2),geo1.state_t0(:,3),'k','LineWidth',2)
       plot3(geo1.state_t0(end,1),geo1.state_t0(end,2),geo1.state_t0(end,3),'*','LineWidth',5)
       
        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title('Space Debris Targets at Mission Start time, {t_0}'); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 10) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','LEO1 orbit','LEO1','Object II orbit', 'Object II', 'MEO1 orbit', 'MEO1','GEO1 Orbit','GEO1','location','best')
       
%        % annotation
%        ob3 = text(crazy.r_epoch(1),crazy.r_epoch(2),crazy.r_epoch(3),'MEO1');
%        set(ob3,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','right') 


%% Calculating new Me, E, and TA for Object II (crazy) orbit after specified time interval
% Time interval: [start(t0) end(duration requested)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%duration = 10.871124728243107; % in HOURS (specified by user) TRAVIS LIKES
%THIS NUMBER A LOT
duration = 69.4555; % in HOURS (specified by user)

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

disp("----------------------")
disp("----------------------")
disp("----------------------")

disp("**Object II stats after " + duration + " hours:")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Altitude is: " + (norm(crazy.state_afterphase(1:3)) - rearth) + " km")
disp("Velocity is: " + norm(crazy.state_afterphase(4:6)) + " km/s")

disp("----------------------")
disp("----------------------")

%% Find location of MEO1 at the moment our s/c burns out from Object II. 

% *IMPORTANT: Don't forget that Me (and thus E) changes with time and can
% be used to calculate theta, True Anomaly

tstart = 0;
tstop = meo1.t0 + durationsec;
tspan = [tstart tstop];
ecc = meo1.ecc;
n = meo1.n;
Me0 = meo1.Me; % from TLEs
state = [meo1.r_epoch';meo1.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

meo1.state_afterphase = newstate(end,1:6);
meo1.time_aferphase = newtime(end);

disp("-------------------------")
disp("-------------------------")
disp("MEO1 after " + duration + " hours:")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(meo1.state_afterphase(1:3)) + " km")
disp("Velocity is: " + norm(meo1.state_afterphase(4:6)) + " km/s")

%% Plot MEO and Crazy at the moment our s/c burns away from Crazy

% Location of CRAZY and our s/c at end of Crazy part of mission
crazy.final_r_p = [-5.395919869552229e+03; -3.931015381482056e+03; 11.951231023118638]; % km
crazy.final_rp_MAG = norm(crazy.final_r_p);
disp("Check for RP crazy (should be 298 km) is: " + (crazy.final_rp_MAG - rearth) + " km")

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % Object II "Crazy" at mission time = 69.6323 hours (beginning of Mission Phase 3 (MP3)
       plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2) % orbit traj
       plot3(crazy.state_afterphase(end,1),crazy.state_afterphase(end,2),crazy.state_afterphase(end,3),'*','LineWidth',5) % current location

   % Object III "MEO1" at mission time = 69.6323 hours (beginning of MP3)
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) % orbit traj
       plot3(meo1.state_afterphase(end,1),meo1.state_afterphase(end,2),meo1.state_afterphase(end,3),'*','LineWidth',5) % current location

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Orbits and locations of S/C (and Object II) and MEO1 at end of Mission Phase 2.']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 8) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','Object II Orbit','Object II (and spacecraft) position', 'MEO1 Orbit','MEO1 location','location','best')

%% Mission Phase I time calculation. 
%{
 We need to know how much time has passed since the end of mission phase 1
(5 orbital periods with LEO1) and the beginning first burn location for
mission phase 2 (after apogee in LEO1 orbit, our s/c needs to travel to
a new location on the LEO orbit such that this new point is on the line
of intersection between the LEO and OBJECT II orbital planes.
%}

% Mission Phase 1: 5 Orbital periods with LEO1. 
phase1.r_t0 = leo1.r_t0;
phase1.v_t0 = leo1.v_t0;
t0 = 0;

% Propogate s/c forward 5 periods. 

% Coast to end of Phase 1 (5 periods) from t0 (mission start)
tspan = [t0 5*leo1.T]; % seconds
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state = [phase1.r_t0; phase1.v_t0]; % leo1 at TLE epoch

% call ode here
[phase1.final_time, phase1.final_state] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % gives us time, r, v at end of mission

% Start tracking mission elapsed time.
% RECORD ALL THESE IN SECONDS! 
% SECONDS
% SECONDS
missiontime.endp1 = phase1.final_time(end); % IN SECONDS

disp("~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~")
disp("MISSION PHASE 1.")
disp("~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~")
disp("After five orbital periods with the LEO1 object, our spacecraft is located at (in km): ")
disp(phase1.final_state(end,1:3))
disp("Velocity (km/s) is: ")
disp(phase1.final_state(end,4:6))
disp("Mission Elapsed Time (MET) is " + hours*missiontime.endp1 + " hours")

disp("Heart check these: t0 radius and velocity match these. Good.")

%% MISSION PHASE 2: S/c is at LEO1 apogee; just finished 5 periods with object. 
% Now, we need to calculate the time to travel between apogee and burn 1
% location. 

phase2.initialstate = phase1.final_state(end,1:6)'; % r, v vectors (col vect)
phase2.initialtime = missiontime.endp1; % in seconds

% COAST TO PERIGEE
t1 = phase2.initialtime; % seconds
t2 = 0.5*leo1.T; % time to go from apogee to perigee = half a period
tspan = [t1 t2]; % seconds
state = phase2.initialstate; % initial conditions

% call ode here
[phase2.perigee1time, phase2.perigee1state] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % gives us time, r, v at end of mission
% this checks out fine; we are indeed at perigee; good. 

% Mission Elapsed time
missiontime.phase2.endwait = t1 + t2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
disp("MISSION ELAPSED TIME (MET)")
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
% Travis' code did this next part

phase2.circtime = 0.7950*3600; % hours to seconds.

% Mission Elapsed time
missiontime.phase2.aftercirc = missiontime.phase2.endwait + phase2.circtime; % in sec

disp("MET after circularization burn is: " + (hours*missiontime.phase2.aftercirc) + " hours.")

% Mission elapsed time after Travis's wait time from circular orbit
% location to node line between object 2 plane and circular plane

phase2.waittime_planes = 0.4749 * 3600; % hours to seconds.
missiontime.phase2.afterwait_planes = missiontime.phase2.aftercirc + phase2.waittime_planes; % in seconds
disp("MET after wait time for planes to line up (circ and Ob 2) is: " + hours* missiontime.phase2.afterwait_planes + " hours.")

% Next, the s/c will wait until perigee for Object 2 burn. Only need to
% change speed, since perigee Ob2 = perigee circular parking orbit!

% Mission Elapsed time just after wait time in parking orbit
phase2waittime_crazyburn = 1.4405*3600; % hours to seconds; got this from Travis

missiontime.phase2.afterwait_crazyburn = missiontime.phase2.afterwait_planes + phase2waittime_crazyburn;
disp("MET after parking orbit just before burn into Object II orbit is: " + hours*missiontime.phase2.afterwait_crazyburn + " hours")

% Good, this agrees with Travis' work so far. 

% We are now on Crazy (Object II) orbit!!! Phasing to catch up with the
% object is next. 

% PHASING MANEUVER
phase2.phasingmaneuver_time = 8.0556*3600; % hours to seconds

% Mission Elapsed time AFTER PHASING MANEUVER
missiontime.phase2.phasing = missiontime.phase2.afterwait_crazyburn + phase2.phasingmaneuver_time; % seconds
disp("MET after phasing maneuver to catch up to Object II is: " + hours*missiontime.phase2.phasing + " hours")

% MAIN MISSION OBJECTIVE 2: FLY WITH OBJECT II FOR 5 PERIODS.
phase2.fiveperiods = 5*crazy.T; % seconds
disp("MISSION OBJECTIVE 2: FLY FIVE ORBITS WITH OBJECT II.")
disp("Five orbits with Object II takes " + phase2.fiveperiods*hours + " hours")

% Mission Elapsed time AFTER 5 ORBITS with OBJECT 2
missiontime.phase2.fiveperiods = missiontime.phase2.phasing + phase2.fiveperiods; % seconds
disp("MET at the end of Phase 2 (after 5 periods with Object 2) is: " + hours*missiontime.phase2.fiveperiods + " hours")

%% Phase 2 complete.
disp("~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~")
disp("Mission Objective #2 Complete: Five orbits with Object II.")
disp("Total Mission Elapsed Time: " + hours*missiontime.phase2.fiveperiods + " hours")
disp("Total Delta-V used so far: 5.0334 km/s")
disp("Next Objective: Phase III: Rendezvous with MEO object.")
disp(" ")

%% MISSION PHASE 3: MEO object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait on object 2 orbit for 3/4 period to line up for Lambert's burn to
% MEO1.
k = 1; % for k = 0, this is the end of Mission Phase 2. k = Ob2 periods.
% DO NOT CHANGE THIS VALUE FROM k = 1.0

k1 = 4/5; % Based on the trajectory plots, this seems a wise location on the orbit to burn from.
k2 = 1.0; % Based on the trajectory plots, this seems a reasonable target location for the MEO rendezvous
phase3.deltat_lamberts = (k2 - k1)*crazy.T; % in seconds
disp("Deltat Lamberts for Phase III is: " + phase3.deltat_lamberts + " seconds")

ob2_wait_time = k*crazy.T /3600; % lining up orbits for LAMBERT'S; this is in hours
duration = hours*missiontime.phase2.fiveperiods  + ob2_wait_time; % in HOURS (specified by user)

durationsec = duration *3600; % in seconds now

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes to Mission Elapsed Time
missiontime.phase3_prelamberts = k1*crazy.T + missiontime.phase2.fiveperiods; % in seconds
missiontime.phase3_postlamberts = k2*crazy.T + missiontime.phase2.fiveperiods; % in seconds

% this is displayed below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

disp("----------------------")
disp("----------------------")
disp("----------------------")

disp("**Object II stats after s/c has travelled " + k + " orbits")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Distance is: " + (norm(crazy.state_afterphase(1:3)) - rearth) + " km")
disp("Speed is: " + norm(crazy.state_afterphase(4:6)) + " km/s")

disp("-------")
disp("-------")

disp("Position vector (km) is: ")
disp(crazy.state_afterphase(1)) % position x
disp(crazy.state_afterphase(2)) % position y
disp(crazy.state_afterphase(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(crazy.state_afterphase(4)) % veloc x
disp(crazy.state_afterphase(5)) % veloc y
disp(crazy.state_afterphase(6)) % veloc z

disp("----------------------")
disp("----------------------")

%% Find location of MEO1 at the moment our s/c burns out from Object II. 

% *IMPORTANT: Don't forget that Me (and thus E) changes with time and can
% be used to calculate theta, True Anomaly

tstart = 0;
tstop = meo1.t0 + durationsec;
tspan = [tstart tstop];
ecc = meo1.ecc;
n = meo1.n;
Me0 = meo1.Me; % from TLEs
state = [meo1.r_epoch';meo1.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

meo1.state_afterphase = newstate(end,1:6);
meo1.time_aferphase = newtime(end);

disp("-------------------------")
disp("-------------------------")
disp("MEO1 after s/c has travelled " + k + " orbits on Object II trajectory")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(meo1.state_afterphase(1:3)) + " km")
disp("Speed is: " + norm(meo1.state_afterphase(4:6)) + " km/s")

disp("-------")
disp("-------")

disp("Position vector (km) is: ")
disp(meo1.state_afterphase(1)) % position x
disp(meo1.state_afterphase(2)) % position y
disp(meo1.state_afterphase(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(meo1.state_afterphase(4)) % veloc x
disp(meo1.state_afterphase(5)) % veloc y
disp(meo1.state_afterphase(6)) % veloc z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot it

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % Object II "Crazy" at mission time = 69.6323 hours (beginning of Mission Phase 3 (MP3)
       plot3(crazy.state_t0(:,1),crazy.state_t0(:,2),crazy.state_t0(:,3),'b','LineWidth',2) % orbit traj
       plot3(crazy.state_afterphase(end,1),crazy.state_afterphase(end,2),crazy.state_afterphase(end,3),'*','LineWidth',2) % current location

   % Object III "MEO1" at mission time = 69.6323 hours (beginning of MP3)
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) % orbit traj
       plot3(meo1.state_afterphase(end,1),meo1.state_afterphase(end,2),meo1.state_afterphase(end,3),'*','LineWidth',2) % current location

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Orbits and locations of S/C upon rendezvous with MEO1 after Lamberts Burn 2 (Arrival).']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 8) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','Object II Orbit','Object II position', 'MEO1 Orbit','S/C and MEO1 rendezvous','location','best')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Total Mission Elapsed Time (Phase III Pre-Lambert's): " + hours*missiontime.phase3_prelamberts + " hours")
disp("Total Mission Elapsed Time (Phase III Post-Lambert's): " + hours*missiontime.phase3_postlamberts + " hours")
disp("Lambert's transfer took: " + abs(missiontime.phase3_prelamberts - missiontime.phase3_postlamberts)*hours + " hours")
disp("Total Delta-V for this Lambert's Transfer is: 3.2360 km/s") % Shaniya's code

% Final r/v vectors of our s/c once we rendezvous with MEO1. 
phase3.postlamberts_state = meo1.state_afterphase(end,1:6);

% this is our s/c after rendezvous

%% Now, MAIN MISSION OBJECTIVE: 5 PERIODS WITH MEO1 <<-----------

tstart = missiontime.phase3_postlamberts;
tstop = missiontime.phase3_postlamberts + 5*meo1.T; % in seconds
tspan = [tstart tstop];

state = phase3.postlamberts_state';

% Call function that does the work
[newtime,newstate] = ode45(@non_impulsive_COAST,tspan,state,options,mu);

phase3.after5periods_state = newstate(end,1:6);
phase3.after5periods_time = newtime(end);

disp("-------------------------")
disp("-------------------------")
disp("State of S/C and MEO1 after 5 orbital periods:")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(phase3.after5periods_state(1:3)) + " km")
disp("Speed is: " + norm(phase3.after5periods_state(4:6)) + " km/s")

disp("-------")
disp("-------")

disp("Position vector (km) is: ")
disp(phase3.after5periods_state(1)) % position x
disp(phase3.after5periods_state(2)) % position y
disp(phase3.after5periods_state(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(phase3.after5periods_state(4)) % veloc x
disp(phase3.after5periods_state(5)) % veloc y
disp(phase3.after5periods_state(6)) % veloc z

disp("---")
disp("Mission Elapsed Time once 5 periods with MEO1 are complete is: " + hours*phase3.after5periods_time + " hours")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Preparation for Mission Phase 4. 
% Now that we have completed Phase 3, we need to determine the optimal
% transfer trajectory from MEO to GEO. 

% Here we plot the positions of both our S/C(and MEO) with GEO at the exact moment 
% Mission Phase 3 is concluded for perspective.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *IMPORTANT: Don't forget that Me (and thus E) changes with time and can
% be used to calculate theta, True Anomaly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FOR MEO - GEO LAMBERTS TRANSFER %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do we need to wait (park) on the MEO trajectory for some period of time
% before burning out to GEO? Let's find out. 

currenttime = phase3.after5periods_time; % This is current mission elapsed time

k = 0; % for k = 0, this is the end of Mission Phase 2. k = Ob2 periods.

k1 = 0; % Based on the trajectory plots, this seems a wise location on the orbit to burn from.
k2 = 0; % Based on the trajectory plots, this seems a reasonable target location for the GEO rendezvous

phase4.geo_wait_time = k*meo1.T; % lining up orbits for LAMBERT'S; in SECONDS HERE
duration = currenttime*hours + phase4.geo_wait_time*hours; % in HOURS (specified by user)

durationsec = duration *3600; % in seconds now for calculation

tstart = 0; % propgate orbit from GEO TLE retrieval
tstop = geo1.t0 + durationsec;
tspan = [tstart tstop];
ecc = geo1.ecc;
n = geo1.n;
Me0 = geo1.Me; % from TLEs
state = [geo1.r_epoch';geo1.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

phase4.start_geo1_state = newstate(end,1:6);
phase4.start_geo1.time = newtime(end);

disp("-------------------------")
disp("-------------------------")
disp("GEO1 position at MET = " + currenttime*hours + "hours")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(phase4.start_geo1_state(1:3)) + " km")
disp("Speed is: " + norm(phase4.start_geo1_state(4:6)) + " km/s")

disp("-------")
disp("-------")

disp("Position vector (km) is: ")
disp(phase4.start_geo1_state(1)) % position x
disp(phase4.start_geo1_state(2)) % position y
disp(phase4.start_geo1_state(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(phase4.start_geo1_state(4)) % veloc x
disp(phase4.start_geo1_state(5)) % veloc y
disp(phase4.start_geo1_state(6)) % veloc z


% Plot S/C and GEO at the moment Phase 3 ends.

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % MEO1 and S/C orbit after 5 orbital periods
       % orbit traj
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) 
       % current location
       plot3(phase3.after5periods_state(end,1),phase3.after5periods_state(end,2),phase3.after5periods_state(end,3),'*','LineWidth',2) 

    % GEO1 and S/C orbit after 5 orbital periods
       % orbit traj
       plot3(geo1.state_t0(:,1),geo1.state_t0(:,2),geo1.state_t0(:,3),'k','LineWidth',2) 
       % current location of GEO
       plot3(phase4.start_geo1_state(end,1),phase4.start_geo1_state(end,2),phase4.start_geo1_state(end,3),'*','LineWidth',2) 

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Orbit and location of S/C and MEO1 after just completing 5 orbital periods']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 9) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','MEO1 Orbital Trajectory','S/C and MEO1 position after 5 periods','GEO1 Orbital Trajectory','GEO1 Position','location','best')

disp("This is not an optimal location for each object for a transfer, since this would require a major direction change.")
disp("We shall propgate the orbit forward in time until an opportune path exists.")

%% MISSION PHASE 4: Wait on MEO orbit until optimal flight path opens up to transfer to GEO. 

currenttime = phase3.after5periods_time; % This is current mission elapsed time in SECONDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 4.6; % for k = 0, this is the end of Mission Phase 2. k = Ob2 periods. <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHANGE THIS TO OPTIMIZE LAMBERTS Ph4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these k1 and k2 values are just for bookkeeping for the Lambert's code.
k1 = 4.35; % Based on the trajectory plots, this seems a wise location on the orbit to burn from.
k2 = 4.6; % Based on the trajectory plots, this seems a reasonable target location for the GEO rendezvous

phase4.deltaT_lamberts = (k2 - k1)*meo1.T;

phase4.geo_wait_time = k*meo1.T; % lining up orbits for LAMBERT'S; in SECONDS HERE
duration = currenttime*hours + phase4.geo_wait_time*hours; % in HOURS (specified by user)

durationsec = duration *3600; % in seconds now for calculation

tstart = 0; % propgate orbit from GEO TLE retrieval
tstop = meo1.t0 + durationsec;
tspan = [tstart tstop];
ecc = meo1.ecc;
n = meo1.n;
Me0 = meo1.Me; % from TLEs
state = [meo1.r_epoch';meo1.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

phase4.waiting_meo1_state = newstate(end,1:6);
phase4.waiting_meo1.time = newtime(end); % current mission time after waiting period

disp("-------------------------")
disp("-------------------------")
disp("MEO1 position at MET = " + duration + "hours")
% disp("TA is: " + ta_degrees + " degrees")
% disp("Me is: " + Me_degrees + " degrees")
% disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(phase4.waiting_meo1_state(1:3)) + " km")
disp("Speed is: " + norm(phase4.waiting_meo1_state(4:6)) + " km/s")

disp("-------")
disp("-------")

disp("Position vector (km) is: ")
disp(phase4.waiting_meo1_state(1)) % position x
disp(phase4.waiting_meo1_state(2)) % position y
disp(phase4.waiting_meo1_state(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(phase4.waiting_meo1_state(4)) % veloc x
disp(phase4.waiting_meo1_state(5)) % veloc y
disp(phase4.waiting_meo1_state(6)) % veloc z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROPOGATE GEO FORWARD SAME TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tstart = 0; % propgate orbit from GEO TLE retrieval
tstop = geo1.t0 + durationsec;
tspan = [tstart tstop];
ecc = geo1.ecc;
n = geo1.n;
Me0 = geo1.Me; % from TLEs
state = [geo1.r_epoch';geo1.v_epoch'];

% Call function that does the work
[ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0);

phase4.waiting_geo1_state = newstate(end,1:6);
phase4.waiting_geo1.time = newtime(end); % current mission time after waiting period
missiontime.phase4_prelamberts = phase4.waiting_geo1.time; % current mission time after waiting period.

% disp("-------------------------")
% disp("-------------------------")
disp("GEO1 position at MET = " + duration + "hours")
% disp("TA is: " + ta_degrees + " degrees")
% disp("Me is: " + Me_degrees + " degrees")
% disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(phase4.waiting_geo1_state(1:3)) + " km")
disp("Speed is: " + norm(phase4.waiting_geo1_state(4:6)) + " km/s")
% 
% disp("-------")
% disp("-------")

disp("Position vector (km) is: ")
disp(phase4.waiting_geo1_state(1)) % position x
disp(phase4.waiting_geo1_state(2)) % position y
disp(phase4.waiting_geo1_state(3)) % position z

disp("Velocity vector (km/s) is: ")
disp(phase4.waiting_geo1_state(4)) % veloc x
disp(phase4.waiting_geo1_state(5)) % veloc y
disp(phase4.waiting_geo1_state(6)) % veloc z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot S/C and GEO at the moment Phase 3 ends.

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on

   % MEO1 and S/C orbit after GEO WAIT TIME
       % orbit traj
       plot3(meo1.state_t0(:,1),meo1.state_t0(:,2),meo1.state_t0(:,3),'g','LineWidth',2) 
       % current location
       plot3(phase4.waiting_meo1_state(end,1),phase4.waiting_meo1_state(end,2),phase4.waiting_meo1_state(end,3),'*','LineWidth',2) 

    % GEO1 after GEO WAIT TIME
       % orbit traj
       plot3(geo1.state_t0(:,1),geo1.state_t0(:,2),geo1.state_t0(:,3),'k','LineWidth',2) 
       % current location of GEO
       plot3(phase4.waiting_geo1_state(end,1),phase4.waiting_geo1_state(end,2),phase4.waiting_geo1_state(end,3),'*','LineWidth',2) 

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Orbit and location of S/C and MEO1 after waiting on MEO orbit']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 9) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','MEO1 Orbital Trajectory','S/C and MEO1 position after wait time','GEO Orbit','GEO location after waiting','location','best')

disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("GEO wait time is: " + hours*phase4.geo_wait_time + " hours.")
disp("Mission Elapsed Time after WAITING PERIOD (pre MEO -> GEO transfer) " + hours*phase4.waiting_geo1.time + " hours")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")

% We have determined the Phase IV Lambert's burn. 
% Update mission times. 

% Need
% i. MET after Lambert's burn = MET preburn + delta-T burn
% ii. MET after 5 periods with GEO = final MET! 

missiontime.pre_geo5 = missiontime.phase4_prelamberts + phase4.deltaT_lamberts;

disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("Lambert's Phase IV burn time is: " + hours*phase4.deltaT_lamberts + " hours.")
disp("Mission Elapsed Time after Phase IV Lambert's transfer) " + hours*missiontime.pre_geo5 + " hours")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")

%% FINAL PHASE: Five periods with GEO to mission completion.

tstart = missiontime.pre_geo5; % seconds
tstop = missiontime.pre_geo5 +  5*geo1.T; % in seconds
tspan = [tstart tstop];

state = phase4.waiting_geo1_state';

% Call function that does the work
[newtime,newstate] = ode45(@non_impulsive_COAST,tspan,state,options,mu);

phase4.final_state = newstate(end,1:6);
phase4.final_time = newtime(end);

disp("-------------------------")
disp("-------------------------")
disp("State of S/C and GEO1 after 5 orbital periods:")
disp("TA is: " + ta_degrees + " degrees")
disp("Me is: " + Me_degrees + " degrees")
disp("E is: " + E_degrees + " degrees")
disp("Position is: " + norm(phase4.final_state(1:3)) + " km")
disp("Speed is: " + norm(phase4.final_state(4:6)) + " km/s")

disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("TOTAL Mission Elapsed Time once 5 periods with GEO1 are complete is: " + hours*phase4.final_time + " hours")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~~~~")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot S/C // GEO after 5 orbital periods with GEO. End of mission!

 figure
   h1 = gca;
   earth_sphere(h1)
   hold on


    % GEO end of mission
       % orbit traj
       plot3(geo1.state_t0(:,1),geo1.state_t0(:,2),geo1.state_t0(:,3),'k','LineWidth',2) 
       % current location of GEO
       plot3(phase4.final_state(end,1),phase4.final_state(end,2),phase4.final_state(end,3),'*','LineWidth',2) 

        % Graph pretty 
        ylim padded 
        xlim padded
        zlim padded
        xLab = xlabel('{\itx [km]}'); 
        yLab = ylabel('{\ity [km]}'); 
        zLab = zlabel('{\itz [km]}'); 
        plotTitle = title(['Orbit and location of S/C and GEO at end of mission']); 
        set(plotTitle,'FontSize',14,'FontWeight','bold') 
        set(gca,'FontName','Palatino Linotype') 
        set([xLab, yLab],'FontName','Palatino Linotype') 
        set(gca,'FontSize', 9) 
        set([xLab, yLab, zLab],'FontSize', 14) 
        grid on 
        legend('Earth','GEO1 Orbital Trajectory','S/C and GEO1 position after 5 periods','location','best')

        % that's it!

