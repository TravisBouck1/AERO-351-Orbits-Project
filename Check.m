%% Checks for Aro 351 Project
clc; clear all; close all;
% Check for Lambert's Problem
% example 5.2
% r1 = [5000 10000 2100]; % km
% r2 = [-14600 2500 7000]; % km
% del_t = 3600; % sec
% string = 1;
% tm = 0;
% mu = 398600;
% [v1, v2] = lamberts1(r1, r2, del_t,tm, mu)
% [v1, v2] = lambert_ASp(r1, r2, del_t, string)

% % Problem 5.6
% r1 = [3600 4600 3600];
% r2 = [-5500 6240 5200];
% del_t = 20*3600;
% string = 1;
% [v1, v2] = lambert_ASp(r1, r2, del_t, string)
% norm(v1)
% norm(v2)

%% MEO Final Position
% What do I need to complete the debri collection?
% Do I need plane change? No, Lambert's take care of that
% Initial position vector of the last orbit - good
% Position vector of the GEO object at the time where MEO object ends its 5 periods
% k = 4.35 & 4.6
% s1 = [-2.076365594621423e+04,-3.311926067976723e+04,-1.098068723150617e+02,1.854403312935766,-1.851364048199687,0.046102973984870];
% s2 = [2.832193658123055e+04,-3.103146710543583e+04,-8.667956771495835e+02,2.272395360666405,2.083358235065258,-0.063752690948722];

% k = 4.8 & 4.9
% s1 = [2.510109373371945e+04,-1.344224270032110e+04,5.363665048720463e+02,0.357711778980122,3.786651897006075,-0.022380516588993];
% s2 = [4.037242363321498e+04,1.155249679821297e+04,-1.164919366044516e+03,-0.843844403762679,2.966678272065547,0.028413292792097];

% k = 4.2 & 4.4;
% s1 = [-2.919115534741421e+04,-1.313332685807126e+04,-4.066633999925554e+02,0.110379959283146,-3.371154563565643,0.027355119747521];
% s2 = [97.393094313320260,-4.209928506935755e+04,-56.185424104511770,3.076510943661527,0.016330006685397,-0.089865372505747];

%s1 = [2.510109373371945e+04,-1.344224270032110e+04,5.363665048720463e+02,0.357711778980122,3.786651897006075,-0.022380516588993];
%s2 = [-7.609020648844741e+03,4.147372775119805e+04,2.748591905925336e+02,-3.022716103177036,-0.546351140314283,0.087622116777757];

% k = 4.35, 4.6
s1 = [-1.473446403861598e+04,-3.762629130736880e+04,28.673587265657297,2.189260581538516,-1.186591955263678,0.046887766408053];
s2 = [3.433878743797304e+04,-2.417490245087891e+04,-1.033902024430127e+03,1.771289141948836,2.524834624894050,-0.048552577088648];

meo_r1 = [s1(1), s1(2), s1(3)]; % Position where MEO debri finishes 5 orbits
meo_v1 = [s1(4), s1(5), s1(6)]; % Velocity of meo debri at the start of the burn

meo_v1_mag = norm(meo_v1);
disp('The velocity of MEO s/c in km/s at the start of the burn is:  ')
disp(meo_v1_mag)

geo_r2 = [s2(1), s2(2), s2(3)]; % Position of GEO debri at the time when MEO finishes 5 orbits
geo_v2 = [s2(4), s2(5), s2(6)]; % velocity of GEO debri at the rendezvous

geo_v2_mag = norm(geo_v2);
disp('The velocity of GEO in km/s at the rendezvous is:  ')
disp(geo_v2_mag)

del_t = 1.258729719024421e+04; % in seconds, value used for now
string = 1; % prograde, CCW

% v_1: velocity we need to travel at initially to arrive at r2 wrt Lamberts
% v_2: velocity at r2 wrt Lamberts

[v_1, v_2] = lambert_ASp(meo_r1, geo_r2, del_t, string);
% v1 = norm(v_1)
% v2 = norm(v_2)
delta_v1 = norm(v_1 - meo_v1)
delta_v2 = norm(v_2 - geo_v2)
total_delV = delta_v1 + delta_v2;
disp('Total delta-v required for a Lamberts burn from MEO to GEO in km/s is:')
disp(total_delV)

% geo wait time: at the end of phase 3, amount of time our s/c stays with meo debri befor it
% reaches the optimal position to burn to reach GEO debri 

% Mission Elapsed Time after WAITING PERIOD: Total mission time right
% before the burn





