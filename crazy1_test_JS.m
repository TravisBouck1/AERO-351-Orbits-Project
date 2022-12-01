% "Crazy" object in MEO/LEO
% Self, Justin
clear all; close all; clc; 

% This script is to test the new r,v vector retrieval function from input
% TLEs. 

% Epoch (UTC):	20 November 2022 15:00:57
% Eccentricity:	0.3883671
% inclination:	1.0827°
% perigee height:	11647 km
% apogee height:	34538 km
% right ascension of ascending node:	346.4374°
% argument of perigee:	293.5371°
% revolutions per day:	1.71597106
% mean anomaly at epoch:	202.1083°

% The function takes inputs in DEGREES and ALTITUDES.
% The function takes UTC as: 
%   UTC = date and time in UTC in [yyyy mm dd hh mm ss] format.

%% Prepare for function call
% define function inputs from TLEs

UTC = [2022 11 20 15 00 57];
ecc = 0.3883671;
inc = 1.0827;
peri = 11647;
apo = 34538;
RAAN = 346.4374;
w = 293.5371;
Me = 202.1083;

% Call function
[r_epoch,v_epoch] = rvepoch_fromTLEs(UTC,ecc,inc,peri,apo,RAAN,w,Me);

% Plot current location at TLE epoch
figure
   h1 = gca;
   earth_sphere(h1)
   hold on
% entire orbit
   plot3(r_epoch(1),r_epoch(2),r_epoch(3),'*','LineWidth',5)
  
% annotation
    th = text(r_epoch(1),r_epoch(2),r_epoch(3),'Current Location');
    set(th,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','right')

