% AERO 351 ODE45 help
% by everyone who was honest (LOL) about not understanding ODE45
% This is posted on CANVAS

% Housekeeping
clear all; close all; clc; 

%% Define variables

r_vect = [3450 -1700 7750]; % position vector, km
velocity_vect = [5.4 -5.4 1]; % km/s
mu_earth = 398600; % mu of Earth in km3/s2
timespan = [0 24*3600]; % 24h in seconds

% set tolerances for ode45 (SUPER IMPORTANT FOR ORBITS)
options = odeset('RelTol',1e-8, 'AbsTol', 1e-8); % USE THROUGHOUT ORBITS 

% call ode function
[time_output, state_new] = ode45(@pig, timespan, [r_vect velocity_vect],options, mu_earth); % pig is our function

% Plot it

figure
plot(state_new(:,1),state_new(:,2),'Linewidth',2)
xlabel('x, km')
ylabel('y, km')


