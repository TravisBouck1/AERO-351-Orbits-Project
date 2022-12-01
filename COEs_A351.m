% Self, Justin 
% Orbits AERO 351, Fall 2022
% Common Orbital Elements Code; in class w Dr. A

close all; clear all; clc; 

%% COEs
% What planet?
muearth = 398600; % km3/s2
% Given state vectors
r0 = [9031.5 -5316.9 -1647.2]; % km
v0 = [-2.8640 5.1112 -5.0805]; % km/s

% create state vector for function call
state = [r0 v0];

% Call function
[r,v,h,inc,Omega,e,w,theta,epsilon,a,u,tau] = COEs_Justin(state,muearth);