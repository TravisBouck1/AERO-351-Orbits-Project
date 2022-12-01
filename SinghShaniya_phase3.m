clear all; close all; clc

global mu
deg = pi/180;
radiusEarth = 6378;

mu = 398600; % km3/s2 

r1 = [1.096035657606564e+04 2.698333735244205e+04 1.194164375722560e+02]; 
r2 = [-1.783230987211009e+04 1.214512502799535e+04 -4.006338669676189e+02]; 

dt = 7.113944104937772e+03; 

string = 'pro'; 

[v1, v2] = lambert(r1, r2, dt, string);

v_crazy = [-2.985109599699163 -1.186895083866727 0.014082710185299];
v_meo = [-3.772645841393155 -3.030358641398373 -0.042496438145642]; 


% Find Delta-V 

%subtracting the vectors then taking the norm
delta_v_crazy = v1 - v_crazy; 
delta_v_meo = v2- v_meo; 

v_crazy_norm = sqrt(delta_v_crazy(1)^2 + delta_v_crazy(2)^2 + delta_v_crazy(3)^2);
v_meo_norm = sqrt(delta_v_meo(1)^2 + delta_v_meo(2)^2 + delta_v_meo(3)^2);
deltaV_total = v_meo_norm + v_crazy_norm; 

fprintf('\n Total delta-V (km/s) = %.4f (km/s)', deltaV_total)
   
