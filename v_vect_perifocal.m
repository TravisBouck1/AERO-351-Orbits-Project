function v = v_vect_perifocal(mu, h,theta,ecc)

% Self, Justin
% Fall 2022
% This function determines the velocity vector of a hyperbolic orbit

% INPUTS:
% mu: grav param
% h = specific angular momentum (scalar)
% theta = true anomaly, degrees

v = (mu/h) * [-sin(theta); ecc + cos(theta); 0]; % p, q, w coordinate vects.

end