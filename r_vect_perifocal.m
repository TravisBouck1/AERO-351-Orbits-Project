function r = r_vect_perifocal(mu, h, theta,ecc)

% Self, Justin
% Fall 2022
% This function determines the position vector of a hyperbolic orbit

% INPUTS:
% mu: grav param
% h = specific angular momentum (scalar)
% theta = true anomaly, degrees

r = (h^2/mu) * (1 + ecc * cos(theta))^-1 *  [cos(theta); sin(theta); 0]; % p, q, w coordinate vects.

end