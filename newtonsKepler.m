function E = newtonsKepler(ecc, Me, TOL)

% Self, Justin
% Aero351: Orbital Mechanics I

% This function uses Newton's method to solve Kepler's equation: 
% " Me = E - ecc*sin(E) " for E. 

% Inputs:
% ecc = eccentricity of the orbit
% Me = mean anomaly, Me = nt
% TOL = tolerance desired (error); use 1e-8

% Output:
% E = eccentric anomaly

% ---------------------------------------------------
% Figure out the best initial guess

if Me < pi
    Eguess = Me + .5*ecc;
else
    Eguess = Me - 0.5*ecc;
end

% Define f and fprime functions
f = @(E) Me - E + ecc*sin(E);
fprime = @(E) -1 + ecc*cos(E);

% Newton's Method
x0 = Eguess;
x1 = x0 - f(x0) / fprime(x0);
x = [x0; x1];
err = abs(x1 - x0);

while err > TOL
    x0 = x1;
    x1 = x0 - f(x0) / fprime(x0);
    err = abs(x1 - x0);
    x = [x; x1];
end

E = x(end);
end % function