function [r,v,h,inc,Omega,ecc,w,theta,epsilon,a,u,tau] = COEs_Justin(state,mu)

% Description: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:                          
% [state] = r0 vector and v0 vector
% mu = mu of massive body about which the spacecraft is orbiting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:                                                  
% r = position at state [r]; norm of the r_vector
% v = speed at state [v]; norm of the v_vector
% h = specific angular momentum, km2/s
% inc = inclination of orbit, report in degrees
% Omega = RAAN = right ascension of ascending node
% ecc = eccentricity
% w = omega = 
% theta = 
% epsilon = 
% a = 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    % Function begins below this line.
% ---------------------------------------------------------------------% 


% Norm the state vector(s)
r_vect = state(1:3);
v_vect = state(4:6);

r = norm(r_vect);
v = norm(v_vect);

disp("Position relative to center of massive body is: " + r + " km")
disp("Speed is: " + v + " km/s")

% radial velocity
veloc_rad = dot(r_vect,v_vect) / r;
disp("Radial velocity is: " + veloc_rad + " km/s")

% Direction check using rad_veloc
    if veloc_rad > 0 
    disp("Flying away from periapse")
    else 
    disp("Flying toward periapse")
    end

%% h
h_vect = cross(r_vect,v_vect);
h = norm(h_vect);
disp("h is: " + h + " km^2 / s")

%% Inclination, report in degrees
hz = h_vect(3);
inc = acosd(hz/h);

% Quadrant check for inc
if (inc > 90) && (inc < 180)
    disp("Object is in retrograde")
else
    disp("Object is in prograde")
end
disp("Inclination (inc) is: " + inc + " degrees")

% Node line, N
K = [0 0 1];
N_vect = cross(K,h_vect);
N = norm(N_vect);

%% RAAN, capital OMEGA
Omega = acosd(N_vect(1)/N);

disp("Omega (RAAN) is: " + Omega + " degrees")

% Quadrant check RAAN
    if (N_vect(1) / N) > 0
        disp("RAAN lies in QI or QIV")
        if N_vect(2) > 0
            disp("NAAN lies in QI")
        else
        end
    else 
        disp("RAAN lies in QII, QIII")
        if N_vect(2) > 0 
            disp("RAAN is in QII")
        else
            Omega = 2*pi - acos(N_vect(1) / N);
            Omega = rad2deg(Omega);
        end
    end

%% Eccentricity
ecc_vect = (mu^-1)* (cross(v_vect,h_vect) - mu* (r_vect / r)  );
ecc = norm(ecc_vect);

% Ecc heart check
    if ecc == 0
        disp("Orbit is circular")
    elseif (ecc > 0) && (ecc < 1)
        disp("Orbit is elliptical")
    elseif ecc == 1
        disp("Orbit is parabolic")
    elseif ecc > 1
        disp("Orbit is hyperbolic")
    end

%% Argument of perigee, omega (little w)
w = acos( dot(N_vect,ecc_vect) / (N*ecc) );
w = rad2deg(w);

% Quadrant check for w
 
    if ecc_vect(3) < 0
        w = 2*pi - acos(dot(N_vect,ecc_vect) / (N*ecc) );
        w = rad2deg(w);
    else
        if dot(N_vect,ecc_vect) > 0
            disp("w is in QI, QIV")
        else 
            disp("w is in QII, QIII")
        end % inner
    end % out
disp("Argument of periapse, w, is: " + w + " degrees")

%% True Anomaly, theta

% Check first.
if veloc_rad >= 0
    theta = acos( dot(ecc_vect, r_vect) / (ecc*r) );
    theta = rad2deg(theta);
else 
    theta = 2*pi - acos( dot(ecc_vect, r_vect) / (ecc*r) );
    theta = rad2deg(theta);
end 

disp("True anomaly, theta, is: " + theta + " degrees")

%% Other important orbital elements

% Specific energy, epsilon
epsilon = ((v^2)/2)  - (mu/r);

% specific energy check
if epsilon < 0 
    disp("Specific energy is negative; orbit is bounded. DOES THIS MATCH with Ecc?")
else
    disp("Specific energy is nonnegative; orbit is unbounded. DOES THIS MATCH with Ecc?")
end


% semi- major axis, a
a = -(mu/(2*epsilon));

% Argument of lattitude, u
u = w + theta;

%% Ecc anomaly, Me, and time (tau)
% Check orbit to determine method
T = 2*pi* (sqrt(a^3 / mu)); % period
Thours = T/3600; % period; sec to hours

disp("Period is: " + Thours + " hours")

n = (2*pi) / T;

if ecc == 0
    % circ orbit method
    t = ((h^3) / (mu^2))* theta;
elseif (ecc>0) && (ecc<1)
    % elliptical orbit method
    E = 2 * atan(sqrt((1- ecc)/(1+ecc))  * tan(theta/2) );
    Me = E - ecc*sin(E);
    t = Me / n; % n is mean motion.
end 

% Time of periapsis passage, tau
tau = t/3600; % t in hours; 
tau = -tau; % since it is time SINCE, we naturally have a (-). Let's fix that.

disp("Time since periapse passage, t, is: " + tau + " hours")

end % function
