% Find the true anomaly of an object that is coasting... not burning... on
% a given orbit after time period has passed.

function [ta_degrees,Me_degrees, E_degrees,newtime,newstate] = true_anomaly_from_orbit_times(tspan,state,mu,ecc,n,Me0)

%{
INPUTS: 
        tstart          = start time, seconds, from last event
        duration        = time in seconds of requested timespan
        mu              = grav param
        Me0             = Mean anomaly (Me) in DEGREES associated with t0. (from the end of last event) 
        ecc             = eccentricity of orbit
        postion_at_t0   = 3x1 r vector of object at t0
        
OUTPUTS:
        ta          = true anomaly at specified time, DEGREES
        Me          = mean anomaly at new time DEGREES

%}
  
%% Function
% call ode here
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
[newtime, newstate] = ode45(@non_impulsive_COAST,tspan,state,options,mu); % this should give us r, v at t = X hours (specified by user)

% Find true anomaly iteratively since TA not given in TLEs. 
TOL = 1e-8;

% Find new Mean Anomaly for current timestep to calc current E and TA. **
% whole point of this function!!!!
Me = deg2rad(Me0) + n* newtime(end); % Me is associated with TLEs

% find E
[E] = newtonsKepler(ecc,Me,TOL);

% Find TA
ta = 2* atan(  (tan(E / 2)) / ( sqrt(  (1 - ecc) / (1 + ecc) ) ) ); % in radians now

% Final outputs
Me_degrees = rad2deg(Me);
E_degrees = rad2deg(E);
ta_degrees = rad2deg(ta);
end