% TLEs to epoch R and V vectors
% Self, Justin
% California Polytechnic State University, SLO
% Fall 2022

function [r_epoch,v_epoch,T,a,TA] = rvepoch_fromTLEs(UTC,ecc,inc,peri,apo,RAAN,w,Me)
%{
This function takes information from Two Line Elements and outputs the
associated position (r) and velocity (v) vectors for the given epoch.

INPUTS
    UTC = date and time in UTC in [yyyy mm dd hh mm ss] format.
    ecc = eccentricity
    Me = mean anomaly in degrees (straight from TLEs; will convert later)
    rp = perigee ALTITUDE, km (straight from TLEs; will convert to radius)
    ra = apogee ALTITUDE, km (straight from TLEs; will convert to radius)
    RAAN = in degrees; will convert to rad later
    inc = in degrees; will convert to rad later
    w = in degrees; will convert to rad later

OUTPUTS
    [r] = 3x1 position vector
    [v] = 3x1 velocity vector
    T = period, seconds
    a = semi-major axis, km
    TA = true anomaly, radians.

AT TIME OF EPOCH. 

%}

% Constants
rearth = 6378; % km
mu = 398600; % km2/s

% Convert TLEs to radians, radii, etc.
UT = [UTC(4:6)]; % time only
inc = deg2rad(inc);
rp = peri + rearth;
ra = apo + rearth;
w = deg2rad(w);
RAAN = deg2rad(RAAN);
Me = deg2rad(Me);

% Additional elements by hand
a = 0.5*(rp + ra); % km
T = ( (2*pi) / (sqrt(mu)) ) * a^(3/2); % sec

% convert to Julian date
[juliandate_epoch,UT_decimal,j0] = juliandate_justin(UTC(1),UTC(2),UTC(3),UT);

% Find true anomaly iteratively since TA not given in TLEs. 
TOL = 1e-8;
[E] = newtonsKepler(ecc,Me,TOL);

TA = 2* atan(  (tan(E / 2)) / (   sqrt(  (1 - ecc) / (1 + ecc)     )    )       ); % in radians now

% Determine epoch location (r,v) and orbital parameters.
h = findh(rp, mu, ecc, 0); % since radius = rp at TA = 0; use 0

% From TLE into perifocal frame
rperifocal = r_vect_perifocal(mu, h,TA,ecc);
vperifocal = v_vect_perifocal(mu, h,TA,ecc);

% Coordinate transform into ECI
[Qperi2geo] = Q_peri2geo_313_rotation_matrix(RAAN,inc,w);

% Q transpose * r_peri = r_eci
r_epoch = Qperi2geo' * rperifocal;
v_epoch = Qperi2geo' * vperifocal;

% convert to col vects
r_epoch = r_epoch';
v_epoch = v_epoch';

test_rad = norm(r_epoch);
test_alt = test_rad - rearth;

end

