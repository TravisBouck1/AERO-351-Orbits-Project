%% Lambert's Function, Curits Algorithm 5.2
% Anmol Sharma
% This function uses known two position vectors on the path
% of mass around M and determines velocity at both points

function[v1, v2] = lambert_ASp(r1, r2, del_t, string)
r1vec = r1;
r2vec = r2;

% INPUT
% r1 & r1vec : initial position vector
% r2 & r2vec : final position vector
% del_t : time of flight from r1 to r2, a constant
% string: Prograde or Retrograde check

% OUTPUT
% v1: Initial velocity vector
% v2: Final velocity vector

mu = 398600; % Gravitational Parameter of Earth in km^3/s^2
% 1. Calculate r1 and r2 
r1_mag = norm(r1vec);
r2_mag = norm(r2vec);

% 2. Determine Prograde or Retrograde
% Setting conditions to check for -grade
cross_prod = cross(r1vec, r2vec);
theta = acos(dot(r1vec,r2vec)/(r1_mag*r2_mag)); % angle between r1 and r2

if string == 1 % Prograde
    if cross_prod(3) >= 0
        del_theta = theta;
    else 
        del_theta = 2*pi - theta;
    end
elseif string == 2 % Retrograde
    if cross_prod < 0 
        del_theta =  theta;
    else
        del_theta = 2*pi - theta;
    end
end

% 3. Calculate A using Eq 5.35; 'A' is constant to determine

A = sin(del_theta) * sqrt((r1_mag*r2_mag)/(1 - cos(del_theta)));

% 4. Solve for z using iteration
% Using Newton's method, we form the function
% sign change; z's sign tells us if the orbit is hyperbola (z<0)
% parabola (z=0) or ellipse (z=0)
%z = -100;

% Function of z given by equation 5.38
y = @(z) r1_mag + r2_mag + A*( (z*S(z) - 1)/sqrt(C(z)) );

% Eqn (5.37)
X = @(z) sqrt( y(z) / C(z) );

% Deltateqn (From lecture)
deltateqn = @(z) (( (X(z))^(3)) *S(z)) / (sqrt(mu)) + (A * sqrt( y(z)) ) / (sqrt(mu));

% Finally, we can use the bisection method to find z.
a = -4 * pi^2; % lower bound, from lecture
b = 4 * pi^2; % upper bound, from lecture
TOL = 1e-8; % tolerance on percision of convergence

z = (a + b) / 2;
while abs(deltateqn(z) - del_t) > TOL
    if deltateqn(z) <= del_t
        a = z; % reset to z lower
    else
        b = z; % reset to z upper
    end % if
    z = (a+b) / 2;
end % while

%{
while F(z,del_t) < 0
    z = z+0.1;
end

% error tolerance
tol = 1e-8;
nmax = 5000;
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
    n = n+1;
    ratio = F(z,del_t)/dFdz(z);
    z = z-ratio;
end

if n >= nmax
    disp('Max number of iterations exceeded.')
end
% 5. Its made into function 'y'
%}

% 6. Calculate Lagrange f,g and gdot
f = 1- (y(z)/r1_mag);
g = A*sqrt(y(z)/mu);
l = 1/g;
gdot = 1 - y(z)/r2_mag;
%fdot = (sqrt(mu) / (r1*r2)) * sqrt( y(z) / C(z)) * ( z * S(z) - 1 );

v1 = (1/g)*(r2vec - f*r1vec);
v2 = (1/g)*(gdot*r2vec - r1vec);

return











