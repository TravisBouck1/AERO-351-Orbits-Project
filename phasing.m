function [deltaV,Time,h2,T2,ecc2,rD] = phasing(rp,ra,mu,theta,divisions)
% Phasing Manuever with a given number of revoltions
% inputs:
% 1. rp: Radius of Periapse
% 2. ra: Radius of Apoapse
% 3. mu
% 4. theta: True Anomaly of Prey
% 6. Divisions (How many periods before rendevous)
% outputs:
% 1. deltaV
% 2. Time 


% TRUE ANOMALY OF CHASER ASSUMED TO BE PERIGEE

% Initial Orbital Elements
ecc1 = (ra - rp)/(ra+rp);
h1 = sqrt(2*mu)*sqrt((ra*rp)/(ra+rp));
a1 = (rp+ra)/2;
T1 = 2*pi*sqrt((a1^3)/mu); 

% Prey Ecentric Anomomaly
E = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tan(theta/2));

% Flight time from Chaser to Prey
tAB = (T1/(2*pi))*(E - ecc1*sin(E));

% Phasing maneuver period
T2 = T1 - tAB/divisions;

% semi major axis of orbit 2
a2 = ((sqrt(mu)*T2)/(2*pi))^(2/3);

% Apogee of orbit 2
rD = 2*a2 - rp;

% ecc of orbit 2
ecc2 = (rD-rp)/(rD+rp);

% h of orbit 2
h2 = sqrt(2*mu)*sqrt((rD*rp)/(rD + rp));

% inital velocity change
Va1 = h1/rp;

% final velocity change
Va2 = h2/rp;

% Delta V
deltaV = 2*(abs(abs(Va1) - abs(Va2)));

% Total Time
Time = T2*divisions;



end