function [deltaV, time,Vp1,Vp2,deltaV2,h1,h2,h3] = Hohmann(p1,a1,p3,a3,mu)
% Calculates the delta V and time for a Hohmannn Transfer

% Inputs: 
% Initial radius periapse
% Initial radius periapse
% Final radius periapse
% Final radius periapse
% mu
% Outputs: deltaV, time


% h for initial(1), transfer(2), and final(3) orbits
h1 = sqrt(2*mu)*sqrt((p1*a1)/(p1 + a1));
h2 = sqrt(2*mu)*sqrt((p1*a3)/(p1 + a3));
h3 = sqrt(2*mu)*sqrt((p3*a3)/(p3 + a3));

% speed at periapse of orbit 1
Vp1 = h1/p1;

% speed at periapse of orbit 2
Vp2 = h2/p1;

% delta V of first burn 
deltaV1 = Vp2 - Vp1;

% speed at apoapse of orbit 2
Va2 = h2/a3;

% speed at apoapse of orbit 3
Va3 = h3/a3;

% delta V of second burn 
deltaV2 = Va2 - Va3;

% Total Delta V
deltaV = abs(deltaV1) + abs(deltaV2);

% Time of Transfer
a = 0.5*(p1+a3);
T = (2*pi*a^(3/2))/sqrt(mu);
time = T/2;


end