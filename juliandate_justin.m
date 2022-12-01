function [juliandate_justin, UT_decimal, j0] = juliandate_justin(y,m,d,UT)
%   This function outputs the julian date (JD) with inputs year, month,
%   day, and GMT (UT)

% y = four digit year
% m = two digit month
% d = two digit day
% UT = array that contains the UT in the form [hour min sec]

% floor(x) rounds x to negative infinity; or rounds "down" so no decimals
j0 = 367*y - floor((7*(y+floor((m+9)/12)))/4) + floor((275*m)/9) + d + 1721013.5;
% from Curtis book and Dr. A lecture 9/19/2022

% j0 = Julian date at 0 hours UT (noon)
% Now ww need to convert HH:MM:SS to HH.HH (decimal hours)

UT_decimal = UT(1) + UT(2)/60 + UT(3)/3600; % 60 min / hour, 3600 sec / hour
juliandate_justin = j0 + (UT_decimal/24);

end