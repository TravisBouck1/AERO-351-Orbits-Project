function [Qgeo2peri] = Q_geo2peri_313_rotation_matrix(RAAN,inc,w)

% This function produces a transformation matrix, Q, from the PERIFOCAL to
% the GEOCENTRIC frames.

% Self, Justin
% Fall 2022

Cz_RAAN =        [  cos(RAAN)     sin(RAAN)     0;
                   -sin(RAAN)     cos(RAAN)     0;
                    0             0             1]; % 3

Cx_inc =        [   1             0             0;
                    0             cos(inc)      sin(inc);
                    0             -sin(inc)     cos(inc)]; % 1

Cz_w =        [     cos(w)        sin(w)        0;
                    -sin(w)       cos(w)        0;
                    0             0             1]; % 3

Qgeo2peri = Cz_w * Cx_inc *  Cz_RAAN;

end