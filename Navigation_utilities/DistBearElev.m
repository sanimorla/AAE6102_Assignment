function [bearing,elevation,distance]=DistBearElev(x,y,z,lat,lon,X,Y,Z)
% The function (DistBearElev) calculates (returns) the bearing/Azimuth, elevation and
% distance of the satellites.
% The INPUT ARGS include the initial receiver positions (x,y,z), the
% latitude and longitudes (lat, lon) in radians, and the satellite positions
% (X,Y,Z).
% OUTPUT ARGS include bearing, elevation, distance.
%
% Reference: Code modified from Mohammed Abougalala (2021)

rlat=lat-pi/2.0;
rlon=lon-pi;
srlat=sin(rlat);
crlat=cos(rlat);
srlon=sin(rlon);
crlon=cos(rlon);

delta_x=X-x; % change in x
delta_y=Y-y; % change in y
delta_z=Z-z; % change in z

du=crlat*crlon*delta_x+ crlat*srlon*delta_y -srlat*delta_z;
dv=srlon*delta_x      - crlon*delta_y;
dw=srlat*crlon*delta_x+ srlat*srlon*delta_y+ crlat*delta_z;


bearing  = atan2(dv,du)/pi*180;
distance= sqrt(du^2+dv^2+dw^2);
elevation  = asin(dw/distance)/pi*180;


