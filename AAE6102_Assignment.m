clc;
clear all;
tic;
% Define working directory
addpath('Data','Navigation_utilities')
% Load data
load('rcvr.dat')
load('eph.dat')
%% Define constants
parameters = 4;
c = 299792458.0;                % speed of light in (m/s)
Wedot = 7.2921151467*10^(-5);     % WGS 84 value of earth's rotation rate (r/s)
mu = 3.986005*10^14;            % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
F = -4.442807633*10^(-10);        %Relativistic correction term constant

%% Calling columns from the raw ranging information (rcvr.dat)

rcvr_tow = rcvr(:,1);       % receiver time of week(s)
svid = rcvr(:,2);           % satellite PRN number(1-32)
pr = rcvr(:,3);             % pseudorange(m)
cycles = rcvr(:,4);         % number of accumulated cycles
phase = rcvr(:,5);          % to conver to (0 - 359.99), multiply by 360/2048
slp_dtct = rcvr(:,6);       % no. of cycle slip detected; non 0 = cycle slip
snr_dbhz = rcvr(:,7);       % signal to noise ration

% Calling columns from the ephemeris data from a GPS receiver (eph.dat)
tow = eph(:,1);         % receiver time of the week(s)
svid_e = eph(:,2);      % satellite PRN number (1 - 32)
toc = eph(:,3);         % reference time of clock parameters (s)
toe = eph(:,4);         % reference time of ephemiris parameters (s)
af0 = eph(:,5);         % clock correction coefficient -- group delay (s)
af1 = eph(:,6);         % clock correction coefficient (s/s)
af2 = eph(:,7);         % clock corrections coefficient (s/s/s)
ura = eph(:,8);         % user range accuracy (m)
e = eph(:,9);           % eccentricity (-)
sqrta = eph(:,10);      % square root of semi-major axis a (m**1/2)
dn = eph(:,11);         % mean motion correction (r/s)
m0 = eph(:,12);         % mean anomaly at reference time (r)
w = eph(:,13);          % argument of perigee (r)
omg0 = eph(:,14);       % right ascension (r)
i0 = eph(:,15);         % inclination angle at reference time (r)
odot = eph(:,16);       % rate of right ascension (r/s)
idot = eph(:,17);       % rate of inclination angle (r/s)
cus = eph(:,18);        % argument of latitude correction, sine (r)
cuc = eph(:,19);        % argument of latitude correction, cosine (r)
cis = eph(:,20);        % inclination correction, sine (r)
cic = eph(:,21);        % inclination correction, cosine (r)
crs = eph(:,22);        % radius correction, sine (r)
crc = eph(:,23);        % radius correction, cosine (r)
iod = eph(:,24);        % issue of data number 

%% Determining satellite positions and clocks

sv = length(svid);
satellite_x = zeros(8,7);
satellite_t = zeros(8,1);
for i = 1:sv
    id = svid(i);
    point = find(svid_e==id);
        [satellite_x(i,1:3), satellite_t(i,1)] = satellite_position((tow(point,1)),pr(point,1),toc(point,1),af0(point,1),af1(point,1),af2(point,1),...
            toe(point,1), m0(point,1),sqrta(point,1)^2,dn(point,1),e(point,1),w(point,1), omg0(point,1),...
            odot(point,1),i0(point,1), idot(point,1),cus(point,1),cuc(point,1),cis(point,1),cic(point,1),...
            crs(point,1),crc(point,1),toe(point,1));
end

%% Initial user coordinate position and clock bias

bias0 = 0;  % initial bias
position = [-2694685.473; -4293642.366; 3857878.924];   % initial position
[lat, lon, alt] = Wgsxyz2lla(position); % conversion of WGS XYZ coordinates to geodetic coordinates
n = length(pr);

while  true
    num = 1;
    for satellite_num = 1:n % determining the elevation, bearings and distance
        [bearing(satellite_num,1),elevation(satellite_num,1),distance(satellite_num,1)]=DistBearElev(position(1),position(2),position(3),...
            lat*pi/180,lon*pi/180,satellite_x(satellite_num,1),satellite_x(satellite_num,2),satellite_x(satellite_num,3));
% Defining the Design Matrix
        M(num,1) = (position(1) - satellite_x(satellite_num,1))/distance(satellite_num,1); % unknown X coordinate parameter
        M(num,2) = (position(2) - satellite_x(satellite_num,2))/distance(satellite_num,1); % unknown Y coordinate parameter
        M(num,3) = (position(3) - satellite_x(satellite_num,3))/distance(satellite_num,1); % unknown Z coordinate parameter
        M(num,4) =  1;  % receiver clock offset parameter 
        
        R(num,1) = distance(satellite_num,1) + bias0 - c*satellite_t(satellite_num,1); % Actual range
        R_observed(num,1) = pr(satellite_num,1) ;          % Observed range
        
        % COV Matrix
        Q(num,num) = 1 ;
        num=num+1;
        x=1;
    end
    % Misclosure
    Err = R_observed - R; % Error (which is the difference)
    % Adjustment and Computation
    W=diag((diag(Q).^(-1))); % W is the wieght of the matrix
    Normal = (M'*(W)*M); % Normal Matrix
    if (n >= parameters) 
        x   = (Normal^-1)*M'*(W)*(R_observed - R); % Solution
        residuals = (M*x + R) - R_observed; % Error estimation
    end
    position   = position   + x(1:3);
    bias0 = bias0 + x(4);
    delta_x = norm(x(1:3));
    if (delta_x < 1e-4)
        break;
    end
    
end
toc;


