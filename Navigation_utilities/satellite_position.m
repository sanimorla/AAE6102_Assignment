function [rs, dts] = satellite_position(tow,pr,toc,af0,af1,af2, toe, m0,A,dn,e,w,omg0,odot,i0,idot,cus,cuc,cis,cic,crs,crc,toes)
%The function satellite_position calculates the positions of satellites at
%time (tow) using the ephemeris parameters and the satellite clock bias.
% The Input Args include the receiver time, and ephemeris parameters as
% defined below;
    %tow = the receiver time of the week
    % pr = pseudorange
    %svid_e = satellite PRN number
    %toc  = reference time of clock parameters 
    %toe =  reference time of ephemiris parameters 
    %af0 = clock correction coefficient -- group delay 
    %af1 =  clock correction coefficient
    %af2 = clock corrections coefficient 
    %ura = user range accuracy
    %e = eccentricity (-)
    %sqrta = square root of semi-major axis a 
    %dn = mean motion correction 
    %m0 = mean anomaly at reference time
    %w = argument of perigee 
    %omg0 = right ascension 
    %i0 = inclination angle at reference time 
    %odot = rate of right ascension 
    %idot = rate of inclination angle 
    %cus = argument of latitude correction, sine
    %cuc = argument of latitude correction, cosine
    %cis = inclination correction, sine 
    %cic = inclination correction, cosine 
    %crs = radius correction, sine 
    %crc = radius correction, cosine
    %iod = issue of data number 
% The Output Args [rs, dts] gives the X,Y,Z of the satellite and Clock bias
% repectively.

% Reference: Code modified from Mohammed Abougalala(2021)

C = 299792458.d0;
% transmission time by satellite clock 
tow=tow-(pr/C);
% satellite clock bias by broadcast ephemeris 
dt = eph2clk(tow, toc, af0, af1, af2);
tow=tow-dt;

% satellite position and clock at transmission time 
[rs, dts] = eph2pos(tow,toe,m0,A,dn,e,w,omg0,odot,i0,idot,cus,cuc,cis,cic,crs,crc,toes,af0,af1,af2);
end 

% !----------------------------------------------------------------------------
% ! broadcast ephemeris to satellite clock bias -------------------------------
% ! compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
% ! args   : gtime_t time     I   time by satellite clock (gpst)
% !          eph_t *eph       I   broadcast ephemeris
% ! return : satellite clock bias (s) without relativeity correction
% ! notes  : see ref (1),(7),(8)
% !          satellite clock does not include relativity correction and tdg
% !----------------------------------------------------------------------------
function eph2clk = eph2clk(tow, toc, af0, af1, af2)
% implicit none
% % type(gtime_t),intent(in) :: time
% type(eph_t), intent(in) :: eph
% real*8 t
% integer*4 i
t=tow-toc;
for i=1:2
    t=t-(af0+af1*t+af2*t*t);
end
eph2clk=af0+af1*t+af2*t*t;
end

%--------------------------------------------------------------------------
% ! broadcast ephemeris to satellite position and clock bias --------------------
% ! compute satellite position and clock bias with broadcast ephemeris (gps,
% ! galileo, qzss)
% ! args   : gtime_t time     I   time (gpst)
% !          eph_t *eph       I   broadcast ephemeris
% !          real*8 *rs       O   satellite position (ecef) {x,y,z} (m)
% !          real*8 *dts      O   satellite clock bias (s)
% !          real*8 *var      O   satellite position and clock variance (m^2)
% ! return : none
% ! notes  : see ref (1),(7),(8)
% !          satellite clock includes relativity correction without code bias
% !          (tgd or bgd)
% !-----------------------------------------------------------------------------
function [rs, dts] = eph2pos(tow,toe,m0,A,dn,e,w,omg0,odot,i0,idot,cus,cuc,cis,cic,crs,crc,toes,af0,af1,af2)

C = 299792458.d0; % speed of light
MU_GPS = 3.9860050d14;          % gravitational constant
OMGE   = 7.292115d-5;           % earth angular velocity (rad/s) ref (2) 
RTOL_KEPLER = 1d-13;            % relative tolerance for Kepler equation 
MAX_ITER_KEPLER = 30;     % max number of iteration of Kelpler 

if(A<=0.d0)
    rs=0.d0; dts=0.d0;
    return;
end

tk=(tow-toe);

mu=MU_GPS; omge1=OMGE;

M=m0+(sqrt(mu/(A*A*A))+dn)*tk;
n=0;E=M;Ek=0.d0;

while(abs(E-Ek)>RTOL_KEPLER && n<MAX_ITER_KEPLER)
    Ek=E; E=E-(E-e*sin(E)-M)/(1.d0-e*cos(E));
    n=n+1;
end

if(n>=MAX_ITER_KEPLER), return; end
sinE=sin(E); cosE=cos(E);

u=atan2(sqrt(1.d0-e*e)*sinE,cosE-e)+w;
r=A*(1.d0-e*cosE);
i=i0+idot*tk;
sin2u=sin(2.d0*u); cos2u=cos(2.d0*u);
u=u+cus*sin2u+cuc*cos2u;
r=r+crs*sin2u+crc*cos2u;
i=i+cis*sin2u+cic*cos2u;
x=r*cos(u); y=r*sin(u); cosi=cos(i);


O=omg0+(odot-omge1)*tk-omge1*toes;
sinO=sin(O); cosO=cos(O);
rs(1)=x*cosO-y*cosi*sinO;
rs(2)=x*sinO+y*cosi*cosO;
rs(3)=y*sin(i);

tk=(tow-toe);
dts=af0+af1*tk+af2*tk*tk;

% relativity correction 
dts=dts-2.d0*sqrt(mu*A)*e*sinE/(C*C);

% position and clock error variance 
% var=var_uraeph(sys,eph%sva)
end 

%--------------------------------------------------------------------------
