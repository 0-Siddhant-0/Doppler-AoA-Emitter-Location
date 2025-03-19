function [Xp,Yp,Zp,Vx,Vy,Vz,long_vect]=weave(g,T,del_T,alt_kft,vel)
%
% USAGE: [Xp,Yp,Zp,Vx,Vy,Vz,long_vect]=weave(g,T,del_T,alt_kft,vel);
%
% Inputs: g = acceleration of platform  (g>0 curves down; g<0 curves up)
%         T = Final Time in Seconds (Initial Time = 0)
%         del_T = measurement time spacing in seconds
%         alt_kft = platform altitude in thousands of feet
%         vel = platform velocity in m/s
%
%
% Outputs: Xp, Yp, Zp = XYZ positions of antenna #1 (in meters)
%          Vx, Vy, Vz = velocities of platform (in m/s)
%          long_vector = unit vector pointing along longitudinal baseline
%          trans_vector = unit vector pointing along transversal baseline



tend=T;


dt = del_T; dtt=del_T/20;
%%% vel = 200; 
accelh = g*9.81; accelv=0.5*9.81;  alt0 = alt_kft*1000*0.3048;

dtr=pi/180;

t = 0:dt:tend ;   % vector of times of measurements
maxturn = 30*dtr ;
vhoriz = vel ;   % assume cos(gamma) close to 1 & vert axis can be decoupled
radturn = (vhoriz^2)/accelh ;   % turn radius
turnrate = vel/radturn;
turndur = 2*maxturn/turnrate ;   % turn duration
turndur = fix(turndur/dtt)*dtt ;
gamma0 = -asin(min(accelv*turndur/vel,0.5))/2 ; %accelv & gamma0 assumed small
dtt = 1/20 ;   % time increment for turn generation

vx =  vel*cos(gamma0)*cos(maxturn) ;   % initialize velocities
vy = -vel*cos(gamma0)*sin(maxturn) ;
vz =  vel*sin(gamma0) ;

XYZA(1,1) = 0 ;
px = 0 ;  py = 0 ;  pz = alt0 ;
itt=0 ;   % count of trajectory integration times
im=0  ;   % count of measurement times

n=0;
for tt = 0 : dtt : tend
n=n+1;
trem = rem(tt,2*turndur) ;
 if trem < turndur
  turnangl = maxturn - turnrate*trem ;
  turn_angle(n)=turnangl;
  ahoriz = accelh ;
  az = accelv ;  %Note: vert & horiz decoupled; only good for small az & gamma
 else
  turnangl = -maxturn + turnrate*(trem-turndur) ;
  turn_angle(n)=turnangl;
  ahoriz = -accelh ;
  az = -accelv ;
 end   % if trem
vx =  vhoriz*cos(turnangl) ;
vy = -vhoriz*sin(turnangl) ;
vz = vz + az*dtt ;
 if itt > 0
  px = px + vx*dtt ;
  py = py + vy*dtt ;
  pz = pz + vz*dtt ;
 end   % if itt

 if rem(1e-7 + itt*dtt,dt) <1e-6   % if remainder=0 with tolerance
  im = im+1 ;
  XYZADOT(1,im) = vx ;
  XYZADOT(2,im) = vy ;
  XYZADOT(3,im) = vz ;
  XYZA(1,im) = px ;
  XYZA(2,im) = py ;
  XYZA(3,im) = pz ;
  ah(im) = ahoriz ;
 end   % if rem
itt = itt + 1 ;
end   % for tt
NMEAS = im;

hdg = pi/2 - atan2(XYZADOT(2,:),XYZADOT(1,:)) ;   % (1,n)
roll = atan(-ah/9.81) ;   % approx. for small pitch
pitch = atan(XYZADOT(3,:)./sqrt(XYZADOT(1,:).^2 + XYZADOT(2,:).^2)) ;  % (AOA=0)
XYZhpr=[XYZA;hdg;pitch;roll];

Xp=XYZA(1,:);Yp=XYZA(2,:);Zp=XYZA(3,:);
Vx=XYZADOT(1,:);Vy=XYZADOT(2,:);Vz=XYZADOT(3,:);

long_vect=[sin(hdg);cos(hdg);sin(pitch)];
vect_norms=sqrt(sum(long_vect.^2));    %%% Compute norms of vectors
long_vect=long_vect./vect_norms(ones(1,3),:);  %%% Normalize the vectors


trans_vect=[cos(hdg);-sin(hdg);-sin(roll)];
vect_norms=sqrt(sum(trans_vect.^2));    %%% Compute norms of vectors
trans_vect=trans_vect./vect_norms(ones(1,3),:);  %%% Normalize the vectors
