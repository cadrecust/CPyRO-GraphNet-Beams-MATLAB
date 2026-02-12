function [s1,s2,s3,d1,d2]=shearDesignBeams(span,b,h,hrec,fcu,fy,...
                          shearbeam,dvs,rho)
                      
%------------------------------------------------------------------------
% Syntax:
% [s1,s2,s3,d1,d2]=shearDesignBeams(span,b,h,rec,fcu,fy,shearbeam,diams,rho)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: N,mm
%
%------------------------------------------------------------------------
% PURPOSE: To design the separation of the transversal reinforcement along
% the whole length of a beam element according to the mechanic shear
% forces distribution.
% 
% OUTPUT: s1:                   is the separation of the transversal 
%                               reinforcement at the left part/section of
%                               the beam element.
%
%         s2:                   is the separation of the transversal 
%                               reinforcement at the middle part/section 
%                               of the beam element.
%
%         s3:                   is the separation of the transversal 
%                               reinforcement at the right part/section 
%                               of the beam element.
%                               
%         d1:                   is the length along which the transversal
%                               reinforcement is separated by s1 (cm).
%                               Left part length.
%                               
%         d2:                   is the length along which the transversal
%                               reinforcement is separated by s3 (cm).
%                               Right part length.
%
% INPUT:  span:                 is the length of the beam element
%
%         b,h:                  are the cross-section dimensions: width 
%                               and height, respectively.
%
%         hrec:                 is the concrete cover
%
%         shear_beam:           is the array containing the shear forces
%                               distribution from left to right.
%
%         fcu, fy:              is the compressive strength of the concrete
%                               and the yield stress of the reinforcing
%                               steel.
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%------------------------------------------------------------------------

dtable=[125 160 175 200 225 250 300 400];

vcdata=[0.15 0.45 0.43 0.41 0.40 0.39 0.38 0.36 0.34;
        0.25 0.53 0.51 0.49 0.47 0.46 0.45 0.43 0.40;
        0.50 0.67 0.64 0.62 0.60 0.58 0.56 0.54 0.50;
        0.75 0.77 0.73 0.71 0.68 0.66 0.65 0.62 0.57;
        1.00 0.84 0.81 0.78 0.75 0.73 0.71 0.68 0.63;
        1.50 0.97 0.92 0.89 0.86 0.83 0.81 0.78 0.72;
        2.00 1.06 1.02 0.98 0.95 0.92 0.89 0.86 0.80;
        3.00 1.22 1.16 1.12 1.08 1.05 1.02 0.98 0.91];
    
np=length(shearbeam(1,:));

%% Cross-section shear resistance of the beam (Vcr)

if fcu<=40
    vr=0.4;
elseif 40<fcu && fcu<=80
    vr=0.4*(fcu/40)^(2/3);
end

%% Stirrups separation
d=h-hrec;

maxsv=0.75*d;

ab=dvs^2*pi*0.25;

k1=100*rho;
i=1;
while k1>vcdata(i,1) && k1>vcdata(i+1,1)
    
    if vcdata(i,1)<=k1 && vcdata(i+1)>=k1
        break;
    end
    if i==length(vcdata(:,1))-1
        break;
    end
    i=i+1;
end

j=1;
while d>dtable(1,j) && d>dtable(1,j+1)
    
    if d>=dtable(1,j) && d<=dtable(1,j+1)
        break;
    end
    if j==length(dtable)-1
        break;
    end
    j=j+1;
end

vc=vcdata(i,j+1);

vmax=max(abs(shearbeam(1,:)))/(b*d);
if vmax>(vc+vr)
    disp('shear reinforcement is needed');
end

%% Left part
v1=abs(shearbeam(1,1))/(b*d);
if 0.5*vc<v1 && v1<=(vc+vr) 
    s1=2*ab*0.87*fy/(vr*b);
elseif (vc+vr)<v1 && v1<min([0.8*sqrt(fcu),7])
    s1=2*ab*0.87*fy/((v1-vc)*b);
elseif vmax<=0.5*vc
    s1=2*ab*0.87*fy/(vr*b);
elseif 0.5*vc>=v1
    s1=maxsv; % mm
end

s1=s1-mod(s1,5);

if s1<60
    s1=60;
elseif s1>maxsv
    s1=maxsv;
end

%% Right part
v3=abs(shearbeam(1,np))/(b*d);

if 0.5*vc<v3 && v3<=(vc+vr)
    s3=2*ab*0.87*fy/(vr*b);
elseif (vc+vr)<v3 && v3<min([0.8*sqrt(fcu),7])
    s3=2*ab*0.87*fy/((v3-vc)*b);
elseif vmax<=0.5*vc
    s3=2*ab*0.87*fy/(vr*b);
else
    s3=maxsv;
end

s3=s3-mod(s3,5);

if s3<60
    s3=60;
elseif s3>maxsv
    s3=maxsv;
end

%% Middle part
s2=2*ab*0.87*fy/(0.4*b);

s2=s2-mod(s2,5);
vr2=(2*ab*0.87*fy*d)/s2;

if s3>s2
    s3=s2;
end
Vus2=(vc*b*d+vr2);

% To determine the position in which the separation of the stirrups 
% will change from s1 to s2 along the beam beam length, starting from
% the left

ps21=[];
for i=1:np-1
    if abs(shearbeam(1,i))>=Vus2 && abs(shearbeam(1,i+1))<=Vus2
        
        ps21=i+1;
        break;
    end
end

ps22=[];
for i=1:np-1
    if abs(shearbeam(1,np+1-i))>=Vus2 && abs(shearbeam(1,np-i))<=Vus2
        ps22=np-i;
        break;
    end
end

lenbeam=linspace(0,span,np);
if isempty(ps21)==0
    d1=lenbeam(ps21);
else
    d1=ceil(span/2);
end
if isempty(ps22)==0
    d2=span-lenbeam(ps22);
else
    d2=span-d1;
end