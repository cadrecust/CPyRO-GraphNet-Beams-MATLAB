% Optimal_Design_SSFSRecBeams_Complete_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing volumes)
%    a beam element for all its three critical cross-sctions (left,middle
%    right)
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%----------------------------------------------------------------

clc
clear all

%% Geometry 
span=5000; % mm
b=300; % width (mm)
h=600; % height (mm)
        
hrec=50; % 
brec=50; % lateral concrete cover

%% Materials
fcu=35; % N/mm2

fy=500; % Yield stress of steel reinforcement (N/mm2)
wac=7.85e-6; % unit volume weight of the reinforcing steel (Kg/mm3)

%% Numerical model for analysis
dL=100; % mm
Ec=(3.46*sqrt(fcu)+3.21)*1e3;

%% Structural analysis FEM
% Supports
supportsLoc=[0 span]; % location of supports
wrange=[0, span];

Ac=b*h;
Ic=b*h^3/12;

%% Loads

W=[25,30]; % Uniformly Distributed Load over the whole beam  N / mm
nspans=length(W(:,1));

[R,U,V,M]=MSFSFEMBeams(span,Ac,Ic,Ec,supportsLoc,W,dL,wrange,0);

ne=zeros(1,nspans);
neSum=0;
load_conditions=[];
for i=1:nspans
    
    spanLengths(i)=supportsLoc(i+1)-supportsLoc(i);
    ne(i)=(supportsLoc(i+1)-supportsLoc(i))/dL;
    i1=neSum+1;
    Mleft=M(1,i1);
    
    neSum=neSum+ne(i);
    Mright=M(1,neSum);
    
    [Mmid,mp]=max(M(1,i1:neSum));
    load_conditions=[load_conditions;
                    i Mleft Mmid Mright]; %Kg-cm (flexure)
    
    %% Cut location ( local coordinates)
    cutxLocBeam(i,:)=cutLocationSSRecBeam(M(:,i1:neSum),dL);
end
load_conditions

%% Rebar data
% Available commercial rebar diameters (in eight-of-an-inch)
                %type diam
rebarAvailable=[1 6;
                2 8;
                3 10;
                4 12;
                5 16;
                6 20;
                7 25;
                8 32;
                9 40]; % mm^2

dvs=10;
pmin=0.003;
pmax=0.025;

%% OPTIMAL DESIGN 

hagg=20;

WUB=[0.6,0.7];
WND=[0.6];
WNB=[0.4];
WNC=[0.6,0.7];
Wcs1=2;
Wcs2=0.5;
WfaConstr=[WUB,WND,WNB,WNC,Wcs1,Wcs2];

nfig=1;
[volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db9Spans,EffSpans,...
MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
CFAspans]=OptimMSFSBeamsRebar3DSec(b,h,spanLengths,brec,hrec,hagg,...
pmin,pmax,rebarAvailable,fcu,load_conditions,fy,wac,cutxLocBeam,WfaConstr,1,1,...
[nfig],60,40);

for i=1:nspans
    % Average percentage of cross-section reinforcement 
    rho=sum(tenbLMRspan(i,:).*(db9Spans(i,:).^2*pi/4))./(b*h)/3;  
    
    [s1(i,1),s2(i,1),s3(i,1),d1(i,1),d2(i,1)]=shearDesignBeams(span,b,h,...
                                                    hrec,fcu,fy,V,dvs,rho);
    
    dvsBeams=dvs(i,1);
end

%% Side rebars (if necessary)
if h>=750
    [dSb,nSb,sepSb,distrSideBars]=sideBarsRecBeams3SecSpan(b,h,fy,...
            brec,hrec,tenbLMRspan,db9Spans,dvs,hagg,rebarAvailable);
    
    beamNSb(1,1)=2*nSb;
    diamlistdSb=zeros(2*nSb,1)+dSb;
    
    plotBeamSideBar3sec(b,h,-DistrRebarLeft,ListRebarDiamLeft,...
            DistrRebarMid,ListRebarDiamMid,-DistrRebarRight,...
            ListRebarDiamRight,diamlistdSb,distrSideBars,nfig);
else
    distrSideBars=[];
    diamlistdSb=[];
    beamNSb=[];
end

directionData='C:\Users\lfver\OneDrive\Desktop\OneDrive\CALDRECUST\Software\Package\Visual_CALDRECUST\Design_Data\RebarBeams\';

ExportDesignSSRecBeam(directionData,[b,h,span,brec,hrec],db9Spans,LenRebarL,...
    LenRebarM,LenRebarR,-DistrRebarLeft,DistrRebarMid,-DistrRebarRight,totnbSpan,...
    tenbLMRspan,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
    diamlistdSb,distrSideBars,beamNSb,[s1,s2,s3,d1,d2,dvsBeams]);
