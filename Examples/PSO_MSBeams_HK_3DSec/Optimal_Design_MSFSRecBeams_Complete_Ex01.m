% Optimal_Design_MSFSRecBeams_Test_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing weight)
%    a two-span beam element.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2025-07-03
% Copyright (c)  School of Engineering
%                HKUST
%----------------------------------------------------------------

clc
clear all

%% Insert folder path to store rebar design results
directionData='C:\Users\lfver\OneDrive\Desktop\OneDrive\CALDRECUST\Software\Package\Visual_CALDRECUST\Design_Data\RebarMSBeams\';

% If results are to be saved, input saveResults=true
saveResults=true

%% Insert geometry of beam element
bv=[300;300]; % mm
hv=[600;600]; % mm

Lv=[4000;5000]; % mm
Totalspan = sum(Lv); % mm

%% Insert materials' properties
fcuv=[35;35]; % MPa

%% Loads
Wv=[20,20;
    25,25]; % N/mm
Wrange=[0, 4000;
        4000, 9000]; % mm

%% FEM model
dL=100;

supportsv=[0,4000,9000]; % mm
 
nspans=length(Wv(:,1));

Ac = bv .* hv;
Ic = bv .* hv .^ 3 / 12;
Ec = (3.46 * sqrt(fcuv) + 3.21) .* 1e3;

[R,U,V,M]=MSFSFEMBeams(Totalspan,Ac,Ic,Ec,supportsv,Wv,dL,Wrange,1);

ne=zeros(1,nspans);
neSum=0;
load_conditions=[];
cutxLocMSBeams=[];
for i=1:nspans
    
    ne(i)=(supportsv(i+1)-supportsv(i))/dL;
    i1=neSum+1;
    Mleft(i,1)=M(1,i1);
    
    neSum=neSum+ne(i);
    Mright(i,1)=M(1,neSum);
    
    [Mmid(i,1),mp]=max(M(1,i1:neSum));
    load_conditions=[load_conditions;
                    i Mleft(i,1) Mmid(i,1) Mright(i,1)]; %Kg-cm (flexure)
    
    %% Cut location ( local coordinates)
    cutxLocBeam(i,:)=cutLocationSSRecBeam(M(:,i1:neSum),dL);
    cutxLocMSBeams=[cutxLocMSBeams,cutxLocBeam(i,:)];
end

nspans=size(cutxLocBeam,1);

%% Materials and parameters for reinforcing steel
%% Concrete cover
hrec=50; % mm
brec=50; % lateral concrete cover

%% Materials
fy=500; % Yield stress of steel reinforcement (N/mm2)
wac=7.85e-6; % unit volume weight of the reinforcing steel (Kg/mm3)

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
                9 40;
                10 50]; % mm^2
 
dvs=10; % default diameter for stirrups 
pmin=0.003;
pmax=0.025;
hagg=20; % mm

%% Construction performance
WUB=[0.6,0.7];
WND=[0.6];
WNB=[0.4];
WNC=[0.6,0.7];
Wcs1=2;
Wcs2=0.5;
WfaConstr=[WUB,WND,WNB,WNC,Wcs1,Wcs2];

%% OPTIMAL DESIGN 
[volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db9Spans,EffSpans,...
MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
CFAspans]=OptimMSFSBeamsRebar3DSec(bv,hv,Lv,brec,hrec,hagg,pmin,pmax,...
rebarAvailable([3:10]',:),fcuv,load_conditions,fy,wac,cutxLocBeam,WfaConstr,1,1,...
[1,2],70,40);

if sum(db9Spans)==0
    return;
end

%% Shear design
dbLMRMSBeam=[];
ShearDesignMSBeam=[];
tenbLMRMSBeam=[];
totnbSpanMSBeams=[];
neSum=0;
for i = 1: nspans
    i1=neSum+1;
    
    neSum=neSum+ne(i);
    
    dbLMRMSBeam=[dbLMRMSBeam,db9Spans(i,:)];
    
    % Average percentage of cross-section reinforcement 
    rho=sum(tenbLMRspan(i,:).*(db9Spans(i,:).^2*pi/4))./(bv(i)*hv(i))/3;  
    
    % Execute shear design
    [s1(i,1),s2(i,1),s3(i,1),d1(i,1),d2(i,1)]=shearDesignBeams(Lv(i),bv(i),hv(i),...
                                            hrec,fcuv(i),fy,V(1,i1:neSum),dvs,rho);
    ShearDesignMSBeam=[ShearDesignMSBeam,...
                        s1(i,1),s2(i,1),s3(i,1),d1(i,1),d2(i,1),dvs];
    dvsBeams=dvs(1,1);
    
    % Decompose NbSpans matrix of size [ nspans x 9 ] into a vector of
    % size [ 1 x nspans * 9 ]       
    tenbLMRMSBeam=[tenbLMRMSBeam,tenbLMRspan(i,:)];
    totnbSpanMSBeams=[totnbSpanMSBeams,totnbSpan(i,:)];
end
ShearDesignMSBeams{1}=ShearDesignMSBeam;


beamNSb=zeros(nMSbeams,1);
%% Side rebars
if max(hv)>=750
    [dSb,nSb,sepSb,distrSideBars]=sideBarsRecBeams3SecSpan(max(bv),max(hv),fy,...
            brec,hrec,tenbLMRspan,db9Spans,dvs,hagg,rebarAvailable);
    
    beamNSb(1,1)=2*nSb;
    diamlistdSb=zeros(2*nSb,1)+dSb;
    
    plotBeamSideBar3sec(b,h,-DistrRebarLeft,ListRebarDiamLeft,...
            DistrRebarMid,ListRebarDiamMid,-DistrRebarRight,...
            ListRebarDiamRight,diamlistdSb,distrSideBars,nfig);
else
    distrSideBars=[];
    diamlistdSb=[];
end

diamlistdSbMSBeam{1}=diamlistdSb;
distrSbMSBeam{1}=distrSideBars;

%% Export results
if saveResults
    nMSbeams=1;
    eleMSBeams{1}=[1,2];
    NbMSBeams{1}=tenbLMRMSBeam;
    DimMSBeams{1}=[bv',hv',Lv',brec,hrec];
    dbMSBeams{1}=dbLMRMSBeam;
    DecompCutxLocMSBeams{1}=cutxLocMSBeams;
    ListDiamLeftMSBeams{1}=ListRebarDiamLeft;
    ListDiamMidMSBeams{1}=ListRebarDiamMid;
    ListDiamRightMSBeams{1}=ListRebarDiamRight;
    LenRebarLeftMSBeams{1}=LenRebarL;
    LenRebarMidMSBeams{1}=LenRebarM;
    LenRebarRightMSBeams{1}=LenRebarR;
    DistrDiamLeftMSBeams{1}=DistrRebarLeft;
    DistrDiamMidMSBeams{1}=DistrRebarMid;
    DistrDiamRightMSBeams{1}=DistrRebarRight;
    NTotalMSBeams{1}=totnbSpanMSBeams;
    
    ExportDesignMSRecBeam(directionData,nMSbeams,eleMSBeams,DimMSBeams,...
        dbMSBeams,DecompCutxLocMSBeams,DistrDiamLeftMSBeams,DistrDiamMidMSBeams,...
        DistrDiamRightMSBeams,NTotalMSBeams,NbMSBeams,ListDiamLeftMSBeams,...
        ListDiamMidMSBeams,ListDiamRightMSBeams,LenRebarLeftMSBeams,...
        LenRebarMidMSBeams,LenRebarRightMSBeams,diamlistdSbMSBeam,distrSbMSBeam,...
        beamNSb,ShearDesignMSBeams)
end

%------------------------------- end -----------------------------------