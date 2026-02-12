% Optimal_Design_SSFSRecBeams_Complete_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To design optimally (with respect to saving in reinforcing volumes)
%    a beam element for all its three critical cross-sctions (left,middle
%    right)
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%----------------------------------------------------------------

clc
clear all

A=importdata('');

DR=A.data;
nsamples=size(DR,1);

%% Geometry 
n1=1;
n2=nsamples;
arrayDataCollection=[];
for i=n1:n2
    fcu=DR(i,3);
    b=DR(i,1);
    h=DR(i,2);
    span=DR(i,4);

    Mleft=DR(i,5);
    Mmid=DR(i,6);
    Mright=DR(i,7);

    W1=DR(i,8);
    W2=DR(i,9);
    
    hrec=50; % 
    brec=50; % lateral concrete cover

    %% Materials
    fy=500; % Yield stress of steel reinforcement (N/mm2)
    wac=7.85e-6; % unit volume weight of the reinforcing steel (N/mm3)

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
    W(i,:)=[W1,W2]; % Uniformly Distributed Load over the whole beam  N / mm
    
    [R,U,V,M]=MSFSFEMBeams(span,Ac,Ic,Ec,supportsLoc,W(i,:),dL,wrange,0);

    %% Estimate cut locations
    % for a more proper computation of rebar volumes
    
    ne=(supportsLoc(2)-supportsLoc(1))/dL;
    
    Mleft=M(1,1);
    Mright=M(1,ne);
    [Mmid,mp]=max(M(1,1:ne));
    load_conditions=[1 Mleft Mmid Mright]; %Kg-cm (flexure)
    
    %% Cut location ( local coordinates)
    cutxLoc=cutLocationSSRecBeam(M(:,1:ne),dL);

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
                    9 40
                    10 50]; % mm^2

    dvs=10;
    pmin=0.003;
    pmax=0.025;

    %% OPTIMAL DESIGN 

    hagg=20;
    Wunb=[1.3,1.4,1];
    Wnd=[1.2,0.8];
    Wcut=[1.3,1.6,2];
    Wfac=[Wunb,Wnd,Wcut];

    nfig=1;
    [volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db9Spans,EffSpans,...
    MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
    DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
    CFAspans]=OptimMSFSBeamsRebarPar(b,h,span,brec,hrec,hagg,...
    pmin,pmax,rebarAvailable([2:10]',:),fcu,load_conditions,fy,wac,cutxLoc,Wfac,0,0,[]);

    arrayi=[b,h,fcu,span,Mleft,Mmid,Mright,W1,W2,tenbLMRspan,db9Spans];
    arrayDataCollection=[arrayDataCollection;
                         arrayi];

end
