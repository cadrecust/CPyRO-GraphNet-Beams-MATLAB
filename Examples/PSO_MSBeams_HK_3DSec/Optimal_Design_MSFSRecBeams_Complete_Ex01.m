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

dL=100;

bv=[300;300];
hv=[600;600];
fcuv=[35;35];
Lv=[4000;5000];

Totalspan = 9000;
Wv=[20,20;
    30,30];

supportsv=[0,4000,9000];
Wrange=[0, 4000;
        4000, 9000];
nspans=length(Wv(:,1));

Ac = bv .* hv;
Ic = bv .* hv .^ 3 / 12;
Ec = (3.46 * sqrt(fcuv) + 3.21) .* 1e3;

[R,U,V,M]=MSFSFEMBeams(Totalspan,Ac,Ic,Ec,supportsv,Wv,dL,Wrange,0);

ne=zeros(1,nspans);
neSum=0;
load_conditions=[];
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
    
end

nBeams=size(cutxLocBeam,1);

%% Materials and parameters for reinforcing steel
%% Concrete cover
hrec=50; % 
brec=50; % lateral concrete cover

%% Materials
fy=500; % Yield stress of steel reinforcement (N/mm2)
wac=7.85e-6; % unit volume weight of the reinforcing steel (N/mm3)

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

dvs=10;
pmin=0.003;
pmax=0.025;
hagg=20;

%% Construction performance
ucfactors=zeros(1,11);
ucwire=4.5;
uclinks=5;
ucrebar=6;
quantwire=0.05;
quantlinks=0.105;

ucLorry=880;
quantLorry=1/1000;

perfManCutBend=1/600;
ucManCutBend=1945;
ucManTyeAssem=2091;
perfManTyeAssem=1/200;

ucfactors(1)=ucwire;
ucfactors(2)=uclinks;
ucfactors(3)=ucrebar;
ucfactors(7)=quantwire;
ucfactors(8)=quantlinks;
ucfactors(4)=ucManCutBend;
ucfactors(6)=ucLorry;
ucfactors(9)=quantLorry;
ucfactors(10)=perfManCutBend;
ucfactors(5)=ucManTyeAssem;
ucfactors(11)=perfManTyeAssem;

WUB=[0.6,0.7];
WND=[0.6];
WNB=[0.4];
WNC=[0.6,0.7];
Wcs1=2;
Wcs2=0.5;
WfaConstr=[WUB,WND,WNB,WNC,Wcs1,Wcs2];

fcu=fcuv;
span=Lv;

b=bv;
h=hv;

%% OPTIMAL DESIGN 
[volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db9Spans,EffSpans,...
MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
CFAspans]=OptimMSFSBeamsRebar3DSec(b,h,span,brec,hrec,hagg,pmin,pmax,...
rebarAvailable,fcu,load_conditions,fy,wac,cutxLocBeam,WfaConstr,1,1,...
[1,2],30,50);

for i=1:nBeams
    bestAbi2=db9Spans.^2*pi/4;
    bestAbTL1=sum(tenbLMRspan(1,1:3).*bestAbi2(1,1:3));
    bestAbBM1=sum(tenbLMRspan(1,4:6).*bestAbi2(1,4:6));
    bestAbTR1=sum(tenbLMRspan(1,7:9).*bestAbi2(1,7:9));

    Ast3(i,:)=[bestAbTL1,bestAbBM1,bestAbTR1];
end

%% Function appendix

function [r,u,esbarsShear,esbarsMoment]=MSFSFEMBeams(L,Az,Iz,Ee,...
        supportsLoc,w,dL,wrange,plMVdiag)

%% Geometry
np=fix(L/dL)+1;
nnodes=np;
nbars=np-1;
nspans=length(supportsLoc)-1;
A=zeros(nbars,1);
I=zeros(nbars,1);
E=zeros(nbars,1);
neSum=0;
for i=1:nspans
    nbarsi=fix((supportsLoc(i+1)-supportsLoc(i))/dL);
    i1=neSum+1;
    neSum=neSum+nbarsi;
    A(i1:neSum)=Az(i);
    I(i1:neSum)=Iz(i);
    E(i1:neSum)=Ee(i);
end

% Mesh
coordxy(1:np,1)=0:dL:L;
coordxy(1:np,2)=0;

% Node connectivity
ni(1:nbars,1)=1:1:(np-1);
nf(1:nbars,1)=2:1:nnodes;

%% Boundary conditions
nsupports=length(supportsLoc);
for i=1:nsupports
    nodeSupports(i)=fix(supportsLoc(i)/dL)+1;
end

%% Topology
% Fixed supports at left end

ndofSupports(1)=nodeSupports(1)*3-2;
ndofSupports(2)=nodeSupports(1)*3-1;
ndofSupports(3)=nodeSupports(1)*3;
    
% Simply supported in between ends
for i=2:length(nodeSupports)-1
    ndofSupports(i*3-2)=nodeSupports(i)*3-2;
    ndofSupports(i*3-1)=nodeSupports(i)*3-1;
end

% Fixed supports at right end
ndofSupports(nsupports*3-2)=nodeSupports(nsupports)*3-2;
ndofSupports(nsupports*3-1)=nodeSupports(nsupports)*3-1;
ndofSupports(nsupports*3)=nodeSupports(nsupports)*3;

ndofs=length(ndofSupports);
ndofSupports2=[];
for i=1:ndofs
    if  ndofSupports(i)>0
        ndofSupports2=[ndofSupports2,ndofSupports(i)];
    end
end
ndofSupports=ndofSupports2;

bc(:,1)=ndofSupports';
bc(:,2)=0;

Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end

%% Loads
qbarray(1:nbars,1)=1:1:nbars;

NDistLoads=length(wrange(:,1));
for i=1:NDistLoads
    ew1=fix(wrange(i,1)/dL)+1;
    ew2=ceil(wrange(i,2)/dL);
    
    W1=w(i,1);
    W2=w(i,2);
    dW=(W2-W1)/(np-1);
    k=1;
    for j=ew1:ew2
        qbarray(j,2)=w(i,1)+(k-1)*dW;
        k=k+1;
    end
end

Kglobal=zeros(3*nnodes);
fglobal=zeros(3*nnodes,1);

%% Matrix assembly and solver
ep_bars=zeros(nbars,3); 
eq_bars=zeros(nbars,2);
for i=1:nbars 
    ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
    ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
    ep=[E(i) A(i) I(i)];
    eq=[0 -qbarray(i,2)];

    ep_bars(i,:)=ep;
    eq_bars(i,:)=eq;
    [Ke_barra,fe_barra]=beam2e(ex,ey,ep,eq);

    [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Ke_barra,fglobal,fe_barra);

end
[u,r]=solveq(Kglobal,fglobal,bc); % solving system of equations

Ed=extract(Edof,u);
ex=coordxy(:,1);
ey=coordxy(:,2);

Ex=zeros(nbars,2);
Ey=zeros(nbars,2);

for j=1:nbars
    Ex(j,1)=ex(Edof(j,4)/3);
    Ex(j,2)=ex(Edof(j,7)/3);

    Ey(j,1)=ey(Edof(j,4)/3);
    Ey(j,2)=ey(Edof(j,7)/3);

end

%% Forces diagrams
nfep=2;
esbarsNormal=zeros(nfep,nbars);
esbarsShear=zeros(nfep,nbars);
esbarsMoment=zeros(nfep,nbars);
for i=1:nbars
    es_bar=beam2s(Ex(i,:),Ey(i,:),ep_bars(i,:),Ed(i,:),eq_bars(i,:),nfep);
    esbarsNormal(:,i)=es_bar(:,1);
    esbarsShear(:,i)=es_bar(:,2);
    esbarsMoment(:,i)=es_bar(:,3);
end

if plMVdiag==1
    
    %----Undeformed mesh-----%
    figure(6)
    xlabel('X [m]')
    ylabel('Y [m]')
    title('Deformed structure');
    plotpar=[2 1 0];
    eldraw2(Ex,Ey,plotpar);
    hold on
    
    %-----Deformed mesh-----%
    figure(6)
    plotpar=[1 2 1];
    eldisp2(Ex,Ey,Ed,plotpar,100);  
    hold on
    
    sfac=scalfact2(Ex(1,:),Ey(1,:),esbarsShear(:,1),1);
    for i=1:nbars
        figure(4)
        plotpar=[2 1];
        eldia2(Ex(i,:),Ey(i,:),esbarsShear(:,i),plotpar,sfac*10);
        title('Shear Force')
    end
    
    sfac=scalfact2(Ex(1,:),Ey(1,:),esbarsMoment(:,1),1);
    for i=1:nbars
         figure(5)
         plotpar=[4 1];
         eldia2(Ex(i,:),Ey(i,:),esbarsMoment(:,i),plotpar,sfac*10);
         title('Bending Moment')
         xlabel('X [m]')
         ylabel('Y [KN-m]')
    end
end
end
%------------------------------- end -----------------------------------