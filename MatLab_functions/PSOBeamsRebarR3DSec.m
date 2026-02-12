function [bestPerformance,bestEffR,bestMrR,bestcRight,bestnb,bestsepRebar,...
    bestListDiam,bestRebarDistr,bestdbc]=PSOBeamsRebarR3DSec...
    (Mur,fc,h,b,hagg,brec,hrec,pmin,pmax,sepRebarleft,distrRebarComp,...
    listRebarDiamComp,nbAfterCut3,nb3l,dbc,rebarAvailable)

%------------------------------------------------------------------------
% Syntax:
% [c_best,bestMr,bestEf,best_area,tbest,h]=PSO3layerBeamsRebar(b,h,duct,...
%    b_rec,h_rec,vSep,fc,Mu,fy,graphConvergencePlot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement area for a given beam 
% cross-section with specified initially dimensions (b,h) through the SGD 
% method.
% 
% OUTPUT: c_best,bestMr,bestEf: The neutral axis depth for the optimal 
%                               design, the resistant bending moment for 
%                               the optimal design,
%
%         best_area,tbest:      The optimal reinforcement area, the optimal 
%                               t width of the ISR
%
%         h:                    The final cross-section height in case it 
%                               is modified from the given initial proposal 
%                               value
%
% INPUT:  load_conditions:      vector as [nload,Mu] size: nloads x 2
%
%         factor_fc:            is determined by de applicable design code. 
%                               The ACI 318-19 specifies it as 0.85
%
%         duct:                 is the ductility demand parameter, with 
%                               possible values of 1,2 or 3, for low 
%                               ductility, medium ductility or high 
%                               ductility respectively
%
%         h_rec,b_rec:          is the concrete cover along the height dimension
%                               and the width cross-section dimension, respectively
%                               (cm)
%
%         h,b:                  cross-section dimensions (cm)
%
%         E:                    is the Elasticity Modulus of reinforcing steel
%                               (Kg/mm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%------------------------------------------------------------------------

Es=200e3;
dvs=10;
bp=b-2*brec-2*dvs-max(dbc);

% Identifying valid rebar diameters
ndiam=length(rebarAvailable(:,1));
indexLayerDiam=[];
for j=1:3
    for i=1:ndiam
        if dbc(j)==rebarAvailable(i,2)
            indexLayerDiam=[indexLayerDiam,i];
        end
    end
end

% Identifying valid layers of rebars for optimization of right section
indexLayerOptim=[];
for i=1:3
    if nbAfterCut3(i)==0
        indexLayerOptim=[indexLayerOptim,i];
    end
end

id=indexLayerDiam(indexLayerOptim); % available sublayers for optim.
ndoptim=length(id);
indexDiamLayBegin=3-ndoptim+1;

%% Number of rebars
dbmin=rebarAvailable(1,2);

[sepMindbmin,sepMax1m]=sepMinMaxHK13(dbmin,hagg,0);
nbmaxi=fix((b-2*brec-2*10+sepMindbmin)/(dbmin+sepMindbmin));

for i=1:3
    if nb3l(i)-nbAfterCut3(i)==0 && nbAfterCut3(i)>2 % if there were cuts
        nbmax(i)=0;
    else
        nbmax(i)=nbmaxi;
    end
end

xmax=[id,             nbmax(1), nbmax(2), nbmax(3)];
xmin=[ones(1,ndoptim),0,        0,        0];

% Algorithmic parameters
alpha=1;
c1=2; % cognitive component
c2=2; % social component
dt=1;
inertiaWeight=1.3;
beta=0.99;

maxVelocity=(xmax-xmin)/dt;

%% Generate position and velocity vector of each particle
numberOfParticles=50;
numberOfDimensionSpace=ndoptim+3;    

PositionMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
velocityMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
for i=1:numberOfParticles
    for j=1:numberOfDimensionSpace
        r=rand;
        PositionMatrix(i,j)=xmin(j)+fix(r*(xmax(j)-(xmin(j)-1)));
        velocityMatrix(i,j)=alpha/dt*(-(xmax(j)-xmin(j))*0.5+r*(xmax(j)-xmin(j)));
    end
end

nMaxIter=40;
iteration=0;
bestPerformance=1e10;
bestEffR=0;
bestPosition=zeros(1,numberOfDimensionSpace);

performance=zeros(numberOfParticles,1);
bestPositionSwarmMatrix=PositionMatrix;
bestPerformanceSwarm=inf(numberOfParticles,1);
for j=1:nMaxIter
    iteration=iteration+1;

    % Determine the best position and best performance

    for i=1:numberOfParticles
        position=PositionMatrix(i,:); % diameter combo
        for k=1:ndoptim
            dbc(indexDiamLayBegin-1+k)=rebarAvailable(position(1,k),2);
        end
        [sepMin1,sepMax1]=sepMinMaxHK13(dbc(1),hagg,0);
        [sepMin2,sepMax2]=sepMinMaxHK13(dbc(2),hagg,0);
        [sepMin3,sepMax3]=sepMinMaxHK13(dbc(3),hagg,0);

        sepMin=[sepMin1,sepMin2,sepMin3];
        nbr(1)=position(1,ndoptim+1)+nbAfterCut3(1);
        nbr(2)=position(1,ndoptim+2)+nbAfterCut3(2);
        nbr(3)=position(1,ndoptim+3)+nbAfterCut3(3);
        
        %% Rebar separation
        
        sepRebar1=(bp)/(nbr(1)-1);
        if nbr(2)>1
            sepRebar2=(bp)/(nbr(2)-1);
        else
            sepRebar2=1e5;
        end
        if nbr(3)>1
            sepRebar3=(bp)/(nbr(3)-1);
        else
            sepRebar3=1e5;
        end
        sepRebarRight=[sepRebar1,sepRebar2,sepRebar3];
        
        for ii=1:3
            if nb3l(ii)-nbAfterCut3(ii)==0 && nbAfterCut3(ii)>2
                sepRebarRight(ii)=sepRebarleft(ii);
            end
        end
        
		[isNbRFeasible]=continuityConstrNbSec(nbAfterCut3,nbr);
		if isNbRFeasible
			[Abr,EffR,MrR,cRight,xBest,ListDiam,RebarDistr,isfeasibleRight]=...
			PlainDistrEffConstrR3DSec(Mur,fc,Es,h,b,hagg,brec,hrec,pmin,pmax,sepMin,...
			sepRebarRight,nbAfterCut3,distrRebarComp,listRebarDiamComp,nb3l,dbc,nbr);
		else
			Abr=1e10;
			EffR=0;
			MrR=0;
			cRight=0;
			xBest=zeros(1,3);
			ListDiam=0;
			RebarDistr=[];
			isfeasibleRight=false;
		end
        if (Abr<bestPerformance && isfeasibleRight)
            
            bestdbc=dbc;
            bestPerformance=Abr;
            bestEffR=EffR;
            bestMrR=MrR;
            bestcRight=cRight;
            bestnb=xBest;
            bestsepRebar=sepRebarRight;
            bestListDiam=ListDiam;
            bestRebarDistr=RebarDistr;
        end
        if (performance(i)<bestPerformanceSwarm(i) && EffR>0)
            bestPerformanceSwarm(i)=performance(i);
            bestPositionSwarmMatrix(i,:)=position;
        end
    end
    performanceHist(j)=bestPerformance;
    
    %% Global best position
    globalBestPosition=bestPosition;
    globalBestPerformance=bestPerformance;
    
    %% Velocity and position uptdating 

    for i=1:numberOfParticles
        q=rand;
        r=rand;
        %% Velocity updating
        for j=1:numberOfDimensionSpace
            velocityMatrix(i,j)=inertiaWeight*velocityMatrix(i,j)+...
                c1*q*((bestPositionSwarmMatrix(i,j)-PositionMatrix(i,j))/dt)+...
                c2*r*((bestPosition(j)-PositionMatrix(i,j)));
            
            absolouteValueVelocity(j)=abs(velocityMatrix(i,j));
            valueVelocity(j)=velocityMatrix(i,j);
            
            if (absolouteValueVelocity(j)>maxVelocity(j))
                    velocityMatrix(i,j)=maxVelocity(j);
            end
            if (valueVelocity(j)<-maxVelocity(j))
                    velocityMatrix(i,j)=-maxVelocity(j);
            end
        end

        %% Position uptdating

        for j=1:numberOfDimensionSpace
            tj=fix(PositionMatrix(i,j)+velocityMatrix(i,j)*dt);
            if tj<xmin(j)
                tj=xmin(j);
            elseif tj>xmax(j)
                tj=xmax(j);
            end
            PositionMatrix(i,j)=tj;
        end
    end
    
    if(inertiaWeight*beta<0.3)
        inertiaWeight=inertiaWeight;
    else
        inertiaWeight=inertiaWeight*beta;
    end
end

if bestEffR==0
    bestPerformance=1e10;
    bestEffR=0;
    bestMrR=0;
    bestcRight=0;
    bestnb=zeros(1,3);
    bestsepRebar=zeros(1,3);
    bestListDiam=0;
    bestRebarDistr=[];
    bestdbc=zeros(1,3);
end