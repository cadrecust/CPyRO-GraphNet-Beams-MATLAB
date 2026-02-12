function [bestPerformance,bestLenRebar,bestsepRebar,bestPosition,bestdbc9,...
    bestEfL,bestEffMid,bestEffRight,bestMrLeft,bestMrMid,bestMrRight,cbestLeft,...
    cBestMid,cBestRight,bestListRebarDiamLeft,bestListRebarDiamMid,...
    bestListRebarDiamRight,bestDistrRebarLeft,bestRebarDistrMid,...
    bestDistrRebarRight,bestnbtLMR,bestnbcut3sec,bestnblowRight,bestdblowRight,...
    bestCFA]=PSOBeamsRebar3DSec(b,h,span,Ccx,Ccy,hagg,pmin,pmax,rebarAvailable,...
    fcu,load_conditions,fy,wac,cutLoc,dbcc,nbcc,dblow,nblow,Wfac,graphConvPlot,ispan,...
	numberOfParticles,nMaxIter)

%------------------------------------------------------------------------
% Syntax:
% [c_best,bestMr,bestEf,best_area,tbest,h]=PSO3layerBeamsRebar3sec(b,h,duct,...
%    b_rec,h_rec,vSep,fc,Mu,fy,graphConvergencePlot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: N,mm
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
%                               (Kg/cm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%------------------------------------------------------------------------
Es=fy/0.00217;
dvs=10;

% Max and min reinforcement area
ndiam=length(rebarAvailable(:,1));

% Identifying existing rebar diameters in left section

indexLayerDiam=[];
for j=1:3
    for i=1:length(rebarAvailable(:,1))
        if dbcc(j)==rebarAvailable(i,2)
            indexLayerDiam=[indexLayerDiam,i];
        end
    end
end

% Identifying layers of rebar valid for optimization for left section

indexLayerOptim=zeros(1,3);
xmaxL=zeros(1,3);
xminL=zeros(1,3);
for i=1:3
    if nbcc(i)==0
        indexLayerOptim(i)=i;
        xmaxL(1,i)=ndiam;
        xminL(1,i)=1;
    else
        indexLayerOptim(i)=i;
        xmaxL(1,i)=indexLayerDiam(i);
        xminL(1,i)=indexLayerDiam(i);
    end
end

% Identifying existing rebar diameters in right section
indexLayerDiam=[];
for j=1:3
    for i=1:length(rebarAvailable(:,1))
        if dblow(j)==rebarAvailable(i,2)
            indexLayerDiam=[indexLayerDiam,i];
        end
    end
end

% Identifying layers of rebar valid for optimization for right section

indexLayerOptim=zeros(1,3);
xmaxR=zeros(1,3);
xminR=zeros(1,3);
for i=1:3
    if nblow(i)==0
        indexLayerOptim(i)=i;
        xmaxR(1,i)=ndiam;
        xminR(1,i)=1;
    else
        indexLayerOptim(i)=i;
        xmaxR(1,i)=indexLayerDiam(i);
        xminR(1,i)=indexLayerDiam(i);
    end
end

%% Number of rebars
dbmin=rebarAvailable(1,2);

[sepMindbmin,sepMax1m]=sepMinMaxHK13(dbmin,hagg,0);
nbmax=fix((b-2*Ccx-2*dvs+sepMindbmin)/(dbmin+sepMindbmin));

if nbcc(1)==0
    nbminl=2;
else
    nbminl=nbcc(1);
end

if nblow(1)==0
    nbminm=2;
else
    nbminm=nblow(1);
end

xmax=[xmaxL,xmaxR,nbmax,  nbmax,  nbmax,  nbmax,   nbmax,   nbmax];
xmin=[xminL,xminR,nbminl,nbcc(2),nbcc(3),nbminm,nblow(2),nblow(3)];

% Algorithmic parameters
alpha=1;
c1=2; % cognitive component
c2=2; % social component
dt=1;
inertiaWeight=1.3;
beta=0.99;

maxVelocity=(xmax-xmin)/dt;

%% Generate position and velocity vector of each particle
    
numberOfDimensionSpace=12;
PositionMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
velocityMatrix=zeros(numberOfParticles,numberOfDimensionSpace);
for i=1:numberOfParticles
    for j=1:numberOfDimensionSpace
        r=rand;
        PositionMatrix(i,j)=xmin(j)+fix(r*(xmax(j)-(xmin(j)-1)));
        velocityMatrix(i,j)=alpha/dt*(-(xmax(j)-xmin(j))*0.5+r*(xmax(j)-xmin(j)));
    end
end

cbestLeft=[];

bestPerformance=1e10; % Assumes non-negative fitness values!
bestPosition=zeros(1,numberOfDimensionSpace);

performance=zeros(numberOfParticles,1);
bestPositionSwarmMatrix=PositionMatrix;
bestPerformanceSwarm=inf(numberOfParticles,1);
for j=1:nMaxIter

    % Determine the best position and best performance

    for i=1:numberOfParticles
        position=PositionMatrix(i,:); % diameter combo + number of rebars
        
        db1l=rebarAvailable(position(1),2);
        db2l=rebarAvailable(position(2),2);
        db3l=rebarAvailable(position(3),2);
        db1m=rebarAvailable(position(4),2);
        db2m=rebarAvailable(position(5),2);
        db3m=rebarAvailable(position(6),2);
        
        [sepMin1l,~]=sepMinMaxHK13([db1l],hagg,0);
        [sepMin2l,~]=sepMinMaxHK13([db2l],hagg,0);
        [sepMin3l,~]=sepMinMaxHK13([db3l],hagg,0);
        
        [sepMin1m,~]=sepMinMaxHK13([db1m],hagg,0);
        [sepMin2m,~]=sepMinMaxHK13([db2m],hagg,0);
        [sepMin3m,~]=sepMinMaxHK13([db3m],hagg,0);
        
        sepMin=[sepMin1l,sepMin2l,sepMin3l,...
                sepMin1m,sepMin2m,sepMin3m];
        
        dbc=[db1l,db2l,db3l,db1m,db2m,db3m];
        
        nb1l=position(7);
        nb2l=position(8);
        nb3l=position(9);
        
		[isNbLFeasible]=continuityConstrNbSec(nbcc,[nb1l,nb2l,nb3l]);
		
        nb1m=position(10);
        nb2m=position(11);
        nb3m=position(12);
        
		[isNbMFeasible]=continuityConstrNbSec(nblow,[nb1m,nb2m,nb3m]);
		
        nblm=[nb1l,nb2l,nb3l,nb1m,nb2m,nb3m];
        
        if all([dbc(1)>=dbc(2),dbc(2)>=dbc(3),...
                dbc(4)>=dbc(5),dbc(5)>=dbc(6),...
                isNbLFeasible,isNbLFeasible])
			
			[performance(i),LenRebar{i},sepRebar(:,:,i),NbCombo9(i,:),EffLeft(i),EffMid(i),...
			EffRight(i),MrLeft(i),MrMid(i),MrRight(i),cLeft(i),cMid(i),cRight(i),ListRebarDiamLeft{i},...
			ListRebarDiamMid{i},ListRebarDiamRight{i},DistrRebarLeft{i},RebarDistrMid{i},...
			DistrRebarRight{i},dbcRight(i,:),nbcut3sec(:,:,i),nblowLeft{i},dblowLeft{i},nbTopMid{i},...
			dbTopMid{i},nblowRight{i},dblowRight{i},CFA(i),const(i)]=CutRedistrOptimRecBeam3DSec(load_conditions,fcu,...
			Es,h,b,span,dbc,hagg,Ccx,Ccy,pmin,pmax,sepMin,rebarAvailable,cutLoc,...
			Wfac,nblm);
		else 
			performance(i)=1e10;
			LenRebar{i}=0;
			CFA(i)=0;
			sepRebar(:,:,i)=zeros(3,3);
			EffLeft(i)=0;
			EffMid(i)=0;
			EffRight(i)=0;
			MrLeft(i)=0;
			MrMid(i)=0;
			MrRight(i)=0;
			cLeft(i)=0;
			cMid(i)=0;
			cRight(i)=0;
			ListRebarDiamLeft{i}=[];
			ListRebarDiamMid{i}=[];
			ListRebarDiamRight{i}=[];
			DistrRebarLeft{i}=[];
			RebarDistrMid{i}=[];
			DistrRebarRight{i}=[];
			NbCombo9(i,:)=zeros(1,9);
			nbcut3sec(:,:,i)=zeros(3,3);
			dbcRight(i,:)=zeros(1,3);
			nblowRight{i}=zeros(1,3);
			dblowRight{i}=zeros(1,3);
			const(i)=10;
        end
        dbc9(i,:)=[dbc,dbcRight(i,:)];
    %end
    %for i=1:numberOfParticles
		if (performance(i)<bestPerformance && const(i)==0)
            bestPerformance=performance(i); % rebar volume (mm^3)
            bestCFA=CFA(i);
            bestLenRebar=LenRebar{i};
            bestsepRebar=sepRebar(:,:,i);
            bestPosition=PositionMatrix(i,:);
            bestdbc9=dbc9(i,:);
            bestEfL=EffLeft(i);
            bestEffMid=EffMid(i);
            bestEffRight=EffRight(i);
            bestMrLeft=MrLeft(i);
            bestMrMid=MrMid(i);
            bestMrRight=MrRight(i);
            cbestLeft=cLeft(i);
            cBestMid=cMid(i);
            cBestRight=cRight(i);
            bestListRebarDiamLeft=ListRebarDiamLeft{i};
            bestListRebarDiamMid=ListRebarDiamMid{i};
            bestListRebarDiamRight=ListRebarDiamRight{i};
            bestDistrRebarLeft=DistrRebarLeft{i};
            bestRebarDistrMid=RebarDistrMid{i};
            bestDistrRebarRight=DistrRebarRight{i};
            bestnbtLMR=NbCombo9(i,:);
            bestnbcut3sec=nbcut3sec(:,:,i);
            bestnblowRight=nblowRight{i};
            bestdblowRight=dblowRight{i};
        end
        
        if (performance(i)<bestPerformanceSwarm(i))
            bestPerformanceSwarm(i)=performance(i);
            bestPositionSwarmMatrix(i,:)=PositionMatrix(i,:);
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

% Optima Convergence Graph
if graphConvPlot==1 && isempty(cbestLeft)==0
    figure(4)
    if sum(nbcc)==0
        plot(1:1:nMaxIter,...
           performanceHist.*wac,'-o')
        hold on
        legend('Span 1')
    else
        plot(1:1:nMaxIter,...
           performanceHist.*wac,'-o','DisplayName',strcat('Span',num2str(ispan)))
        hold on
    end
    xlabel('Iteration')
    ylabel('Rebar weight (Kg)')
    title({'Optimum convergence of reinforcement weight for a concrete beam';
            'PSO'})
    hold on
end

if isempty(cbestLeft)==1
    bestPerformance=1e10; % rebar volume (mm^3)
    bestLenRebar=0;
    bestCFA=0;
    bestsepRebar=zeros(3,3);
    bestPosition=[];
    bestdbc9=zeros(1,9);
    bestEfL=0;
    bestEffMid=0;
    bestEffRight=0;
    bestMrLeft=0;
    bestMrMid=0;
    bestMrRight=0;
    cbestLeft=0;
    cBestMid=0;
    cBestRight=0;
    bestListRebarDiamLeft=[];
    bestListRebarDiamMid=[];
    bestListRebarDiamRight=[];
    bestDistrRebarLeft=[];
    bestRebarDistrMid=[];
    bestDistrRebarRight=[];
    bestnbtLMR=zeros(1,9);
    bestnbcut3sec=[];
    
    bestnblowRight=zeros(1,3);
    bestdblowRight=zeros(1,3);
end