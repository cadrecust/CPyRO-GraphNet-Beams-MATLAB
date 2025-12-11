% Optim_Design_2SFSpan_Beams_Simple_GNNvsCPyRO_1MLOT_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To obtain the optimum PF from the NSGA-II assisted by
%    CPyRO-GraphNet-Beams and plain GNN, along with the PF evolution in
%    time, to compare performance of both approaches.
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

bv=[250;250];
hv=[500;500];
fcuv=[30;30];
Lv=[4000;5000];

Totalspan = 9000;
Wv=[35,35;
    40,40];

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
cutLoc=[];
for i=1:nspans
    
    ne(i)=(supportsv(i+1)-supportsv(i))/dL;
    i1=neSum+1;
    Mleft=M(1,i1);
    
    neSum=neSum+ne(i);
    Mright=M(1,neSum);
    
    [Mmid,mp]=max(M(1,i1:neSum));
    load_conditions=[load_conditions;
                    i Mleft Mmid Mright]; %Kg-cm (flexure)
    
    %% Cut location ( local coordinates)
    cutxLocBeam(i,:)=cutLocationSSRecBeam(M(:,i1:neSum),dL);
    
    xMmid = (mp-1) * dL ;

    %% Location of design moments Mu
    x0L(:,i)=[10,xMmid,Lv(i,1)-10]';
    
    %% Data for PIGNN
    DR(i,:)=[bv(i),hv(i),fcuv(i),Lv(i),Mleft,Mmid,Mright];
end

nBeams=size(x0L,2);

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

%% Constructability
hagg=20;

WUNb=[0.5, 0.6];
WNb=[0.5];
WND=[0.5]; 
Wcut=[0.7, 0.8]; 

Wcs1 = 0.8;
Wcs2 = 1.5;

Wfac=[WUNb,WND,WNb,Wcut,Wcs1,Wcs2];

%% Generalization parameters
Nv=[nBeams];
ndl=length(Nv);

N=Nv(1);

[idxTrain,idxValidation,idxTest] = trainingPartitions(N,[0 0 1]);
nTrain=length(idxTrain);
nValidation=length(idxValidation);
nTest=length(idxTest);

cutxLocBeamTrain=cutxLocBeam(idxTrain,:);
cutxLocBeamValidation=cutxLocBeam(idxValidation,:);
cutxLocBeamTest=cutxLocBeam(idxTest,:);

% Data for GCNN
meanX=[349.605734767025	633.225806451613	37.4534050179212	4505.55555555556	-26339479.7363044];
sigsqX=[4949.66534345653	11216.8343161059	31.3195134954590	996295.300677007	2.10064997682719e+15];

[XTrain,XValidation,XTest,ATrain,AValidation,ATest]=dataTrainNLClass...
    (DR,idxTrain,idxValidation,idxTest,meanX,sigsqX);

XTestOpt = DR(idxTest,1:7);

%% Load surrogate models

% Model for prediction of As per cross-section

nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_4000.mat");
paramPIGCNN=load("PIGCNN_As_Section_4000.mat");

Ao3CPyRO=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,ATest,nheadsparamnGATPIGNN.numHeads);
Ao3CPyRO=extractdata(Ao3CPyRO);

nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
paramGCNN=load("GCNN_As_Section_4000.mat");

Ao3GNN=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,ATest,nheadsparamnGATGNN.numHeads);
Ao3GNN=extractdata(Ao3GNN);

muY1=1000.5;
muY2=766.6692;
muY3=851.6936;

sigsqY1=2.0948e5;
sigsqY2=8.0557e4;
sigsqY3=1.5316e5;

% Denormalizing targets and outputs
ntest=length(idxTest);
for i=1:ntest
    i1=(i-1)*3+1;
    i2=i*3;
    
    Ao3GNN(i1,1)=Ao3GNN(i1,1)*sqrt(sigsqY1) + muY1;
    Ao3GNN(i1+1,1)=Ao3GNN(i1+1,1)*sqrt(sigsqY2) + muY2;
    Ao3GNN(i2,1)=Ao3GNN(i2,1)*sqrt(sigsqY3) + muY3;
    
end
fcu=XTestOpt(:,3);
span=XTestOpt(:,4);

b=XTestOpt(:,1);
h=XTestOpt(:,2);

%% Loads
Mleft=XTestOpt(:,5);
Mmid=XTestOpt(:,6);
Mright=XTestOpt(:,7);

%% OPTIMAL DESIGN 

nbcc=zeros(1,3);
nblow=zeros(1,3);
dbcc=zeros(1,3);
dblow=zeros(1,3);

PF_REF_WEIGHT=[0];
PF_CS_REF=[1];

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

%% Optimization
MIGDconv=0;
pop_size=60;           % Population size
gen_max=80;            % MAx number of generations - stopping criteria

ngen1=ceil(gen_max/4);
ngen2=ceil(gen_max/2);
ngen3=ceil(3*gen_max/4);
ngen4=gen_max;

[extrOptPFCS1,PF_CFA1,PF_VOL1,newPop1,feasibleSol1,genCFA1,genVOL1,...
IGDt1,IGDv1]=SANSGAIIMSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable([2:10]',:),...
fcu,load_conditions,fy,wac,cutxLocBeamTest,dbcc,nbcc,dblow,nblow,PF_REF_WEIGHT,...
PF_CS_REF,Wfac,Ao3CPyRO,ucfactors,MIGDconv,pop_size,gen_max);

model2Use='CPyRO Graph-Net';
plotEvolutionPF(PF_CS_REF,PF_REF_WEIGHT,genCFA1,genVOL1,IGDt1,...
                IGDv1,feasibleSol1,wac,ngen1,ngen2,ngen3,gen_max,model2Use,1,11);

[extrOptPFCS2,PF_CFA2,PF_VOL2,newPop2,feasibleSol2,genCFA2,genVOL2,...
IGDt2,IGDv2]=SANSGAIIMSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable([2:10]',:),...
fcu,load_conditions,fy,wac,cutxLocBeamTest,dbcc,nbcc,dblow,nblow,PF_REF_WEIGHT,...
PF_CS_REF,Wfac,Ao3GNN,ucfactors,MIGDconv,pop_size,gen_max);

model2Use='Plain GCNN';
plotEvolutionPF(PF_CS_REF,PF_REF_WEIGHT,genCFA2,genVOL2,IGDt2,...
                IGDv2,feasibleSol2,wac,ngen1,ngen2,ngen3,gen_max,model2Use,1,12);

%% Analysis of PF Dominance
[CI,PFdom]=dominancePFs(PF_CFA1,PF_VOL1,PF_CFA2,PF_VOL2);

%% Function appendix

function [CI,PFdom]=dominancePFs(CS_PF1,WEIGHT_PF1,CS_PF2,WEIGHT_PF2)
    np1=length(CS_PF1);
    np2=length(CS_PF2);
    
    CI=[];
    PFdom1=0;
    PFdom2=0;
    for i=1:np1
        cip1=0;
        for j=1:np2
            if CS_PF2(j)>CS_PF1(i) && WEIGHT_PF2(j)<WEIGHT_PF1(i)
                PFdom2=PFdom2+1;
                cip1=1;
            else
                PFdom1=PFdom1+1;
            end
        end
        CI=[CI,cip1];
    end
    CI=sum(CI)/np1;

    if PFdom1>PFdom2
        PFdom=1;
    elseif PFdom1<PFdom2

        PFdom=2;
    else
        PFdom=0;
    end
end

function plotEvolutionPF(CFA_PF_REF,WEIGHT_PF_REF,gen_cfa,gen_vol,IGDt,...
    IGDv,feasibleSol,wac,nstep1,nstep2,nstep3,nstep4,surrogate,beamNo,nfig)
%% Baseline pareto for mIGD

% PARETO FRONT COST-CFA
figure(nfig)
xlabel('Constructability Score of Rebar Designs (CS)')
plot(CFA_PF_REF,WEIGHT_PF_REF,'r -',...
    'linewidth',1.8,'displayName',strcat('Reference PF'))
hold on
set(gca, 'Fontname', 'Times New Roman','FontSize',32);

refvol=1e10;
refcfa=0;
[np,ng]=size(gen_vol);
for i=1:ng
    for j=1:np
        if gen_vol(j,i)~=1e10
            refvol=gen_vol(j,i);
        else
            gen_vol(j,i)=refvol;
        end
        if gen_cfa(j,i)~=0 
            refcfa=gen_cfa(j,i);
        else
            gen_cfa(j,i)=refcfa;
        end
    end
end
pastel_gray = [0.663,0.663,0.663]; % #D3D3D3 for scatter points

figure(nfig)
title({strcat(surrogate,'-',' NSGAII');strcat('Double-Span Beam ',num2str(beamNo),'-RDP: Simple')})
xlabel('CS of Rebar Designs (CS)')
ylabel('Rebar Weight (Kgf)')
plot(gen_cfa(:,nstep1),gen_vol(:,nstep1).*wac,'o','linewidth',2.5,'color','#80B3FF','markerFaceColor','#80B3FF','markerSize',9.0)
hold on
plot(gen_cfa(:,nstep2),gen_vol(:,nstep2).*wac,'x','linewidth',2.0,'color','#77AC30','markerSize',13.0)
hold on
plot(gen_cfa(:,nstep3),gen_vol(:,nstep3).*wac,'s','linewidth',2.0,'color','magenta','markerSize',14.0)
hold on
plot(gen_cfa(:,nstep4),gen_vol(:,nstep4).*wac,'k ^','linewidth',1.5,'color','black','markerSize',13)
hold on
plot(1-feasibleSol(:,2),feasibleSol(:,1).*wac,'o','linewidth',0.4,'color',pastel_gray,'markerSize',5)
hold on
legend(strcat('Reference PF'),...
        strcat('PF-Gen-',num2str(nstep1),' (Time: ',num2str(IGDt(nstep1,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep2),' (Time: ',num2str(IGDt(nstep2,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep3), ' (Time: ',num2str(IGDt(nstep3,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep4), ' (Time: ',num2str(IGDt(nstep4,1)),' sec )'),...
        'Feasible Solutions',...
        'Location','northwest')
grid on
set(gca, 'Fontname', 'Times New Roman','FontSize',22);

end

function [XTrain,XValidation,XTest,ATrain,AValidation,ATest]=dataTrainNLClass...
    (DR,idxTrain,idxValidation,idxTest,meanX,sigsqX)

    sigsqX1=sigsqX(1);
    sigsqX2=sigsqX(2);
    sigsqX3=sigsqX(3);
    sigsqX4=sigsqX(4);
    sigsqX5=sigsqX(5);

    muX1=meanX(1);
    muX2=meanX(2);
    muX3=meanX(3);
    muX4=meanX(4);
    muX5=meanX(5);

    numObservations=length(DR(:,1));

    elements=[2 1;
              1 3];

    % Adjacency matrix
    numNodesGNN=3;
    % Adjancency matrix
    adjacency = zeros(numNodesGNN);
    for i = 1:size(elements,2)
        % The following logic specifies each node in an element is connected to
        % each other node in that element.
        nodesForElement = elements(:,i);
        for node = nodesForElement
            adjacency(nodesForElement,node) = 1;
        end
    end

    adjacency=repmat(adjacency,[1,1,numObservations]);

    X = [DR(:,1:7)];

    features1=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,5)];
    features2=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,6)];
    features3=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,7)];

    coulombData1=zeros(numObservations,3,3);
    coulombData2=zeros(numObservations,3,3);
    coulombData3=zeros(numObservations,3,3);
    coulombData4=zeros(numObservations,3,3);
    coulombData5=zeros(numObservations,3,3);
    for i=1:numObservations
        features=[features1(i,:)',features2(i,:)',features3(i,:)'];

        for j=1:numNodesGNN
            coulombData1(i,j,j)=features(1,j);
            coulombData2(i,j,j)=features(2,j);
            coulombData3(i,j,j)=features(3,j);
            coulombData4(i,j,j)=features(4,j);
            coulombData5(i,j,j)=features(5,j);
        end
    end

    coulombData1 = double(permute(coulombData1, [2 3 1]));
    coulombData2 = double(permute(coulombData2, [2 3 1]));
    coulombData3 = double(permute(coulombData3, [2 3 1]));
    coulombData4 = double(permute(coulombData4, [2 3 1]));
    coulombData5 = double(permute(coulombData5, [2 3 1]));

    %% Partition of data

    % node adjacency data
    adjacencyDataTrain = adjacency(:,:,idxTrain);
    adjacencyDataValidation = adjacency(:,:,idxValidation);
    adjacencyDataTest = adjacency(:,:,idxTest);

    % feature data
    coulombDataTrain1 = coulombData1(:,:,idxTrain);
    coulombDataValidation1 = coulombData1(:,:,idxValidation);
    coulombDataTest1 = coulombData1(:,:,idxTest);

    coulombDataTrain2 = coulombData2(:,:,idxTrain);
    coulombDataValidation2 = coulombData2(:,:,idxValidation);
    coulombDataTest2 = coulombData2(:,:,idxTest);

    coulombDataTrain3 = coulombData3(:,:,idxTrain);
    coulombDataValidation3 = coulombData3(:,:,idxValidation);
    coulombDataTest3 = coulombData3(:,:,idxTest);

    coulombDataTrain4 = coulombData4(:,:,idxTrain);
    coulombDataValidation4 = coulombData4(:,:,idxValidation);
    coulombDataTest4 = coulombData4(:,:,idxTest);

    coulombDataTrain5 = coulombData5(:,:,idxTrain);
    coulombDataValidation5 = coulombData5(:,:,idxValidation);
    coulombDataTest5 = coulombData5(:,:,idxTest);


    % Train partition

    [ATrain,XTrain1] = preprocessData(adjacencyDataTrain,coulombDataTrain1);
    [~,XTrain2] = preprocessData(adjacencyDataTrain,coulombDataTrain2);
    [~,XTrain3] = preprocessData(adjacencyDataTrain,coulombDataTrain3);
    [~,XTrain4] = preprocessData(adjacencyDataTrain,coulombDataTrain4);
    [~,XTrain5] = preprocessData(adjacencyDataTrain,coulombDataTrain5);

    % Validation partition
    [AValidation,XValidation1] = preprocessData(adjacencyDataValidation,coulombDataValidation1);
    [~,XValidation2] = preprocessData(adjacencyDataValidation,coulombDataValidation2);
    [~,XValidation3] = preprocessData(adjacencyDataValidation,coulombDataValidation3);
    [~,XValidation4] = preprocessData(adjacencyDataValidation,coulombDataValidation4);
    [~,XValidation5] = preprocessData(adjacencyDataValidation,coulombDataValidation5);

    %% Normalizing training data

    XTrain1 = (XTrain1 - muX1)./sqrt(sigsqX1);
    XTrain2 = (XTrain2 - muX2)./sqrt(sigsqX2);
    XTrain3 = (XTrain3 - muX3)./sqrt(sigsqX3);
    XTrain4 = (XTrain4 - muX4)./sqrt(sigsqX4);
    XTrain5 = (XTrain5 - muX5)./sqrt(sigsqX5);

    XTrain=[XTrain1,XTrain2,XTrain3,XTrain4,XTrain5];

    %% Normalizing validation data
    XValidation1 = (XValidation1 - muX1)./sqrt(sigsqX1);
    XValidation2 = (XValidation2 - muX2)./sqrt(sigsqX2);
    XValidation3 = (XValidation3 - muX3)./sqrt(sigsqX3);
    XValidation4 = (XValidation4 - muX4)./sqrt(sigsqX4);
    XValidation5 = (XValidation5 - muX5)./sqrt(sigsqX5);

    XValidation=[XValidation1,XValidation2,XValidation3,XValidation4,XValidation5];

    %% Normalizing test data
    [ATest,XTest1] = preprocessData(adjacencyDataTest,coulombDataTest1);
    XTest1 = (XTest1 - muX1)./sqrt(sigsqX1);
    XTest1 = dlarray(XTest1);

    [~,XTest2] = preprocessData(adjacencyDataTest,coulombDataTest2);
    XTest2 = (XTest2 - muX2)./sqrt(sigsqX2);
    XTest2 = dlarray(XTest2);

    [~,XTest3] = preprocessData(adjacencyDataTest,coulombDataTest3);
    XTest3 = (XTest3 - muX3)./sqrt(sigsqX3);
    XTest3 = dlarray(XTest3);

    [~,XTest4] = preprocessData(adjacencyDataTest,coulombDataTest4);
    XTest4 = (XTest4 - muX4)./sqrt(sigsqX4);
    XTest4 = dlarray(XTest4);

    [~,XTest5] = preprocessData(adjacencyDataTest,coulombDataTest5);
    XTest5 = (XTest5 - muX5)./sqrt(sigsqX5);
    XTest5 = dlarray(XTest5);

    XTest=[XTest1,XTest2,XTest3,XTest4,XTest5];
end

function Y = PlainGNNmodel1fc1GAT1Conv1fc(parameters,X,A,numHeads)

    Z1 = X * parameters.Embedding.Weights + parameters.Embedding.b;

    weights1 = parameters.attn1.Weights;
    numHeadsAttention1 = numHeads.attn1;
    
    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = relu(Z2);

    ANorm = normalizeAdjacency(A);
    Z3 = single(full(ANorm)) * Z2 * double(parameters.mult1.Weights);
    Z3  = relu(Z3) + Z2;

    Z4 = Z3 * parameters.Decoder.Weights + parameters.Decoder.b;
    
    Y = Z4;
end

function Y = PIGNNmodel1fc1GAT1Conv1fc(parameters,X,A,numHeads)

    Z1 = X * parameters.Embed.Weights + parameters.Embed.b;

    weights1 = parameters.attn1.Weights;
    numHeadsAttention1 = numHeads.attn1;
    
    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = relu(Z2);

    ANorm = normalizeAdjacency(A);
    Z3 = single(full(ANorm)) * Z2 * double(parameters.mult1.Weights);
    Z3  = relu(Z3) + Z2;

    Z4 = Z3 * parameters.Decoder.Weights + parameters.Decoder.b;
    
    Y = Z4;
end

function [outputFeatures,normAttentionCoeff] = graphAttention(inputFeatures,...
    adjacency,weights,numHeads,aggregation)
    
    % Split weights with respect to the number of heads and reshape the matrix to a 3-D array
    szFeatureMaps = size(weights.linearWeights);
    numOutputFeatureMapsPerHead = szFeatureMaps(2)/numHeads;
    linearWeights = reshape(weights.linearWeights,[szFeatureMaps(1), numOutputFeatureMapsPerHead, numHeads]);
    attentionWeights = reshape(weights.attentionWeights,[numOutputFeatureMapsPerHead, 2, numHeads]);
    
    % Compute linear transformations of input features
    value = pagemtimes(inputFeatures,linearWeights);
    
    % Compute attention coefficients
    query = pagemtimes(value, attentionWeights(:, 1, :));
    key = pagemtimes(value, attentionWeights(:, 2, :));
    
    attentionCoefficients = query + permute(key,[2, 1, 3]);
    attentionCoefficients = leakyrelu(attentionCoefficients,0.2);
    
    % Compute masked attention coefficients
    mask = -10e9 * (1 - adjacency);
    attentionCoefficients = attentionCoefficients + mask;
    
    % Compute normalized masked attention coefficients
    normAttentionCoeff = softmax(attentionCoefficients,DataFormat = "BCU");
    
    % Normalize features using normalized masked attention coefficients
    headOutputFeatures = pagemtimes(normAttentionCoeff,value);
    
    % Aggregate features from multiple heads
    if strcmp(aggregation, "cat")
        outputFeatures = headOutputFeatures(:,:);
    else
        outputFeatures =  mean(headOutputFeatures,3);
    end

end


function varargout = trainingPartitions(numObservations,splits)
	%TRAININGPARTITONS Random indices for splitting training data
	%   [idx1,...,idxN] = trainingPartitions(numObservations,splits) returns
	%   random vectors of indices to help split a data set with the specified
	%   number of observations, where SPLITS is a vector of length N of
	%   partition sizes that sum to one.
	%
	%   % Example: Get indices for 50%-50% training-test split of 500
	%   % observations.
	%   [idxTrain,idxTest] = trainingPartitions(500,[0.5 0.5])
	%
	%   % Example: Get indices for 80%-10%-10% training, validation, test split
	%   % of 500 observations. 
	%   [idxTrain,idxValidation,idxTest] = trainingPartitions(500,[0.8 0.1 0.1])
	%{
	arguments
		numObservations (1,1) {mustBePositive}
		splits {mustBeVector,mustBeInRange(splits,0,1,"exclusive"),mustSumToOne}
	end
	%}
	numPartitions = numel(splits);
	varargout = cell(1,numPartitions);

	idx = randperm(numObservations);

	idxEnd = 0;

	for i = 1:numPartitions-1
		idxStart = idxEnd + 1;
		idxEnd = idxStart + floor(splits(i)*numObservations) - 1;

		varargout{i} = idx(idxStart:idxEnd);
	end

	% Last partition.
	varargout{end} = idx(idxEnd+1:end);

end

function mustSumToOne(v)
    % Validate that value sums to one.

    if sum(v,"all") ~= 1
        error("Value must sum to one.")
    end

end


function ANorm = normalizeAdjacency(A)

	% Add self connections to adjacency matrix.
	A = A + speye(size(A));

	% Compute inverse square root of degree.
	degree = sum(A, 2);
	degreeInvSqrt = sparse(sqrt(1./degree));

	% Normalize adjacency matrix.
	ANorm = diag(degreeInvSqrt) * A * diag(degreeInvSqrt);

end

function [adjacency,features] = preprocessData(adjacencyData,coulombData)

    [adjacency, features] = preprocessPredictors(adjacencyData,coulombData);
    
end

function [adjacency,features] = preprocessPredictors(adjacencyData,coulombData)

    adjacency = sparse([]);
    features = [];
    
    for i = 1:size(adjacencyData, 3)
        % Extract unpadded data.
        numNodes = find(any(adjacencyData(:,:,i)),1,"last");
        
        A = adjacencyData(1:numNodes,1:numNodes,i);
        X = coulombData(1:numNodes,1:numNodes,i);
    
        % Extract feature vector from diagonal of Coulomb matrix.
        X = diag(X);
    
        % Append extracted data.
        adjacency = blkdiag(adjacency,A);
        features = [features; X];
    end

end

function [B]=MLR2(D,inter)
    
    n=length(D(:,1));
    p=length(D(1,:));
    if inter==1
        X=[ones(n,1),D(:,1:p-1)];
    elseif inter==0
        X=[D(:,1:p-1)];
    end
    Y=D(:,p);
    
    B=inv(X'*X)*X'*Y;
end


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