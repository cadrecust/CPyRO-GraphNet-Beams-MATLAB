% Optimal_Design_SSFSRecBeams_Complete_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To test the performance of the PSO assisted by CPyRO-GraphNet-Beams
%    and/or plain GNN, for the optimum design of rebar in a single span
%    beam.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2025-07-03
% Copyright (c)  School of Engineering
%                HKUST
%----------------------------------------------------------------

clc
clear all

%% Geometry 
span=4000; % mm
b=350; % width (mm)
h=650; % height (mm)
        
hrec=50; % 
brec=50; % lateral concrete cover

%% Materials
fcu=25; % Kg/mm2

fy=500; % Yield stress of steel reinforcement (N/mm2)
wac=78.5e-6; % unit volume weight of the reinforcing steel (N/mm3)

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

W=[20.8,22.2]; % Uniformly Distributed Load over the whole beam  N / mm
nspans=length(W(:,1));

[R,U,V,M]=MSFSFEMBeams(span,Ac,Ic,Ec,supportsLoc,W,dL,wrange,0);

ne=zeros(1,nspans);
neSum=0;
load_conditions=[];
cutLoc=[];
for i=1:nspans
    
    ne(i)=(supportsLoc(i+1)-supportsLoc(i))/dL;
    i1=neSum+1;
    Mleft=M(1,i1);
    Lv(i,1)=span(i);
    bv(i,1)=b(i);
    hv(i,1)=h(i);
    fcuv(i,1)=fcu(i);

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
model2Use="PIGCNN";

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

%% OPTIMAL DESIGN 

%% Generalization parameters
Nv=[nBeams];
ndl=length(Nv);

N=Nv(1);

[idxTrain,idxValidation,idxTest] = trainingPartitions(N,[0 0 1]);

nTest=length(idxTest);

% Data for GNNs
meanX=[349.605734767025	633.225806451613	37.4534050179212	4505.55555555556	-26339479.7363044];
sigsqX=[4949.66534345653	11216.8343161059	31.3195134954590	996295.300677007	2.10064997682719e+15];

[XTrain,XValidation,XTest,ATrain,AValidation,ATest]=dataTrainNLClass...
    (DR,idxTrain,idxValidation,idxTest,meanX,sigsqX);

XTestOpt = DR(idxTest,1:7);

%% Load surrogate models
% Model for prediction of Number of Layers

% Model for prediction of As per cross-section
if model2Use=="PIGCNN"
    nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_4000.mat");
    paramPIGCNN=load("PIGCNN_As_Section_4000.mat");

    Ao3=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,ATest,nheadsparamnGATPIGNN.numHeads);
elseif model2Use=="GCNN"
    nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
    paramGCNN=load("GCNN_As_Section_4000.mat");

    Ao3=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,ATest,nheadsparamnGATGNN.numHeads);
end

Ao3=extractdata(Ao3);

fcu=XTestOpt(:,3);
span=XTestOpt(:,4);

b=XTestOpt(:,1);
h=XTestOpt(:,2);

%% Loads
Mleft=XTestOpt(:,5);
Mmid=XTestOpt(:,6);
Mright=XTestOpt(:,7);

hagg=20;

nfig=1;

[volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db2Spans,EffSpans,...
MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
CFAspans]=SASOOptimMSFSBeamsRebarBasic(b,h,span,brec,hrec,hagg,...
pmin,pmax,rebarAvailable([2:10]',:),fcu,load_conditions,fy,wac,cutxLocBeam,...
WfaConstr,Ao3,1,1,[1]);

for i=1:nspans
    % Average percentage of cross-section reinforcement 
    rho=(sum(tenbLMRspan(i,1:3).*(db2Spans(i,1).^2*pi/4))+...
        sum(tenbLMRspan(i,4:6).*(db2Spans(i,2).^2*pi/4))+...
        sum(tenbLMRspan(i,7:9).*(db2Spans(i,3).^2*pi/4)))./(b*h)/3;  
    
    [s1(i,1),s2(i,1),s3(i,1),d1(i,1),d2(i,1)]=shearDesignBeams(span,b,h,...
                                                    hrec,fcu,fy,V,dvs,rho);
    
    dvsBeams=dvs(i,1);
end

%% Side rebars (if necessary)
if h>=750
    [dSb,nSb,sepSb,distrSideBars]=sideBarsRecBeams3SecSpan(b,h,fy,...
            brec,hrec,tenbLMRspan,db18Spans,dvs,hagg,rebarAvailable);
    
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
%{
directionData='C:\Users\lfver\OneDrive - HKUST Connect\PhD\Rebar_Design_Material\';

ExportDesignSSRecBeam(directionData,[b,h,span,brec,hrec],[0,0,0],LenRebarL,...
    LenRebarM,LenRebarR,-DistrRebarLeft,DistrRebarMid,-DistrRebarRight,totnbSpan,...
    tenbLMRspan,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
    diamlistdSb,distrSideBars,beamNSb,[s1,s2,s3,d1,d2,dvsBeams]);
%}
%% Function appendix
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
sigsqX=[sigsqX1,sigsqX2,sigsqX3,sigsqX4,sigsqX5];
meanX=[muX1,muX2,muX3,muX4,muX5];

%% Normalizing validation data
XValidation1 = (XValidation1 - muX1)./sqrt(sigsqX1);
XValidation2 = (XValidation2 - muX2)./sqrt(sigsqX2);

XValidation3 = (XValidation3 - muX3)./sqrt(sigsqX3);
XValidation4 = (XValidation4 - muX4)./sqrt(sigsqX4);
XValidation5 = (XValidation5 - muX5)./sqrt(sigsqX5);

XValidation=[XValidation1,XValidation2,XValidation3,XValidation4,XValidation5];

%% Normalizing test data
muX1=meanX(1);
sigsqX1=sigsqX(1);

muX2=meanX(2);
sigsqX2=sigsqX(2);

muX3=meanX(3);
sigsqX3=sigsqX(3);

muX4=meanX(4);
sigsqX4=sigsqX(4);

muX5=meanX(5);
sigsqX5=sigsqX(5);

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

