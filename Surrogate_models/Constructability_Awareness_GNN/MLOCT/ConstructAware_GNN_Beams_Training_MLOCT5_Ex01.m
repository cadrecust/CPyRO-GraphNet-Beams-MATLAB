
%% Plain GNN - Node classification - for optimization design of RC beams
% This example shows how to classify nodes to estimate optimum amounts of
% number of rebar layers for beams' cross-sections.

% The network is trained specifically for optimum rebar designs with the
% largest amount of rebar weight (MLOCT5), according to a database where
% multi-objective optimizaiton was used, and from which 5 different optimum 
% levels along the Pareto Front were extracted.
clc
clear all

%% Load Data
% Give the path of the folder where the data is stored. 
% Adjust it to your own path.
A=importdata('Enhanced_Data_MOO/Enhanced_Data_5LOT_HK_Nb_Db_Simple_4000.xlsx');

DR=A.data(1:4000,:);

%% Training 

doTraining = true;

numObservations=ceil(n);
Y = [DR(1:numObservations,46:54)];

Y1 = [];
for i=1:numObservations
    if sum(Y(i,:))~=0
        nlay1=sum(Y(i,1:3)>0);
        nlay2=sum(Y(i,4:6)>0);
        nlay3=sum(Y(i,7:9)>0);
        
        Y1=[Y1;[nlay1,nlay2,nlay3]];
    end
end

n=length(Y1(:,1));
numObservations=ceil(n);
trainDataN=numObservations;

YN1=Y1(:,1); % Node 1: Left
YN2=Y1(:,2); % Node 2: Middle
YN3=Y1(:,3); % Node 3: Right

%% GNN architecture
elements=[2 1;
          1 3];

%% Adjacency matrix
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

%% Features

features1=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,5)];
features2=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,6)];
features3=[X(:,1),X(:,2),X(:,3),X(:,4),X(:,7)];

DataX1=zeros(numObservations,3,3);
DataX2=zeros(numObservations,3,3);
DataX3=zeros(numObservations,3,3);
DataX4=zeros(numObservations,3,3);
DataX5=zeros(numObservations,3,3);
for i=1:numObservations
    features=[features1(i,:)',features2(i,:)',features3(i,:)'];

    for j=1:numNodesGNN
        DataX1(i,j,j)=features(1,j);
        DataX2(i,j,j)=features(2,j);
        DataX3(i,j,j)=features(3,j);
        DataX4(i,j,j)=features(4,j);
        DataX5(i,j,j)=features(5,j);
    end
end

DataX1 = double(permute(DataX1, [2 3 1]));
DataX2 = double(permute(DataX2, [2 3 1]));
DataX3 = double(permute(DataX3, [2 3 1]));
DataX4 = double(permute(DataX4, [2 3 1]));
DataX5 = double(permute(DataX5, [2 3 1]));

nodesNlay3 = [YN1(:,1),YN2(:,1),YN3(:,1)];

nodesNlay1to3=[nodesNlay3];

%% Categories counting (N layers on all cross section)
[symbol,count] = NRebarLayers(nodesNlay1to3);
figure(4)
histogram(categorical(symbol))
xlabel("Node Label - N layers 1 to 3")
ylabel("Frequency")
title(strcat("Label Counts: Optimum Constructability Target 5. N=",num2str(subsize)))
set(gca, 'Fontname', 'Times New Roman','FontSize',20);

%% Adjacency
adjacency=repmat(adjacency,[1,1,n]);
size(adjacency)

%% Partition of data
numObservations = size(adjacency,3);
[idxTrain,idxValidation,idxTest] = trainingPartitions(numObservations,[0.7 0.15 0.15]);

% node adjacency data
adjacencyDataTrain = adjacency(:,:,idxTrain);
adjacencyDataValidation = adjacency(:,:,idxValidation);
adjacencyDataTest = adjacency(:,:,idxTest);

% feature data
DataXTrain1 = DataX1(:,:,idxTrain);
DataXValidation1 = DataX1(:,:,idxValidation);
DataXTest1 = DataX1(:,:,idxTest);

DataXTrain2 = DataX2(:,:,idxTrain);
DataXValidation2 = DataX2(:,:,idxValidation);
DataXTest2 = DataX2(:,:,idxTest);

DataXTrain3 = DataX3(:,:,idxTrain);
DataXValidation3 = DataX3(:,:,idxValidation);
DataXTest3 = DataX3(:,:,idxTest);

DataXTrain4 = DataX4(:,:,idxTrain);
DataXValidation4 = DataX4(:,:,idxValidation);
DataXTest4 = DataX4(:,:,idxTest);

DataXTrain5 = DataX5(:,:,idxTrain);
DataXValidation5 = DataX5(:,:,idxValidation);
DataXTest5 = DataX5(:,:,idxTest);

% target data
diamDataTrain = nodesNlay1to3(idxTrain,:);
diamDataValidation = nodesNlay1to3(idxValidation,:);
diamDataTest = nodesNlay1to3(idxTest,:);

% Train partition

[ATrain,XTrain1,labelsTrain] = preprocessData(adjacencyDataTrain,DataXTrain1,diamDataTrain);
[~,XTrain2,~] = preprocessData(adjacencyDataTrain,DataXTrain2,diamDataTrain);
[~,XTrain3,~] = preprocessData(adjacencyDataTrain,DataXTrain3,diamDataTrain);
[~,XTrain4,~] = preprocessData(adjacencyDataTrain,DataXTrain4,diamDataTrain);
[~,XTrain5,~] = preprocessData(adjacencyDataTrain,DataXTrain5,diamDataTrain);

% Validation partition
[AValidation,XValidation1,labelsValidation] = preprocessData(adjacencyDataValidation,DataXValidation1,diamDataValidation);
[~,XValidation2,~] = preprocessData(adjacencyDataValidation,DataXValidation2,diamDataValidation);
[~,XValidation3,~] = preprocessData(adjacencyDataValidation,DataXValidation3,diamDataValidation);
[~,XValidation4,~] = preprocessData(adjacencyDataValidation,DataXValidation4,diamDataValidation);
[~,XValidation5,~] = preprocessData(adjacencyDataValidation,DataXValidation5,diamDataValidation);

%% Normalizing training data
muX1 = mean(XTrain1);
sigsqX1 = var(XTrain1,1);
XTrain1 = (XTrain1 - muX1)./sqrt(sigsqX1);
XValidation1 = (XValidation1 - muX1)./sqrt(sigsqX1);

muX2 = mean(XTrain2);
sigsqX2 = var(XTrain2,1);
XTrain2 = (XTrain2 - muX2)./sqrt(sigsqX2);
XValidation2 = (XValidation2 - muX2)./sqrt(sigsqX2);

muX3 = mean(XTrain3);
sigsqX3 = var(XTrain3,1);
XTrain3 = (XTrain3 - muX3)./sqrt(sigsqX3);
XValidation3 = (XValidation3 - muX3)./sqrt(sigsqX3);

muX4 = mean(XTrain4);
sigsqX4 = var(XTrain4,1);
XTrain4 = (XTrain4 - muX4)./sqrt(sigsqX4);
XValidation4 = (XValidation4 - muX4)./sqrt(sigsqX4);

muX5 = mean(XTrain5);
sigsqX5 = var(XTrain5,1);
XTrain5 = (XTrain5 - muX5)./sqrt(sigsqX5);
XValidation5 = (XValidation5 - muX5)./sqrt(sigsqX5);

XTrain=[XTrain1,XTrain2,XTrain3,XTrain4,XTrain5];
XValidation=[XValidation1,XValidation2,XValidation3,XValidation4,XValidation5];

%% DL model

parameters = struct;
numHiddenFeatureMaps = 32;
numHiddenFeatureMaps2 = 16;
numInputFeatures = size(XTrain,2);

%% GAT layers
% Graph attention operation
numHeads=struct;

numHeads.attn1=4;
numHeads.attn2=2;

% First attention layer
sz = [numInputFeatures numHiddenFeatureMaps];
numOut = numHiddenFeatureMaps;
numIn = numInputFeatures;

parameters.attn1.Weights.linearWeights = initializeGlorot(sz,numOut,numIn,'double');
parameters.attn1.Weights.attentionWeights = initializeGlorot([numOut 2],1,2*numOut,'double');

% Second layer
sz = [numHiddenFeatureMaps numHiddenFeatureMaps];
numOut = numHiddenFeatureMaps;
numIn = numHiddenFeatureMaps;

parameters.attn2.weights.linearWeights = initializeGlorot(sz,numOut,numIn);
parameters.attn2.weights.attentionWeights = initializeGlorot([numOut 2],1,2*numOut);

% Multiply operation

classes = categories(labelsTrain);
numClasses = numel(classes);

sz = [numHiddenFeatureMaps numHiddenFeatureMaps2];
numOut = numHiddenFeatureMaps2;
numIn = numHiddenFeatureMaps;

parameters.mult.Weights = initializeGlorot(sz,numOut,numIn,"double");

%% Define training parameters
numEpochs = 120;
learnRate = 0.02;

validationFrequency = 10;
%% Begin training
if doTraining
    trailingAvg = [];
    trailingAvgSq = [];
    XTrain = dlarray(XTrain);
    XValidation = dlarray(XValidation);
    
    if canUseGPU
        XTrain = gpuArray(XTrain);
    end
    TTrain = onehotencode(labelsTrain,2,ClassNames=classes);
    TValidation = onehotencode(labelsValidation,2,ClassNames=classes);
    
    
    monitor = trainingProgressMonitor( ...
        Metrics=["TrainingLoss","ValidationLoss"], ...
        Info="Epoch", ...
        XLabel="Epoch");
    
    groupSubPlot(monitor,"Loss",["TrainingLoss","ValidationLoss"])
    epoch = 0;
    
    while epoch < numEpochs && ~monitor.Stop
        epoch = epoch + 1;
    
        % Evaluate the model loss and gradients.
        [loss,gradients] = dlfeval(@modelLoss2GAT1Conv,parameters,XTrain,ATrain,TTrain,numHeads);
    
        % Update the network parameters using the Adam optimizer.
        [parameters,trailingAvg,trailingAvgSq] = adamupdate(parameters,gradients, ...
            trailingAvg,trailingAvgSq,epoch,learnRate);
    
        % Record the training loss and epoch.
        recordMetrics(monitor,epoch,TrainingLoss=loss);
        updateInfo(monitor,Epoch=(epoch+" of "+numEpochs));
    
        % Display the validation metrics.
        if epoch == 1 || mod(epoch,validationFrequency) == 0
            YValidation = model2GAT1Conv(parameters,XValidation,AValidation,numHeads);
            lossValidation = crossentropy(YValidation,TValidation,DataFormat="BC");
    
            % Record the validation loss.
            recordMetrics(monitor,epoch,ValidationLoss=lossValidation);
        end
    
        monitor.Progress = 100*(epoch/numEpochs);
    end

    save('nLay_GNN__MOConstrucT5_4000.mat', 'parameters');
    save('nHead_nLay_GNN_MOConstrucT5_4000.mat', 'numHeads');
else

    paramNL=load('nLay_GNN_MOConstrucT5_4000.mat');
    nheadsparamnNLGNN=load('nHead_nLay_GNN_MOConstrucT5_4000.mat');
end
%% Test model
[ATest,XTest1,labelsTest] = preprocessData(adjacencyDataTest,DataXTest1,diamDataTest);
XTest1 = (XTest1 - muX1)./sqrt(sigsqX1);
XTest1 = dlarray(XTest1);

[~,XTest2,~] = preprocessData(adjacencyDataTest,DataXTest2,diamDataTest);
XTest2 = (XTest2 - muX2)./sqrt(sigsqX2);
XTest2 = dlarray(XTest2);

[~,XTest3,~] = preprocessData(adjacencyDataTest,DataXTest3,diamDataTest);
XTest3 = (XTest3 - muX3)./sqrt(sigsqX3);
XTest3 = dlarray(XTest3);

[~,XTest4,~] = preprocessData(adjacencyDataTest,DataXTest4,diamDataTest);
XTest4 = (XTest4 - muX4)./sqrt(sigsqX4);
XTest4 = dlarray(XTest4);

[~,XTest5,~] = preprocessData(adjacencyDataTest,DataXTest5,diamDataTest);
XTest5 = (XTest5 - muX5)./sqrt(sigsqX5);
XTest5 = dlarray(XTest5);

XTest=[XTest1,XTest2,XTest3,XTest4,XTest5];
if doTraining
    YTest = model2GAT1Conv(parameters,XTest,ATest,numHeads);
else
    YTest = model2GAT1Conv(paramNL.parameters,XTest,ATest,nheadsparamnNLGNN.numHeads);
end

YTest = onehotdecode(YTest,classes,2);

nTest=size(idxTest,2);
for i=1:nTest
    i2=i*3;
    if double(string(YTest(i2,1)))==3 
        YTest(i2,1)=categorical(2);
    end
end
TTest = onehotencode(labelsTest,2,ClassNames=classes);

accuracy = mean(YTest == labelsTest)

%% Confusion matrix

figure
cm = confusionchart(labelsTest,YTest, ...
    ColumnSummary="column-normalized", ...
    RowSummary="row-normalized");
title(strcat("GCN-NLayer Confusion Chart: Optimum Constructability Target 5. N=",num2str(subsize)))
set(gca, 'Fontname', 'Times New Roman','FontSize',20);

%% Function appendix
function [loss,gradients] = modelLoss2GAT1Conv(parameters,X,A,T,numHeads)
    
    Y = model2GAT1Conv(parameters,X,A,numHeads);
    loss = crossentropy(Y,T,DataFormat="BC");
    gradients = dlgradient(loss, parameters);

end

function Y = model2GAT1Conv(parameters,X,A,numHeads)

    weights1 = parameters.attn1.Weights;
    numHeadsAttention1 = numHeads.attn1;
    
    Z1 = X;
    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = elu(Z2);

    weights2 = parameters.attn2.weights;
    numHeadsAttention2 = numHeads.attn2;
    
    [Z3,~] = graphAttention(Z2,A,weights2,numHeadsAttention2,"cat");
    Z3  = elu(Z3) + Z2;
    
    ANorm = normalizeAdjacency(A);
    Z4 = single(full(ANorm)) * Z3 * double(parameters.mult.Weights);
    
    Z4 = graphPoolingLayer(Z4, 'mean', 3, 3);

    Y = softmax(Z4,DataFormat="BC");

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

function [adjacency,features,labels] = preprocessData(adjacencyData,coulombData,atomData)

    [adjacency, features] = preprocessPredictors(adjacencyData,coulombData);
    labels = [];
    
    % Convert labels to categorical.
    for i = 1:size(adjacencyData,3)
        % Extract and append unpadded data.
        T = nonzeros(atomData(i,:));
        labels = [labels; T];
    end
    
    labels2 = nonzeros(atomData);
    %labels2 = atomData;
    assert(isequal(labels2,labels2))
    
    diamSizeNodes = unique(labels);

    diamSizes =  NRebarLayers(diamSizeNodes);
    labels = categorical(labels, diamSizeNodes, diamSizes);

end


function [symbol,count] = NRebarLayers(atomicNum)
% ATOMICSYMBOL Convert atomic number to symbol
%   symbol = atomicSymbol(atomicNum) returns the atomic symbol of the
%   specified atomic number.
%
%   [symbol,count] = atomicSymbol(atomicNum) also returns the count for
%   each symbol.
%

numSymbols = numel(atomicNum);

symbol = strings(numSymbols, 1);
count = strings(numSymbols,1);


Count1 = 0;
Count2 = 0;
Count3 = 0;
for i = 1:numSymbols
    switch atomicNum(i)
        case 1
            symbol(i) = 1;
            Count1 = Count1 + 1;
            count(i) = "1" + Count1;
        case 2
            symbol(i) = 2;
            Count2 = Count2 + 1;
            count(i) = "2" + Count2;
        case 3
            symbol(i) = 3;
            Count3 = Count3 + 1;
            count(i) = "3" + Count3;
    end
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

function y = elu(x)

    y = max(0, x) + (exp(min(0, x)) -1);

end

function [outputFeatures,normAttentionCoeff] = graphAttention(inputFeatures,adjacency,weights,numHeads,aggregation)
    
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

function pooledFeatures = graphPoolingLayer(nodeFeatures, poolType, outputDim, numNodes)

    % Input validation
    if ~ischar(poolType) || ~ismember(lower(poolType), {'mean', 'max', 'sum'})
        error('poolType must be ''mean'', ''max'', or ''sum''.');
    end
    numFeatures = size(nodeFeatures, 2); % Number of input features per node, e.g., 5
    nGraphs = size(nodeFeatures, 1)/numNodes;
    % If outputDim equals input feature dimension, apply pooling without coarsening
    if outputDim == numFeatures
        switch lower(poolType)
            case 'mean'
                for i=1:nGraphs
                    i1=(i-1)*numNodes+1;
                    i2=i*numNodes;
                    pooledFeatures(i1:i2,:) = mean(nodeFeatures(i1:i2,:), 2) * ones(1, outputDim,'like',nodeFeatures);
                end
            case 'max'
                for i=1:nGraphs
                    i1=(i-1)*numNodes +1;
                    i2=i*numNodes;
                    pooledFeatures(i1:i2,:) = max(nodeFeatures(i1:i2,:), [], 2) * ones(1, outputDim,'like',nodeFeatures);
                end
            case 'sum'
                for i=1:nGraphs
                    i1=(i-1)*numNodes +1;
                    i2=i*numNodes;
                    pooledFeatures(i1:i2,:) = sum(nodeFeatures(i1:i2,:), 2) * ones(1, outputDim,'like',nodeFeatures);
                end
        end
        return;
    end

    % Determine grouping of features for pooling
    if outputDim < numFeatures
        % Group features into approximately equal-sized subsets
        groupSize = ceil(numFeatures / outputDim);
        pooledFeatures = dlarray(zeros(numNodes*nGraphs, outputDim));
        
        for i = 1:outputDim
            startIdx = (i-1) * groupSize + 1;
            endIdx = min(i * groupSize, numFeatures);
            if startIdx > numFeatures
                break; % No more features to process
            end
            featureSubset = nodeFeatures(:, startIdx:endIdx);
            
            % Apply pooling operation
            switch lower(poolType)
                case 'mean'
                    for j=1:nGraphs
                        i1=(j-1)*numNodes+1;
                        i2=j*numNodes;
                        pooledFeatures(i1:i2, i) = mean(featureSubset(i1:i2,:), 2);
                    end
                case 'max'
                    for j=1:nGraphs
                        i1=(j-1)*numNodes +1;
                        i2=j*numNodes;
                        pooledFeatures(i1:i2, i) = max(featureSubset(i1:i2,:), [], 2);
                    end
                case 'sum'
                    for j=1:nGraphs
                        i1=(j-1)*numNodes +1;
                        i2=j*numNodes;
                        pooledFeatures(i1:i2, i) = sum(featureSubset(i1:i2,:), 2);
                    end
            end
        end
    else
        % If outputDim > numFeatures, pad with zeros or repeat features
        pooledFeatures = zeros(numNodes*nGraphs, outputDim);
        switch lower(poolType)
            case 'mean'
                for j=1:nGraphs
                    i1=(j-1)*numNodes+1;
                    i2=j*numNodes;
                    baseFeatures(i1:i2,:) = mean(nodeFeatures(i1:i2,:), 2) * ones(1, numFeatures,'like',nodeFeatures);
                end
            case 'max'
                for j=1:nGraphs
                    i1=(j-1)*numNodes+1;
                    i2=j*numNodes;
                    baseFeatures(i1:i2,:) = max(nodeFeatures(i1:i2,:), [], 2) * ones(1, numFeatures,'like',nodeFeatures);
                end
            case 'sum'
                for j=1:nGraphs
                    i1=(j-1)*numNodes+1;
                    i2=j*numNodes;
                    baseFeatures(i1:i2,:) = sum(nodeFeatures(i1:i2,:), 2) * ones(1, numFeatures,'like',nodeFeatures);
                end
        end
        % Fill output dimensions with base features and pad with zeros if necessary
        pooledFeatures(:, 1:numFeatures) = baseFeatures;
    end
end

function weights = initializeGlorot(sz,numOut,numIn,className)

    arguments
        sz
        numOut
        numIn
        className = 'single'
    end
    
    Z = 2*rand(sz,className) - 1;
    bound = sqrt(6 / (numIn + numOut));
    
    weights = bound * Z;
    weights = dlarray(weights);

end

function parameter = initializeZeros(sz,className)
    
    arguments
        sz
        className = 'single'
    end
    
    parameter = zeros(sz,className);
    parameter = dlarray(parameter);

end