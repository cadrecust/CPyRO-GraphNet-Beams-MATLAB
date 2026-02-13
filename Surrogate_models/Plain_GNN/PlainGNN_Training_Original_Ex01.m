clc
clear all

%% Plain GNN - Node regression - for optimization design of RC beams
% This example shows how to train / test use plain GNNs for node regression 
% of optimum amounts of rebar in beams' cross-sections.

% Note: Training this network is a computationally intensive task. To make
% the example run quicker, this example skips the training step and loads a
% pretrained network. To instead train the network, set the doTraining 
% variable to true.
doTraining = false;

%% Load Data
% Give the path of the folder where the data is stored. 
% Adjust it to your own path.
A=importdata('Datasets/Data_1LOT_HK_Nb_Db_Simple_4000.xlsx');

DR1=A.data;

n=length(DR1(:,1));
subsize=4000;

idxSubData=ceil(rand(subsize,1)*n);

DR=DR1(idxSubData,:);

n=length(DR(:,1));
X = [DR(:,1:7)];

%% Training 

numObservations=n;

Y1 = [];
for i=1:numObservations
    As9(i,:)=[sum(DR(i,10:12).*DR(i,19:21).^2*pi/4),...
              sum(DR(i,13:15).*DR(i,22:24).^2*pi/4),...
              sum(DR(i,16:18).*DR(i,25:27).^2*pi/4)];

    if sum(As9(i,:)) > 0
        Y1=[Y1;As9(i,:)];
    end
end

n=length(Y1(:,1));
numObservations=ceil(n);

YN1=Y1(:,1); % Node 1: Left
YN2=Y1(:,2); % Node 2: Right
YN3=Y1(:,3); % Node 3: Middle

%  Change variable "targetNorm" to true if normalization of targets is 
% required. Recommended only for Bayesian Optimization of hyperparameters
targetNorm=false
if targetNorm
    muY1=mean(Y1(:,1));
    muY2=mean(Y1(:,2));
    muY3=mean(Y1(:,3));
    
    sigsqY1=var(Y1(:,1),1);
    sigsqY2=var(Y1(:,2),1);
    sigsqY3=var(Y1(:,3),1);
else
    muY1=0;
    muY2=0;
    muY3=0;
    
    sigsqY1=1;
    sigsqY2=1;
    sigsqY3=1;
end
YN1=(YN1-muY1)./sqrt(sigsqY1);
YN2=(YN2-muY2)./sqrt(sigsqY2);
YN3=(YN3-muY3)./sqrt(sigsqY3);

% Input features
X = [DR(:,1:7)];
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

%% Process input features

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

% Process Targets
nodesAb1 = [YN1(:,1),YN2(:,1),YN3(:,1)];

nodesAb1to6=[nodesAb1];

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
AsDataTrain = nodesAb1to6(idxTrain,:);
AsDataValidation = nodesAb1to6(idxValidation,:);
AsDataTest = nodesAb1to6(idxTest,:);

% Train partition

[ATrain,XTrain1,labelsTrain] = preprocessData(adjacencyDataTrain,DataXTrain1,AsDataTrain);
[~,XTrain2,~] = preprocessData(adjacencyDataTrain,DataXTrain2,AsDataTrain);
[~,XTrain3,~] = preprocessData(adjacencyDataTrain,DataXTrain3,AsDataTrain);
[~,XTrain4,~] = preprocessData(adjacencyDataTrain,DataXTrain4,AsDataTrain);
[~,XTrain5,~] = preprocessData(adjacencyDataTrain,DataXTrain5,AsDataTrain);

% Validation partition
[AValidation,XValidation1,labelsValidation] = preprocessData(adjacencyDataValidation,DataXValidation1,AsDataValidation);
[~,XValidation2,~] = preprocessData(adjacencyDataValidation,DataXValidation2,AsDataValidation);
[~,XValidation3,~] = preprocessData(adjacencyDataValidation,DataXValidation3,AsDataValidation);
[~,XValidation4,~] = preprocessData(adjacencyDataValidation,DataXValidation4,AsDataValidation);
[~,XValidation5,~] = preprocessData(adjacencyDataValidation,DataXValidation5,AsDataValidation);

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

%% Define model's architecture
parameters = struct;
numHiddenFeatureMaps1 = 16;
numHiddenFeatureMaps2 = 16;

numInputFeatures = size(XTrain,2);

sz = [numInputFeatures numHiddenFeatureMaps1];
numOut = numHiddenFeatureMaps1;
numIn = numInputFeatures;

parameters.Embedding.Weights = initializeGlorot(sz,numOut,numIn,"double");
parameters.Embedding.b = initializeZeros([1,numOut]);

% Graph attention operation
numHeads=struct;
numHeads.attn1=8;

% First attention layer
sz = [numHiddenFeatureMaps1 numHiddenFeatureMaps1];
numOut = numHiddenFeatureMaps1;
numIn = numHiddenFeatureMaps1;

parameters.attn1.Weights.linearWeights = initializeGlorot(sz,numOut,numIn,'double');
parameters.attn1.Weights.attentionWeights = initializeGlorot([numOut 2],1,2*numOut,'double');

% Multiply operation - Convolutional Layer
sz = [numHiddenFeatureMaps1 numHiddenFeatureMaps2];
numOut = numHiddenFeatureMaps2;
numIn = numHiddenFeatureMaps1;

parameters.mult1.Weights = initializeGlorot(sz,numOut,numIn,"double");

% Fully Connected Layer
numClasses = 1;
sz = [numHiddenFeatureMaps2 numClasses];
numOut = numClasses;
numIn = numHiddenFeatureMaps2;

parameters.Decoder.Weights = initializeGlorot(sz,numOut,numIn,"double");
parameters.Decoder.b = initializeZeros([1,numOut]);

%% Define training parameters

nbatches=1;

numEpochs = fix(160/nbatches);

initialLearnRate = 0.005/sqrt(nbatches);
learnRateDecay = 0.001;
validationFrequency = 5;

%% Begin training
if doTraining
    trailingAvg = [];
    trailingAvgSq = [];
    XTrain = dlarray(XTrain);
    XValidation = dlarray(XValidation);
    
    TTrain = labelsTrain;
    TValidation = labelsValidation;
    
    ngraphsTrain=size(XTrain,1)/3;
    ngraphsBatchTrain=fix(ngraphsTrain/nbatches);
    nnodesBatchTrain=ngraphsBatchTrain*3;
    
    dsTrain=cell(nbatches,2);
    
    for i=1:nbatches
        i1=(i-1)*nnodesBatchTrain+1;
        i2=(i)*nnodesBatchTrain;
    
        dsTrain{i,1}=XTrain(i1:i2,:);
        dsTrain{i,2}=TTrain(i1:i2,:);
    end
    
    monitor = trainingProgressMonitor( ...
        Metrics=["TrainingLoss","ValidationLoss"], ...
        Info="Epoch", ...
        XLabel="Epoch");
    
    groupSubPlot(monitor,"Loss",["TrainingLoss","ValidationLoss"])
    epoch = 0;
    iteration = 0;
    learningRate = initialLearnRate;
    while epoch < numEpochs && ~monitor.Stop
        epoch = epoch + 1;
        for j = 1:  nbatches
            iteration = iteration + 1;
            
            XTrain=dsTrain{j,1};
            TTrain=dsTrain{j,2};
            ATrain = sparse([]);
    
            ngraphsBatchTrain = nnodesBatchTrain/3;
            for i=1:ngraphsBatchTrain
                A=adjacencyDataTrain(:,:,i);
                ATrain = blkdiag(ATrain,A);
            end
    
            % Evaluate the model loss and gradients.
            [loss(iteration),gradients] = dlfeval(@modelLoss1fc2Gat1Conv1fc,parameters,XTrain,ATrain,TTrain,numHeads);
        
            % Update the network parameters using the Adam optimizer.
            [parameters,trailingAvg,trailingAvgSq] = adamupdate(parameters,gradients, ...
                trailingAvg,trailingAvgSq,iteration,initialLearnRate);
        
            % Record the training loss and epoch.
            recordMetrics(monitor,iteration,TrainingLoss=loss(iteration));
        
            % Display the validation metrics.
            if iteration == 1 || mod(iteration,validationFrequency) == 0
                YValidation = model1FC2GAT1Conv1FC(parameters,XValidation,AValidation,numHeads);
                
                lossValidation(iteration) = mse(YValidation,TValidation,DataFormat="BC");
                
                % Record the validation loss.
                recordMetrics(monitor,iteration,ValidationLoss=lossValidation(iteration));
            end
        end
    
        idx=randperm(nbatches);
        dsTrain=dsTrain(idx,:);
    
        updateInfo(monitor,Epoch=(epoch+" of "+numEpochs));
        learningRate = initialLearnRate / (1+learnRateDecay*epoch);
        monitor.Progress = 100*(epoch/numEpochs);
    end
    
    save('GCNN_As_Section_4000','parameters')
    save('nHeads_GAT_GCNN_As_Section_4000','numHeads')
else
    nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
    paramPIGNN=load("GCNN_As_Section_4000.mat");
end
%% Test model on unseen data
[ATest,XTest1,labelsTest] = preprocessData(adjacencyDataTest,DataXTest1,AsDataTest);
XTest1 = (XTest1 - muX1)./sqrt(sigsqX1);
XTest1 = dlarray(XTest1);

[~,XTest2,~] = preprocessData(adjacencyDataTest,DataXTest2,AsDataTest);
XTest2 = (XTest2 - muX2)./sqrt(sigsqX2);
XTest2 = dlarray(XTest2);

[~,XTest3,~] = preprocessData(adjacencyDataTest,DataXTest3,AsDataTest);
XTest3 = (XTest3 - muX3)./sqrt(sigsqX3);
XTest3 = dlarray(XTest3);

[~,XTest4,~] = preprocessData(adjacencyDataTest,DataXTest4,AsDataTest);
XTest4 = (XTest4 - muX4)./sqrt(sigsqX4);
XTest4 = dlarray(XTest4);

[~,XTest5,~] = preprocessData(adjacencyDataTest,DataXTest5,AsDataTest);
XTest5 = (XTest5 - muX5)./sqrt(sigsqX5);
XTest5 = dlarray(XTest5);

XTest=[XTest1,XTest2,XTest3,XTest4,XTest5];
if doTraining
    YTest = model1FC2GAT1Conv1FC(parameters,XTest,ATest,numHeads);
else
    YTest = model1FC2GAT1Conv1FC(paramPIGNN.parameters,XTest,ATest,nheadsparamnGATGNN.numHeads);
end

% Denormalizing targets and outputs
ntest=length(idxTest);
for i=1:ntest
    i1=(i-1)*3+1;
    i2=i*3;
    
    YTest(i1,1)=YTest(i1,1)*sqrt(sigsqY1) + muY1;
    YTest(i1+1,1)=YTest(i1+1,1)*sqrt(sigsqY2) + muY2;
    YTest(i2,1)=YTest(i2,1)*sqrt(sigsqY3) + muY3;
    
    labelsTest(i1,1)=labelsTest(i1,1)*sqrt(sigsqY1) + muY1;
    labelsTest(i1+1,1)=labelsTest(i1+1,1)*sqrt(sigsqY2) + muY2;
    labelsTest(i2,1)=labelsTest(i2,1)*sqrt(sigsqY3) + muY3;
end
%% Performance assessment

mre = [];
YAo1=dlarray([]);
Ypvt1=dlarray([]);
YAo2=dlarray([]);
Ypvt2=dlarray([]);
YAo3=dlarray([]);
Ypvt3=dlarray([]);
for i=1:ntest
    i1=(i-1)*numNodesGNN+1;
    i2=(i)*numNodesGNN;

    Ytl(i,:)=labelsTest([i1:i2]',1)';
    Yt(i,1)=sum(Ytl(i,:));
    
    Ypred(i,:)=YTest([i1:i2]',1)';
    Yp(i,1)=sum(Ypred(i,:));
    
    Aopt1=sum([Ytl(i,1)]);
    Apred1=sum(Ypred(i,1));

    Aopt2=sum([Ytl(i,2)]);
    Apred2=sum(Ypred(i,2));

    Aopt3=sum([Ytl(i,3)]);
    Apred3=sum(Ypred(i,3));

    % Mean Relative Error
    MRE1=abs(Aopt1-Apred1)/Aopt1;
    MRE2=abs(Aopt2-Apred2)/Aopt2;
    MRE3=abs(Aopt3-Apred3)/Aopt3;

    MREc = extractdata([MRE1,MRE2,MRE3]);

    % MRE
    mre=[mre; MREc'];

    % Gather predictions and labels
    YAo1=[YAo1; dlarray([Ytl(i,1)],"BC")];
    Ypvt1=[Ypvt1; dlarray([Ypred(i,1)],"BC")];

    YAo2=[YAo2; dlarray([Ytl(i,2)],"BC")];
    Ypvt2=[Ypvt2; dlarray([Ypred(i,2)],"BC")];

    YAo3=[YAo3; dlarray([Ytl(i,3)],"BC")];
    Ypvt3=[Ypvt3; dlarray([Ypred(i,3)],"BC")];
end
MRE=sum(mre)/(3*ntest);
disp('MRE')
disp(MRE)

lossTest = mse(YTest,labelsTest,DataFormat="BC");
disp('MSE')
disp(lossTest)

%% Regression coeficients
YAo1=extractdata(YAo1);
Ypvt1=extractdata(Ypvt1);
[BT1]=MLR2([[YAo1],[Ypvt1]],0);
YPT1=BT1(1).*[YAo1];

if BT1(1)>1
    BT1(1)=1/BT1(1);
end

YAo2=extractdata(YAo2);
Ypvt2=extractdata(Ypvt2);
[BT2]=MLR2([[YAo2],[Ypvt2]],0);
YPT2=BT2(1).*[YAo2];

if BT2(1)>1
    BT2(1)=1/BT2(1);
end

YAo3=extractdata(YAo3);
Ypvt3=extractdata(Ypvt3);
[BT3]=MLR2([[YAo3],[Ypvt3]],0);
YPT3=BT3(1).*[YAo3];

if BT3(1)>1
    BT3(1)=1/BT3(1);
end

disp('R coefficient - Left cross-section')
disp(BT1(1))

disp('R coefficient - Mid cross-section')
disp(BT2(1))

disp('R coefficient - Right cross-section')
disp(BT3(1))

%% Regresion plots
% Define pastel colors
pastel_gray = [0.663,0.663,0.663]; % #D3D3D3 for scatter points
pastel_green = [0.537,0.812,0.941]; % #B5EAD7 for fit line
pastel_orange = [0.992,0.761,0.506]; % #FDC281 for identity line

figure(9)
subplot(1,3,1)
plot([0;YAo1],[0;Ypvt1],'o', 'Color', pastel_gray,'LineWidth',2)
hold on
plot([0;YAo1],[0;YAo1],'-', 'Color', pastel_orange,'LineWidth',2)
hold on
plot([0;YAo1],[0;YPT1],'-', 'Color', pastel_green,'LineWidth',2)
hold on
legend(strcat('Data ','R = ',num2str(BT1(1))),'Y=T','Fit')
xlabel('$Y$',interpreter='latex')
ylabel('$\hat{Y}$',interpreter='latex')
title({strcat('True vs GNN ', ' solution:'),'Optimum rebar area of a Beam',...
    'Left section'},interpreter='latex') 
hold on
grid on
axis([[0 4000],[0 2000]])
set(gca, 'Fontname', 'Times New Roman','FontSize',12);


figure(9)
subplot(1,3,2)
plot([0;YAo2],[0;Ypvt2],'o', 'Color', pastel_gray,'LineWidth',2)
hold on
plot([0;YAo2],[0;YAo2],'-', 'Color', pastel_orange,'LineWidth',2)
hold on
plot([0;YAo2],[0;YPT2],'-', 'Color', pastel_green,'LineWidth',2)
hold on
legend(strcat('Data ','R = ',num2str(BT2(1))),'Y=T','Fit')
xlabel('$Y$',interpreter='latex')
ylabel('$\hat{Y}$',interpreter='latex')
title({strcat('True vs GNN ', ' solution:'),'Optimum rebar area of a Beam',...
    'Mid section'},interpreter='latex') 
hold on
grid on
axis([[0 3000],[0 2000]])
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

figure(9)
subplot(1,3,3)
plot([0;YAo3],[0;Ypvt3],'o', 'Color', pastel_gray,'LineWidth',2)
hold on
plot([0;YAo3],[0;YAo3],'-', 'Color', pastel_orange,'LineWidth',2)
hold on
plot([0;YAo3],[0;YPT3],'-', 'Color', pastel_green,'LineWidth',2)
hold on
legend(strcat('Data ','R = ',num2str(BT3(1))),'Y=T','Fit')
xlabel('$Y$',interpreter='latex')
ylabel('$\hat{Y}$',interpreter='latex')
title({strcat('True vs GNN ', ' solution:'),'Optimum rebar area of a Beam',...
    'Right section'},interpreter='latex') 
hold on
grid on
axis([[0 4000],[0 2000]])
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

%% Function appendix

function [adjacency,features,labels] = preprocessData(adjacencyData,coulombData,atomData)

    [adjacency, features] = preprocessPredictors(adjacencyData,coulombData);
    labels = [];
    
    % Convert labels to categorical.
    for i = 1:size(adjacencyData,3)
        % Extract and append unpadded data.
        T = atomData(i,:);
        labels = [labels; T'];
    end
    
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

function Y = model1FC2GAT1Conv1FC(parameters,X,A,numHeads)

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

function [loss,gradients] = modelLoss1fc2Gat1Conv1fc(parameters,X,A,T,numHeads)
    
    Y = model1FC2GAT1Conv1FC(parameters,X,A,numHeads);
    
    loss = mse(Y,T,DataFormat="BC");
    gradients = dlgradient(loss, parameters);

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