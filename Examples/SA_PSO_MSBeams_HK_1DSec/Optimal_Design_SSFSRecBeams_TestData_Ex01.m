% Optimal_Design_SSFSRecBeams_TestData_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To test the performance of the PSO assisted by CPyRO-GraphNet-Beams
%    and/or plain GNN for the optimu design of rebar in single-span beams.
%    A maximum of 1 rebar diameter size is considered for each beam
%    cross-section.
%
%----------------------------------------------------------------
%
% LAST MODIFIED: L.F.Veduzco    2025-07-03
% Copyright (c)  School of Engineering
%                HKUST
%----------------------------------------------------------------

clc
clear all

A=importdata('Enhanced_Data_1LOT_HK_Nb_Db_Simple_4000.xlsx');

DR2=A.data;

nObservations=size(DR2,1);

subsize=2000;
idxSubData=ceil(rand(subsize,1)*nObservations);
DR=DR2(idxSubData,:);
nObservations=length(DR(:,1));

i=0;
dL=100;
for j =1:nObservations
    As9(j,:)=[sum(DR2(j,10:12).*DR2(j,19:21).^2*pi/4),...
              sum(DR2(j,13:15).*DR2(j,22:24).^2*pi/4),...
              sum(DR2(j,16:18).*DR2(j,25:27).^2*pi/4)];

    A13=As9(j,:)';
    if sum([A13']) > 0
        i = i + 1;
        DR(i,:)=DR2(j,:);
        fcu(1,i)=DR2(j,3);
        b(1,i)=DR2(j,1);
        h(1,i)=DR2(j,2);
        L(1,i)=DR2(j,4);
        Mul(1,i)=abs(DR2(j,5));
        Mum(1,i)=DR2(j,6);
        Mur(1,i)=abs(DR2(j,7));
        Wleft=DR2(j,8);
        Wright=DR2(j,9);

        Ac = b(1,i) * h(1,i);
        Ic=b(1,i)*h(1,i)^3/12;
        Ec=(3.46*sqrt(fcu(1,i))+3.21)*1e3;

        [R,U,V,M]=MSFSFEMBeams(L(1,i),Ac,Ic,Ec,[0,L(1,i)],[Wleft,Wright],dL,[0,L(1,i)],0);
        [Mmid,mp]=max(M(1,:));
        xMmid = (mp-1) * dL ;

        %% LOcation of cuts per beam
        cutxLocBeam(i,:)=cutLocationSSRecBeam(M(:,1:end),dL);
        
        %% Rebar cross-section area quantities

        u1L(:,i)=A13;

        %% Location of design moments Mu
        x0L(:,i)=[10,xMmid,L(1,i)-10]';
    end
end
nObservations=size(x0L,2);

model2Use="PIGCNN";

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
                9 40
                10 50]; % mm^2

dvs=10;
pmin=0.003;
pmax=0.025;

%% Constructability
hagg=20;
Wunb=[1.3,1.4];
Wnd=[0.6];
Wcut=[1.3,1.6];
Wnb=[0.8];
Wcs1=[0.6];
Wcs2=[2];
Wfac=[Wunb,Wnd,Wcut,Wnb,Wcs1,Wcs2];

%% Generalization parameters

MREv=[];
MSEv=[];
R2v=[];
MAEv=[];
accuracyBeaml=[];
Rcoefv=[];
nreps=1;

Nv=[nObservations];
ndl=length(Nv);

for nd=1:ndl  % loop on vector of data size
    
    nodesAb1to3=u1L';
    N=Nv(nd);
    
    [idxTrain,idxValidation,idxTest] = trainingPartitions(N,[0.7 0.1 0.2]);
    nTrain=length(idxTrain);
    nValidation=length(idxValidation);
    nTest=length(idxTest);
    
    % Target data
    AsDataTrain = nodesAb1to3(idxTrain,:);
    AsDataValidation = nodesAb1to3(idxValidation,:);
    AsDataTest = nodesAb1to3(idxTest,:);

    cutxLocBeamTrain=cutxLocBeam(idxTrain,:);
    cutxLocBeamValidation=cutxLocBeam(idxValidation,:);
    cutxLocBeamTest=cutxLocBeam(idxTest,:);
    
    % Data for GNNs
    meanX=[349.605734767025	633.225806451613	37.4534050179212	4505.55555555556	-26339479.7363044];
    sigsqX=[4949.66534345653	11216.8343161059	31.3195134954590	996295.300677007	2.10064997682719e+15];
    
    [XTrain,XValidation,XTest,ATrain,AValidation,ATest]=dataTrainNLClass...
        (DR,idxTrain,idxValidation,idxTest,meanX,sigsqX);

    YTest = AsDataTest;
    XTestOpt = DR(idxTest,1:7);

    %% Load surrogate models
    % Model for prediction of Number of Layers

    % Model for prediction of As per cross-section
    if model2Use=="PIGCNN"

        nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_500.mat");
        paramPIGCNN=load("PIGCNN_As_Section_500.mat");

        Ao3=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,ATest,nheadsparamnGATPIGNN.numHeads);
    else
        nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
        paramGCNN=load("GCNN_As_Section_4000.mat");
        
        Ao3=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,ATest,nheadsparamnGATGNN.numHeads);
    end
    
    Ao3=extractdata(Ao3);
    
    for reps=1:nreps % loop on data
        
        Y_test_At=zeros(nTest,1);
        for j=1:nTest
            Y_test_At(j,1)=mean(YTest(j,:));
        end
        s1=sum(Y_test_At)/nTest;
        
        SST=0;
        for j=1:nTest
            SST=SST+(Y_test_At(j,1)-s1)^2;
        end
        
        Ast3=zeros(nTest,3);

        Aopt3=zeros(nTest,3);
        Aeff_Well=[];
        Aeff_Poor=[];
        SSres=0;
        MAE=0;
        MSE=0;
        Y_test_Aoriginal=[];
        Ypvt=[];
        YAo1=[];
        YAo2=[];
        YAo3=[];
        Ypvt1=[];
        Ypvt2=[];
        Ypvt3=[];
        mre=[];
        n1=1;
        n2=nTest;
        
        for i=n1:n2 % loop on samples
            
            fcu=XTestOpt(i,3);
            b=XTestOpt(i,1);
            h=XTestOpt(i,2);
            span=XTestOpt(i,4);
            
            %% Loads
            Mleft=XTestOpt(i,5);
            Mmid=XTestOpt(i,6);
            Mright=XTestOpt(i,7);

            load_conditions=[1 Mleft Mmid Mright]; %Kg-cm (flexure)
            
            %% Cut location ( local coordinates)
            cutxLoc=cutxLocBeamTest(i,:);

            %% OPTIMAL DESIGN 
            i1=(i-1)*3+1;
            i2=(i)*3;
            
            Aos3=Ao3(i1:i2,1);
            
            [volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db2Spans,EffSpans,...
            MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
            DistrRebarLeft,DistrRebarMid,DistrRebarRight,tenbLMRspan,totnbSpan,...
            CFAspans]=SASOOptimMSFSBeamsRebarBasic(b,h,span,brec,hrec,hagg,...
            pmin,pmax,rebarAvailable([2:10]',:),fcu,load_conditions,fy,wac,cutxLoc,...
            Wfac,Aos3,0,0,[]);
            
            bestAbi2=db2Spans.^2*pi/4;
            bestAbTL1=sum(tenbLMRspan(1,1:3).*bestAbi2(1));
            bestAbBM1=sum(tenbLMRspan(1,4:6).*bestAbi2(2));
            bestAbTR1=sum(tenbLMRspan(1,7:9).*bestAbi2(3));
            
            Ast3(i,:)=[bestAbTL1,bestAbBM1,bestAbTR1];
            
            %Ast3(i,:)=Aos3';
            Aopt3(i,:)=YTest(i,:);
            
            % Mean Relative Error
            A1opt=Aopt3(i,1);
            A2opt=Aopt3(i,2);
            A3opt=Aopt3(i,3);
            
            Apred1=Ast3(i,1);
            Apred2=Ast3(i,2);
            Apred3=Ast3(i,3);
            
            MRE1=abs(A1opt-Apred1)/A1opt;
            MRE2=abs(A2opt-Apred2)/A2opt;
            MRE3=abs(A3opt-Apred3)/A3opt;

            MREC = [MRE1,MRE2,MRE3];

            mre=[mre; MREC'];
            
            % Accuracy
            for j=1:3
                if MREC(j)<=0.15
                    Aeff_Well=Aeff_Well+1;
                else
                    Aeff_Poor=Aeff_Poor+1;
                end
            end
            Ap=mean(Ast3(i,:));
            Ao=mean(Aopt3(i,:));
            
            % R coefficient
            YAo1=[YAo1; A1opt];
            Ypvt1=[Ypvt1; Apred1];

            YAo2=[YAo2; A2opt];
            Ypvt2=[Ypvt2; Apred2];

            YAo3=[YAo3; A3opt];
            Ypvt3=[Ypvt3; Apred3];
            
        end
        
        [BT1]=MLR2([[YAo1],[Ypvt1]],0);
        YPT1=BT1(1).*[YAo1];

        if BT1(1)>1
            BT1(1)=1/BT1(1);
        end

        [BT2]=MLR2([[YAo2],[Ypvt2]],0);
        YPT2=BT2(1).*[YAo2];

        if BT2(1)>1
            BT2(1)=1/BT2(1);
        end

        [BT3]=MLR2([[YAo3],[Ypvt3]],0);
        YPT3=BT3(1).*[YAo3];

        if BT3(1)>1
            BT3(1)=1/BT3(1);
        end

        disp('R coefficient')
        disp(BT1(1))

        disp('R coefficient')
        disp(BT2(1))

        disp('R coefficient')
        disp(BT3(1))
        
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
        title({strcat('True vs ',model2Use, '-PSO solution (1DSec RP):'),'Optimum rebar area of a Beam',...
            'Left section'},interpreter='latex') 
        hold on
        grid on
        axis([0 3000 0 2500])
        set(gca, 'Fontname', 'Times New Roman','FontSize',18);


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
        title({strcat('True vs ',model2Use, '-PSO solution (1DSec RP):'),'Optimum rebar area of a Beam',...
            'Mid section'},interpreter='latex') 
        hold on
        grid on
        axis([0 3000 0 3500])
        set(gca, 'Fontname', 'Times New Roman','FontSize',18);

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
        title({strcat('True vs ',model2Use, '-PSO solution (1DSec RP):'),'Optimum rebar area of a Beam',...
            'Right section'},interpreter='latex') 
        hold on
        grid on
        axis([0 3000 0 2500])
        set(gca, 'Fontname', 'Times New Roman','FontSize',18);
        
        Rcoef(reps,1)=mean([BT1(1),BT3(1),BT3(1)]);
        
    end
    
    Rcoefv=[Rcoefv;Rcoef'];
    
end

YAo=[YAo1;YAo2;YAo3];
Ypv=[Ypvt1;Ypvt2;Ypvt3];
MSE=mse(YAo,Ypv)

%% Function appendix
function [XTrain,XValidation,XTest,ATrain,AValidation,ATest]=dataTrainNLClass...
    (DR,idxTrain,idxValidation,idxTest,AsTrain,meanX,sigsqX)

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
muX1 = mean(XTrain1);
sigsqX1 = var(XTrain1,1);
XTrain1 = (XTrain1 - muX1)./sqrt(sigsqX1);

muX2 = mean(XTrain2);
sigsqX2 = var(XTrain2,1);
XTrain2 = (XTrain2 - muX2)./sqrt(sigsqX2);

muX3 = mean(XTrain3);
sigsqX3 = var(XTrain3,1);
XTrain3 = (XTrain3 - muX3)./sqrt(sigsqX3);

muX4 = mean(XTrain4);
sigsqX4 = var(XTrain4,1);
XTrain4 = (XTrain4 - muX4)./sqrt(sigsqX4);

muX5 = mean(XTrain5);
sigsqX5 = var(XTrain5,1);
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
