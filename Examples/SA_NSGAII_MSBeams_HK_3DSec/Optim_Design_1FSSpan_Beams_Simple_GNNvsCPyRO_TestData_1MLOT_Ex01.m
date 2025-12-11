% Optim_Design_1FSSpan_Beams_Simple_GNNvsCPyRO_TestData_1MLOT_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To test pareto front dominance between NSGAII assisted by 
%    CPyRO-GraphNet-Beams and plain GNNs.
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
B=importdata('Enhanced_Data_5LOT_HK_Nb_Db_Simple_4000.xlsx');
DR3=B.data;


DR2=A.data;
numObservations=size(DR2,1);
numObservations1=33;
numObservations2=33;
beamNo=1;
i=0;
dL=100;
for j =numObservations1:numObservations2
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

        DRB(i,:)=DR3(j,:);
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
numObservations=size(DR,1);

meanXNL=[349.605734767025	633.225806451613	37.4534050179212	4505.55555555556	-26339479.7363044];
sigsqXNL=[4949.66534345653	11216.8343161059	31.3195134954590	996295.300677007	2.10064997682719e+15];

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

pccb=[1.2,0.9,128,190,0.04,0.105,7];

%% Constructability
hagg=20;
Wunb=[1.3,1.4,1];
Wnd=[1.2,0.8];
Wcut=[1.3,1.6,2];
Wfac=[Wunb,Wnd,Wcut];

%% Generalization parameters
NT=numObservations;

MIGD1=[];
tMIGD1=[];

genMIGD1=[];
genMIGD2=[];
MIGD2=[];
tMIGD2=[];
gemMIGD2=[];
DomPF=[];
CIPF=[];
DomPFCI=[];

nreps=1;
Nv=[numObservations];
ndl=length(Nv);

for nd=1:ndl  % loop on vector of data size
    
    nodesAb1to3=u1L';
    N=Nv(nd);
    
    [idxTrain,idxValidation,idxTest] = trainingPartitions(N,[0 0 1]);
    
    % Target data
    AsDataTrain = nodesAb1to3(idxTrain,:);
    AsDataValidation = nodesAb1to3(idxValidation,:);
    AsDataTest = nodesAb1to3(idxTest,:);

    cutxLocBeamTrain=cutxLocBeam(idxTrain,:);
    cutxLocBeamValidation=cutxLocBeam(idxValidation,:);
    cutxLocBeamTest=cutxLocBeam(idxTest,:);
    
    % Data for NLayer classifier and/or GCNN
    [XTrainNL,XValidation,XTest,ATrain,AValidation,ATest,labelsTrain,...
    labelsValidation,labelsTest]=dataGNN...
        (DR,idxTrain,idxValidation,idxTest,AsDataTrain,AsDataValidation,...
        AsDataTest,meanXNL,sigsqXNL);
    
    YTest = AsDataTest;
    XTestOpt = DR(idxTest,1:7);

    %% Load surrogate models
    
    % Model for prediction of As per cross-section
    
    nheadsparamnGATPIGNN=load("C:/Users/luizv/OneDrive/CALDRECUST/Software/Package/CALDRECUST-MatLab/MatLab_functions/AI_CALDRECUST_MatLab/Surrogates_RecBeams/CPyRO_Graph_Net/nHeads_GAT_PIGNN_As_Section_4000.mat");
    paramPIGCNN=load("C:/Users/luizv/OneDrive/CALDRECUST/Software/Package/CALDRECUST-MatLab/MatLab_functions/AI_CALDRECUST_MatLab/Surrogates_RecBeams/CPyRO_Graph_Net/PIGCNN_As_Section_4000.mat");
    [XXTrain,XXValidation,XXTest] = dataPIGCNN(x0L,3,nodesAb1to3,idxTrain,...
        idxValidation,idxTest);

    Ao3CPyRO=model1fc2GAT1Conv1fc(paramPIGCNN.pignn,XTest,ATest,nheadsparamnGATPIGNN.numHeads);
    Ao3CPyRO=extractdata(Ao3CPyRO);

    nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
    paramGCNN=load("GCNN_As_Section_4000.mat");

    nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
    paramGCNN=load("GCNN_As_Section_4000.mat");

    Ao3GNN=GNNmodel1fc2GAT1Conv1fc(paramGCNN.parameters,XTest,ATest,nheadsparamnGATGNN.numHeads);
    Ao3GNN=extractdata(Ao3GNN);
    
    %Ao3=YTest;
    for reps=1:nreps % loop on data
        n1=1;
        n2=numObservations;
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

            Aos3CPyRO=Ao3CPyRO(i1:i2,1);
            Aos3GNN=Ao3GNN(i1:i2,1);
            
            nbcc=zeros(1,3);
            nblow=zeros(1,3);
            dbcc=zeros(1,3);
            dblow=zeros(1,3);
            pop_size=60;           % Population size
            gen_max=60;            % MAx number of generations - stopping criteria
            
            PF_CS_REF=1-DRB(i,100:104);
            PF_REF_WEIGHT=DRB(i,105:109).*wac;   

            MIGDconv=0;
            [extrOptPFArea,extrOptPFCS,PF_CFA1,PF_VOL1,newPop1,feasibleSol1,genCFA1,genVOL1,...
            IGDt1,IGDv1]=SANSGAIIMSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable([2:10]',:),...
            fcu,load_conditions,fy,wac,cutxLoc,dbcc,nbcc,dblow,nblow,PF_REF_WEIGHT,...
            PF_CS_REF,Wfac,Aos3CPyRO,pccb,MIGDconv,pop_size,gen_max,0);

            [extrOptPFArea,extrOptPFCS,PF_CFA2,PF_VOL2,newPop2,feasibleSol2,genCFA2,genVOL2,...
            IGDt2,IGDv2]=SANSGAIIMSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable([2:10]',:),...
            fcu,load_conditions,fy,wac,cutxLoc,dbcc,nbcc,dblow,nblow,PF_REF_WEIGHT,...
            PF_CS_REF,Wfac,Aos3GNN,pccb,MIGDconv,pop_size,gen_max,0);
            
            maxGen1=size(genCFA1,2);
            maxGen2=size(genCFA2,2);
            
            model2Use="CPyRO Graph-Net";
            plotEvolutionPF(PF_CS_REF,PF_REF_WEIGHT,genCFA1,genVOL1,IGDt1,...
                IGDv1,wac,10,20,30,maxGen1,model2Use,beamNo,10);
            
            model2Use="Plain GCNN";
            plotEvolutionPF(PF_CS_REF,PF_REF_WEIGHT,genCFA2,genVOL2,IGDt2,...
                IGDv2,wac,10,20,30,maxGen2,model2Use,beamNo,11);
            
            if all([sum(feasibleSol1(:,22)==0)>0, sum(feasibleSol2(:,22)==0)>0])
                [CI,PFdom]=dominancePFs(PF_CFA1,PF_VOL1,PF_CFA2,PF_VOL2);
                
                if CI>0.6
                    PfDomCI = 2;
                elseif CI < 0.4
                    PfDomCI = 1;
                else
                    PfDomCI = 0;
                end
            elseif all([sum(feasibleSol1(:,22)==0)>0, sum(feasibleSol2(:,22)==0)==0])
                PfDomCI=1;
                PfDom = 1;
                CI=0;
            elseif all([sum(feasibleSol1(:,22)==0)==0, sum(feasibleSol2(:,22)==0)>0])
                PfDomCI=0;
                PfDom = 2;
                CI=1;
            elseif all([sum(feasibleSol1(:,22)==0)==0, sum(feasibleSol2(:,22)==0)==0])
                PfDomCI=0.5;
                PfDom = 0;
                CI=0.5;
            end
            DomPF=[DomPF;PFdom];
            CIPF=[CIPF;CI];
            DomPFCI = [DomPFCI; PfDomCI];
            MIGD1 = [MIGD1;IGDv1(maxGen1,1)];
            tMIGD1 = [tMIGD1;IGDt1(maxGen1,1)];
            genMIGD1=[genMIGD1;maxGen1];

            MIGD2 = [MIGD2;IGDv2(maxGen2,1)];
            tMIGD2 = [tMIGD2;IGDt2(maxGen2,1)];
            genMIGD2=[genMIGD2;maxGen2];

        end
    end
end

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
    IGDv,wac,nstep1,nstep2,nstep3,nstep4,surrogate,beamNo,nfig)
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

figure(nfig)
title({strcat(surrogate,'-',' NSGAII');strcat('Beam ',num2str(beamNo))})
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
legend(strcat('Reference PF'),...
        strcat('PF-Gen-',num2str(nstep1),' (Time: ',num2str(IGDt(nstep1,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep2),' (Time: ',num2str(IGDt(nstep2,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep3), ' (Time: ',num2str(IGDt(nstep3,1)),' sec )'),...
       strcat('PF-Gen-',num2str(nstep4), ' (Time: ',num2str(IGDt(nstep4,1)),' sec )'),...
       'Location','northwest')
grid on
set(gca, 'Fontname', 'Times New Roman','FontSize',22);

end

function [XTrain,XValidation,XTest,ATrain,AValidation,ATest,labelsTrain,...
    labelsValidation,labelsTest]=dataGNN(DR,idxTrain,...
            idxValidation,idxTest,AsTrain,AsValidation,AsDataTest,meanXNL,sigsqXNL)

sigsqX1=sigsqXNL(1);
sigsqX2=sigsqXNL(2);
sigsqX3=sigsqXNL(3);
sigsqX4=sigsqXNL(4);
sigsqX5=sigsqXNL(5);

muX1=meanXNL(1);
muX2=meanXNL(2);
muX3=meanXNL(3);
muX4=meanXNL(4);
muX5=meanXNL(5);

numObservations=length(DR(:,1));

elements=[1 2 3 ;
          2 3 1];

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

[ATrain,XTrain1,labelsTrain] = preprocessData(adjacencyDataTrain,coulombDataTrain1,AsTrain);
[~,XTrain2,~] = preprocessData(adjacencyDataTrain,coulombDataTrain2,AsTrain);
[~,XTrain3,~] = preprocessData(adjacencyDataTrain,coulombDataTrain3,AsTrain);
[~,XTrain4,~] = preprocessData(adjacencyDataTrain,coulombDataTrain4,AsTrain);
[~,XTrain5,~] = preprocessData(adjacencyDataTrain,coulombDataTrain5,AsTrain);

% Validation partition
[AValidation,XValidation1,labelsValidation] = preprocessData(adjacencyDataValidation,coulombDataValidation1,AsValidation);
[~,XValidation2,~] = preprocessData(adjacencyDataValidation,coulombDataValidation2,AsValidation);
[~,XValidation3,~] = preprocessData(adjacencyDataValidation,coulombDataValidation3,AsValidation);
[~,XValidation4,~] = preprocessData(adjacencyDataValidation,coulombDataValidation4,AsValidation);
[~,XValidation5,~] = preprocessData(adjacencyDataValidation,coulombDataValidation5,AsValidation);

%% Normalizing training data

%muX1 = mean(XTrain1);
%sigsqX1 = var(XTrain1,1);
XTrain1 = (XTrain1 - muX1)./sqrt(sigsqX1);

%muX2 = mean(XTrain2);
%sigsqX2 = var(XTrain2,1);
XTrain2 = (XTrain2 - muX2)./sqrt(sigsqX2);

%muX3 = mean(XTrain3);
%sigsqX3 = var(XTrain3,1);
XTrain3 = (XTrain3 - muX3)./sqrt(sigsqX3);

%muX4 = mean(XTrain4);
%sigsqX4 = var(XTrain4,1);
XTrain4 = (XTrain4 - muX4)./sqrt(sigsqX4);

%muX5 = mean(XTrain5);
%sigsqX5 = var(XTrain5,1);
XTrain5 = (XTrain5 - muX5)./sqrt(sigsqX5);

XTrain=[XTrain1,XTrain2,XTrain3,XTrain4,XTrain5];

%meanXNL=[muX1,muX2,muX3,muX4,muX5];
%sigsqXNL=[sigsqX1,sigsqX2,sigsqX3,sigsqX4,sigsqX5];

%% Normalizing validation data
XValidation1 = (XValidation1 - muX1)./sqrt(sigsqX1);
XValidation2 = (XValidation2 - muX2)./sqrt(sigsqX2);

XValidation3 = (XValidation3 - muX3)./sqrt(sigsqX3);
XValidation4 = (XValidation4 - muX4)./sqrt(sigsqX4);
XValidation5 = (XValidation5 - muX5)./sqrt(sigsqX5);

XValidation=[XValidation1,XValidation2,XValidation3,XValidation4,XValidation5];

%% Normalizing test data
[ATest,XTest1,labelsTest] = preprocessData(adjacencyDataTest,coulombDataTest1,AsDataTest);
XTest1 = (XTest1 - muX1)./sqrt(sigsqX1);
XTest1 = dlarray(XTest1);

[~,XTest2,~] = preprocessData(adjacencyDataTest,coulombDataTest2,AsDataTest);
XTest2 = (XTest2 - muX2)./sqrt(sigsqX2);
XTest2 = dlarray(XTest2);

[~,XTest3,~] = preprocessData(adjacencyDataTest,coulombDataTest3,AsDataTest);
XTest3 = (XTest3 - muX3)./sqrt(sigsqX3);
XTest3 = dlarray(XTest3);

[~,XTest4,~] = preprocessData(adjacencyDataTest,coulombDataTest4,AsDataTest);
XTest4 = (XTest4 - muX4)./sqrt(sigsqX4);
XTest4 = dlarray(XTest4);

[~,XTest5,~] = preprocessData(adjacencyDataTest,coulombDataTest5,AsDataTest);
XTest5 = (XTest5 - muX5)./sqrt(sigsqX5);
XTest5 = dlarray(XTest5);

XTest=[XTest1,XTest2,XTest3,XTest4,XTest5];
end


function Y = GNNmodel1fc2GAT1Conv1fc(parameters,X,A,numHeads)

    
    ANorm = normalizeAdjacency(A);
    
    Z1 = X * parameters.Embedding.Weights + parameters.Embedding.b;
    
    weights1 = parameters.attn1.Weights;
    numHeadsAttention1 = numHeads.attn1;

    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = relu(Z2);

    weights2 = parameters.attn2.Weights;
    numHeadsAttention2 = numHeads.attn2;

    [Z3,~] = graphAttention(Z2,A,weights2,numHeadsAttention2,"cat");
    Z3  = relu(Z3) + Z2;

    Z4 = single(full(ANorm)) * Z3 * double(parameters.mult1.Weights);
    Z4 = relu(Z4);

    Z5 = Z4 * parameters.Decoder.Weights + parameters.Decoder.b;
    
    Y = Z5;
end


function Y = model1fc2GAT1Conv1fc(parameters,X,A,numHeads)

    
    ANorm = normalizeAdjacency(A);
    
    Z1 = X * parameters.Embed.Weights + parameters.Embed.b;
    
    weights1 = parameters.attn1.Weights;
    numHeadsAttention1 = numHeads.attn1;

    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = relu(Z2);

    weights2 = parameters.attn2.Weights;
    numHeadsAttention2 = numHeads.attn2;

    [Z3,~] = graphAttention(Z2,A,weights2,numHeadsAttention2,"cat");
    Z3  = relu(Z3) + Z2;

    Z4 = single(full(ANorm)) * Z3 * double(parameters.mult1.Weights);
    Z4 = relu(Z4);

    Z5 = Z4 * parameters.Decoder.Weights + parameters.Decoder.b;
    
    Y = Z5;
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
    normAttentionCoeff = softmax(attentionCoefficients,DataFormat="BCU");
    
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


function Y = modelNL(paramnNL,X,A,numHeadsNL)

    weights1 = paramnNL.attn1.Weights;
    numHeadsAttention1 = numHeadsNL.attn1;
    
    Z1 = X;
    [Z2,~] = graphAttention(Z1,A,weights1,numHeadsAttention1,"cat");
    Z2  = elu(Z2);

    weights2 = paramnNL.attn2.weights;
    numHeadsAttention2 = numHeadsNL.attn2;
    
    [Z3,~] = graphAttention(Z2,A,weights2,numHeadsAttention2,"cat");
    Z3  = elu(Z3) + Z2;
    
    ANorm = normalizeAdjacency(A);
    Z4 = single(full(ANorm)) * Z3 * double(paramnNL.mult.Weights);
    
    Z4 = graphPoolingLayer(Z4, 'mean', 3, 3);

    Y = softmax(Z4,"BC");
    Y = double(gather(extractdata(Y)));

    classes = {'1';'2';'3'};
    
    Y = onehotdecode(Y,classes,2);
    Y = dlarray(double(Y));
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

function [XTrain,XValidation,XTest] = dataPIGCNN(X,numNodesGNN,nodesAb1to3,...
    idxTrain,idxValidation,idxTest)


numObservations=length(idxTrain)+length(idxValidation)+length(idxTest);

elements=[1 2 3 ;
          2 3 1];

% Adjacency matrix
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

coulombData1=zeros(numObservations,numNodesGNN,numNodesGNN);
for i=1:numObservations
    for j=1:numNodesGNN
        coulombData1(i,j,j)=X(j,i);
    end
end
coulombData1 = double(permute(coulombData1, [2 3 1]));

% node adjacency data
adjacencyDataTrain = adjacency(:,:,idxTrain);
adjacencyDataValidation = adjacency(:,:,idxValidation);
adjacencyDataTest = adjacency(:,:,idxTest);

% feature data
coulombDataTrain1 = coulombData1(:,:,idxTrain);
coulombDataValidation1 = coulombData1(:,:,idxValidation);
coulombDataTest1 = coulombData1(:,:,idxTest);

% target data
AsDataTrain = nodesAb1to3(idxTrain,:);
AsDataValidation = nodesAb1to3(idxValidation,:);
AsDataTest = nodesAb1to3(idxTest,:);

% Train partition
[ATrain,XTrain1,labelsTrain] = preprocessData(adjacencyDataTrain,coulombDataTrain1,AsDataTrain);

% Validation partition
[AValidation,XValidation1,labelsValidation] = preprocessData(adjacencyDataValidation,coulombDataValidation1,AsDataValidation);

% Test partition
[ATest,XTest1,labelsTest] = preprocessData(adjacencyDataTest,coulombDataTest1,AsDataTest);

XTrain=[XTrain1];
XValidation=[XValidation1];
XTest=[XTest1];
end


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


function y = elu(x)

y = max(0, x) + (exp(min(0, x)) -1);

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
