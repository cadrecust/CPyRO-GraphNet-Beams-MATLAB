function [extrOptPFCS,PF_CFA,PF_VOL,newPop,feasibleSol,genCFA,genVOL,IGDt,IGDv]=...
    SANSGAIIMSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable,...
    fcu,load_conditions,fy,wac,cutLoc,dbcc,nbcc,dblow,nblow,W_PF_REF,...
    CFA_PF_REF,Wfac,Ao3,ucfactors,MIGDconv,pop_size,gen_max)
NMLOCT=size(Ao3,2);
nspans=length(b);
yp_max=max(W_PF_REF);
np_baseline=length(CFA_PF_REF(:,1));

%% Construction performance 
ucwire=ucfactors(1);
uclinks=ucfactors(2);
ucrebar=ucfactors(3);
quantwire=ucfactors(7);
quantlinks=ucfactors(8);
ucManCutBend=ucfactors(4);
ucLorry=ucfactors(6);
quantLorry=ucfactors(9);
perfManCutBend=ucfactors(10);
ucManTyeAssem=ucfactors(5);
perfManTyeAssem=ucfactors(11);

%% Variables' range

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

%% OPTIMIZATION PARAMETERS

V=6; % number of optimization variables
M=2+7+9*nspans+9*nspans; % number of optimization functions + most important features of 
                         % each individual

% Upper and Lower bound variable range
dimax=[xmaxL,xmaxR];
dimin=[xminL,xminR];

txl=zeros(1,V)+dimin;
txu=zeros(1,V)+dimax;
xl=(txl(1,1:V));            % lower bound vector
xu=(txu(1,1:V));            % upper bound vectorfor 
etac = 20;                  % distribution index for crossover
etam = 20;                  % distribution index for mutation / mutation constant
          
pm=1/V;                     % Mutation Probability

Q=[];

%% Initial population 
xl_temp=repmat(xl, pop_size,1);
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V)); % x
for i =1:pop_size
    %% Design optimization
    % Optimization design
    db1l=rebarAvailable(1+fix(x(i,1)),2);
    db2l=rebarAvailable(1+fix(x(i,2)),2);
    db3l=rebarAvailable(1+fix(x(i,3)),2);

    db1m=rebarAvailable(1+fix(x(i,4)),2);
    db2m=rebarAvailable(1+fix(x(i,5)),2);
    db3m=rebarAvailable(1+fix(x(i,6)),2);

    [sepMin1l,sepMax1l]=sepMinMaxHK13([db1l],hagg,0);
    [sepMin2l,sepMax2l]=sepMinMaxHK13([db2l],hagg,0);
    [sepMin3l,sepMax3l]=sepMinMaxHK13([db3l],hagg,0);

    [sepMin1m,sepMax1m]=sepMinMaxHK13([db1m],hagg,0);
    [sepMin2m,sepMax2m]=sepMinMaxHK13([db2m],hagg,0);
    [sepMin3m,sepMax3m]=sepMinMaxHK13([db3m],hagg,0);

    sepMin=[sepMin1l,sepMin2l,sepMin3l,...
            sepMin1m,sepMin2m,sepMin3m];
    
    rm=ceil(rand*NMLOCT);
    Aos3=Ao3(:,rm);
    dbc=[db1l,db2l,db3l,db1m,db2m,db3m];
    if all([dbc(1)>=dbc(2),dbc(2)>=dbc(3),...
             dbc(4)>=dbc(5),dbc(5)>=dbc(6)])

        [rebarVol(i),volRebarSpans(i,:),LenRebarL{i},LenRebarM{i},LenRebarR{i},sepRebarSpan(:,:,i),...
        EffSpans(:,:,i),MrSpans(:,:,i),cSpans(:,:,i),ListRebarDiamLeft{i},...
        ListRebarDiamMid{i},ListRebarDiamRight{i},DistrRebarLeft{i},RebarDistrMid{i},...
        DistrRebarRight{i},nbLMRspan(:,:,i),totalnbSpan(:,:,i),CSrebarSpans(i,:),...
        CSAssemCutSpan(:,:,i),nbcut3sec(:,:,i),nblowLeft(:,:,i),nbTopMid(:,:,i),nblowRight(:,:,i),bestCFA(i),...
        cc(i,1)]=SAMOOptimMSFSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,...
        sepMin,fcu,load_conditions,fy,cutLoc,Aos3,Wfac,dbc);
    
    else
        rebarVol(i)=1e10;
        volRebarSpans(i,:)=zeros(1,nspans)+1e10;
        sepRebarSpan(:,:,i)=zeros(3*nspans,3);
        EffSpans(:,:,i)=zeros(nspans,3);
        MrSpans(:,:,i)=zeros(nspans,3);

        RebarDistrMid{i}=[];
        ListRebarDiamMid{i}=0;

        DistrRebarLeft{i}=[];
        ListRebarDiamLeft{i}=0;

        DistrRebarRight{i}=[];
        ListRebarDiamRight{i}=0;
        cSpans(:,:,i)=zeros(nspans,3);
        CSrebarSpans(i,:)=zeros(1,nspans);
        CSAssemCutSpan(:,:,i)=zeros(nspans,2);
        nblowLeft(:,:,i)=zeros(nspans,3);
        nbTopMid(:,:,i)=zeros(nspans,3);
        nblowRight(:,:,i)=zeros(nspans,3);
        nbLMRspan(:,:,i)=zeros(nspans,9);
        nbcut3sec(:,:,i)=zeros(3*nspans,3);
        bestCFA(i)=0;
        cc(i,1)=10;
    end
    dbc6(i,:)=dbc;
    
    for k=1:nspans
        i1=(k-1)*9+1;
        i2=(k)*9;

        % number of bars in tension
        NbCombo9span(i,i1:i2)=nbLMRspan(k,:,i);

        % number of bars in compression
        nb2c9span(i,i1:i2)=[nblowLeft(k,:,i),nbTopMid(k,:,i),nblowRight(k,:,i)];

        %% Construction cost, considering constructability with cuts

        CS1=CSAssemCutSpan(k,1,i);
        CS2=CSAssemCutSpan(k,2,i);
        puRecBeam=UCostRebarPrefab1(ucwire,uclinks,ucrebar,quantwire,quantlinks,...
                          ucManCutBend,ucLorry,quantLorry,perfManCutBend/CS2,...
                          ucManTyeAssem,perfManTyeAssem/CS1);
        
        bestCostSpan(i,k)=volRebarSpans(i,k)*wac*puRecBeam;
    end
    bestCost(i,1)=sum(bestCostSpan(i,:));

    %% Objective functions
    ff(i,:)=[rebarVol(i) 1-bestCFA(i)];
    
    %% Evlauation of constraints violation
    err(i,:)=(cc(i,1)>0).*cc(i,1);
end

ff2=[];
for i=1:pop_size
    ff2=[ff2;
        [ff(i,:), bestCost(i,1), NbCombo9span(i,:), dbc6(i,:), cc(i,1)]];
end

error_norm=normalisation(err); % Normalisation of the constraint violation

optimal_features_ff=[bestCost NbCombo9span dbc6 nb2c9span];

population_init=[x ff optimal_features_ff error_norm];
[population front]=NDS_CD_cons(population_init,V,M);    % Non domination Sorting on initial population

IGDv=zeros(gen_max,1);
IGDt=zeros(gen_max,1);

t1=datetime;
%% Generation Starts
for gen_count=1:gen_max
    % selection (Parent Pt of 'N' pop size)
    parent_selected=tour_selection(population);                     % 10 Tournament selection
    %%% Reproduction (Offspring Qt of 'N' pop size)
    child_offspring  = genetic_operator(parent_selected(:,1:V),V,xl,xu,etac,etam,pm);    % SBX crossover and polynomial mutation

    for ii = 1:length(child_offspring(:,1))
        %% Design optimization
        % Optimization design
        db1l=rebarAvailable(1+fix(child_offspring(ii,1)),2);
        db2l=rebarAvailable(1+fix(child_offspring(ii,2)),2);
        db3l=rebarAvailable(1+fix(child_offspring(ii,3)),2);

        db1m=rebarAvailable(1+fix(child_offspring(ii,4)),2);
        db2m=rebarAvailable(1+fix(child_offspring(ii,5)),2);
        db3m=rebarAvailable(1+fix(child_offspring(ii,6)),2);

        [sepMin1l,sepMax1l]=sepMinMaxHK13([db1l],hagg,0);
        [sepMin2l,sepMax2l]=sepMinMaxHK13([db2l],hagg,0);
        [sepMin3l,sepMax3l]=sepMinMaxHK13([db3l],hagg,0);

        [sepMin1m,sepMax1m]=sepMinMaxHK13([db1m],hagg,0);
        [sepMin2m,sepMax2m]=sepMinMaxHK13([db2m],hagg,0);
        [sepMin3m,sepMax3m]=sepMinMaxHK13([db3m],hagg,0);

        sepMin=[sepMin1l,sepMin2l,sepMin3l,...
                sepMin1m,sepMin2m,sepMin3m];
        
        rm=ceil(rand*NMLOCT);
        Aos3=Ao3(:,rm);
        dbc=[db1l,db2l,db3l,db1m,db2m,db3m];
        if all([dbc(1)>=dbc(2),dbc(2)>=dbc(3),...
             dbc(4)>=dbc(5),dbc(5)>=dbc(6)])

            [rebarVol(ii),volRebarSpans(ii,:),LenRebarL{ii},LenRebarM{ii},LenRebarR{ii},...
            sepRebarSpan(:,:,ii),EffSpans(:,:,ii),MrSpans(:,:,ii),cSpans(:,:,ii),ListRebarDiamLeft{ii},...
            ListRebarDiamMid{ii},ListRebarDiamRight{ii},DistrRebarLeft{ii},RebarDistrMid{ii},...
            DistrRebarRight{ii},nbLMRspan(:,:,ii),totalnbSpan(:,:,ii),CSrebarSpans(ii,:),...
            CSAssemCutSpan(:,:,ii),nbcut3sec(:,:,ii),nblowLeft(:,:,ii),nbTopMid(:,:,ii),nblowRight(:,:,ii),bestCFA(ii),...
            cc(ii,1)]=SAMOOptimMSFSBeamsRebarSimple(b,h,span,brec,hrec,hagg,pmin,pmax,...
            sepMin,fcu,load_conditions,fy,cutLoc,Aos3,Wfac,dbc);

        else
            rebarVol(ii)=1e10;
            sepRebarSpan(:,:,ii)=zeros(3*nspans,3);
            EffSpans(:,:,ii)=zeros(nspans,3);
            MrSpans(:,:,ii)=zeros(nspans,3);

            RebarDistrMid{ii}=[];
            ListRebarDiamMid{ii}=0;

            DistrRebarLeft{ii}=[];
            ListRebarDiamLeft{ii}=0;

            DistrRebarRight{ii}=[];
            ListRebarDiamRight{ii}=0;
            cSpans(:,:,ii)=zeros(nspans,3);
            CSrebarSpans(ii,:)=zeros(1,nspans);
            CSAssemCutSpan(:,:,ii)=zeros(nspans,2);
            nblowLeft(:,:,ii)=zeros(nspans,3);
            nbTopMid(:,:,ii)=zeros(nspans,3);
            nblowRight(:,:,ii)=zeros(nspans,3);
            
            nbLMRspan(:,:,ii)=zeros(nspans,9);
            nbcut3sec(:,:,ii)=zeros(3*nspans,3);
            bestCFA(ii)=0;
            cc(ii,1)=10;
        end
        dbc6(ii,:)=dbc;

        for k=1:nspans
            i1=(k-1)*9+1;
            i2=(k)*9;

            % number of bars in tension
            NbCombo9span(ii,i1:i2)=nbLMRspan(k,:,ii);

            % number of bars in compression
            nb2c9span(ii,i1:i2)=[nblowLeft(k,:,ii),nbTopMid(k,:,ii),nblowRight(k,:,ii)];

            %% Construction cost
            
            CS1=CSAssemCutSpan(k,1,ii);
            CS2=CSAssemCutSpan(k,2,ii);
            puRecBeam=UCostRebarPrefab1(ucwire,uclinks,ucrebar,quantwire,quantlinks,...
                              ucManCutBend,ucLorry,quantLorry,perfManCutBend/CS2,...
                              ucManTyeAssem,perfManTyeAssem/CS1);
            
            bestCostSpan(ii,k)=volRebarSpans(ii,k)*wac*puRecBeam;
        end
        bestCost(ii,1)=sum(bestCostSpan(ii,:));

        %% Objective functions
        fff(ii,:)=[rebarVol(ii) 1-bestCFA(ii)];

        %% Evlauation of constraints violation
        err(ii,:)=(cc(ii,1)>0).*cc(ii,1);

    end
    for i=1:pop_size
        ff2=[ff2;
            [fff(i,:),bestCost(i,:),NbCombo9span(i,:), dbc6(i,:), cc(i,1)]];
    end

    optimal_features_fff=[bestCost NbCombo9span dbc6 nb2c9span];

    error_norm=normalisation(err);                                
    child_offspring=[child_offspring fff optimal_features_fff error_norm];

    %%% INtermediate population (Rt= Pt U Qt of 2N size)
    population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)];
    [population_inter_sorted front]=NDS_CD_cons(population_inter,V,M);              % Non domination Sorting on offspring
    %%% Replacement - N
    newPop=replacement(population_inter_sorted, front, pop_size);
    population=newPop;

    t2=datetime;
    tgen=seconds(t2-t1);

    genCFA(:,gen_count)=1-newPop(:,V+2);
    genVOL(:,gen_count)=newPop(:,V+1);

    %% Inverted Generational Distance

    MIGD=0;
    nonValid=0;
    for npi=1:pop_size
        xp=population(npi,V+2);
        yp=population(npi,V+1).*wac;
        minIGD=inf;
        for npj=1:np_baseline
            xp_base=1-CFA_PF_REF(npj);
            yp_base=W_PF_REF(npj);

            IGD=sqrt((xp-xp_base)^2+((yp-yp_base)/yp_max)^2);
            if IGD<minIGD
                minIGD=IGD;
            end
        end
        if minIGD>=1e5
            nonValid=nonValid+1;
        elseif minIGD<1e5
            MIGD=MIGD+minIGD;
        end
    end
    MIGD=MIGD/(pop_size-nonValid);

    if gen_count>1
        IGDv(gen_count,1)=(sum(IGDv(1:gen_count-1,1))+MIGD)/gen_count;
        dIGD=abs(IGDv(gen_count,1)-IGDv(gen_count-1,1));
    else
        IGDv(gen_count,1)=MIGD;
    end
    IGDt(gen_count,1)=tgen;

    if MIGDconv==1
        if gen_count>1
            if dIGD/max(IGDv)<0.01/100
                break;
            end
        end
    end
end
newPop=sortrows(newPop,V+1);

%% Collecting designs along the PF
cscore=newPop(:,8);

min2index=fix(1/4*pop_size);
midindex=fix(pop_size/2);
mid2index=fix(3/4*pop_size);

minOpt=cscore(1,:);
min2Opt=cscore(min2index,:);
midOpt=cscore(midindex,:);
mid2Opt=cscore(mid2index,:);
maxOpt=cscore(pop_size,:);

extrOptPFCS=[minOpt,min2Opt,midOpt,mid2Opt,maxOpt];

[PF_CFA,PF_VOL]=filterSingle(newPop(:,V+2),newPop(:,V+1));
PF_CFA=1-PF_CFA;

%% Results
% PARETO FRONT COST-CFA
feasibleSol=[];
nfeasible=size(ff2,1);
for i=1:nfeasible
    if ff2(i,2)<1
        feasibleSol=[feasibleSol;
                    ff2(i,:)];
    end
end

if sum(ff2(:,9+9*nspans+1)==0)==0
    disp('No solution was found')
end
end

%% Function appendix
function [parent_selected] = tour_selection(pool)
%%% Description

% 1. Parents are selected from the population pool for reproduction by using binary tournament selection
%    based on the rank and crowding distance. 
% 2. An individual is selected if the rank is lesser than the other or if
%    crowding distance is greater than the other.
% 3. Input and output are of same size [pop_size, V+M+3].


%%% Binary Tournament Selection
[pop_size, distance]=size(pool);
rank=distance-1;
candidate=[randperm(pop_size);randperm(pop_size)]';

for i = 1: pop_size
    parent=candidate(i,:);                                  % Two parents indexes are randomly selected
    if pool(parent(1),rank)~=pool(parent(2),rank)              % For parents with different ranks
        if pool(parent(1),rank)<pool(parent(2),rank)            % Checking the rank of two individuals
            mincandidate=pool(parent(1),:);
        elseif pool(parent(1),rank)>pool(parent(2),rank)
            mincandidate=pool(parent(2),:);
        end
        parent_selected(i,:)=mincandidate;                          % Minimum rank individual is selected finally
    else                                                       % for parents with same ranks  
        if pool(parent(1),distance)>pool(parent(2),distance)    % Checking the distance of two parents
            maxcandidate=pool(parent(1),:);
        elseif pool(parent(1),distance)< pool(parent(2),distance)
            maxcandidate=pool(parent(2),:);
        else
            temp=randperm(2);
            maxcandidate=pool(parent(temp(1)),:);
        end 
        parent_selected(i,:)=maxcandidate;                          % Maximum distance individual is selected finally
    end
end
end

function new_pop=replacement(population_inter_sorted, front, pop_size)

%% Description
% The next generation population is formed by appending each front subsequently until the
% population size exceeds the current population size. If When adding all the individuals
% of any front, the population exceeds the population size, then the required number of 
% remaining individuals alone are selected from that particular front based
% on crowding distance.
%% code starts
index=0;
ii=1;
while index < pop_size
    l_f=length(front(ii).fr);
    if index+l_f < pop_size 
        new_pop(index+1:index+l_f,:)= population_inter_sorted(index+1:index+l_f,:);
        index=index+l_f;
    else
            temp1=population_inter_sorted(index+1:index+l_f,:);
            temp2=sortrows(temp1,size(temp1,2));
            new_pop(index+1:pop_size,:)= temp2(l_f-(pop_size-index)+1:l_f,:);
            index=index+l_f;
    end
    ii=ii+1;
end
end

function mutated_child = poly_mutation(y,V, xl, xu, etam, pm)

%% Description
% 1. Input is the crossovered child of size (1,V) in the vector 'y' from 
%    'genetic_operator.m'.
% 2. Output is in the vector 'mutated_child' of size (1,V).
%% Polynomial mutation including boundary constraint
del=min((y-xl),(xu-y))./(xu-xl);
t=rand(1,V);
loc_mut=t<pm;        
u=rand(1,V);
delq=(u<=0.5).*((((2*u)+((1-2*u).*((1-del).^(etam+1)))).^(1/(etam+1)))-1)+...
    (u>0.5).*(1-((2*(1-u))+(2*(u-0.5).*((1-del).^(etam+1)))).^(1/(etam+1)));
c=y+delq.*loc_mut.*(xu-xl);
mutated_child=c;
end

function err_norm  = normalisation(error_pop)
  
%% Description
% 1. This function normalises the constraint violation of various 
%    individuals, since the range of constraint violation of every 
%    chromosome is not uniform.
% 2. Input is in the matrix error_pop with size [pop_size, number of 
%    constraints].
% 3. Output is a normalised vector, err_norm of size [pop_size,1]
 
%% Error Nomalisation
[N,nc]=size(error_pop);
con_max=0.001+max(error_pop);
con_maxx=repmat(con_max,N,1);
cc=error_pop./con_maxx;
err_norm=sum(cc,2);                % finally sum up all violations
end


function [chromosome_NDS_CD front] = NDS_CD_cons(population,V,M) 

%%%%% Description
% 1. This function is to perform Deb's fast elitist non-domination sorting and crowding distance assignment. 
% 2. Input is in the variable 'population' with size: [size(popuation), V+M+1]
% 3. This function returns 'chromosome_NDS_CD' with size [size(population),V+M+3]
% 4. A flag 'problem_type' is used to identify whether the population is fully feasible (problem_type=0) or fully infeasible (problem_type=1) 
%    or partly feasible (problem_type=0.5). 

%%% Reference:
%Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan, " A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II", 
%IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 6, No. 2, APRIL 2002. 


%%% function begins
% Initialising structures and variables
chromosome_NDS_CD1=[];
infpop=[];
front.fr=[];
struct.sp=[];
rank=1;


%%% Segregating feasible and infeasible solutions

if all(population(:,V+M+1)==0)
    problem_type=0;
    chromosome=population(:,1:V+M);                         % All Feasible chromosomes;    
    pop_size1=size(chromosome,1);
elseif all(population(:,V+M+1)~=0)
    problem_type=1;
    pop_size1=0;
    infchromosome=population;                               % All InFeasible chromosomes;       
else
    problem_type=0.5;
    feas_index=find(population(:,V+M+1)==0);
    chromosome=population(feas_index,1:V+M);                % Feasible chromosomes;    
    pop_size1=size(chromosome,1);
    infeas_index=find(population(:,V+M+1)~=0);
    infchromosome=population(infeas_index,1:V+M+1);         % infeasible chromosomes;    
end

%%% Handling feasible solutions 
if problem_type==0 || problem_type==0.5
    pop_size1 = size(chromosome,1);
    f1 = chromosome(:,V+1);   % objective function values
    f2 = chromosome(:,V+2);
    %Non- Domination Sorting 
    % First front
    for p=1:pop_size1
        struct(p).sp=find(((f1(p)-f1)<0 &(f2(p)-f2)<0) | ((f2(p)-f2)==0 &(f1(p)-f1)<0) | ((f1(p)-f1)==0 &(f2(p)-f2)<0)); 
        n(p)=length(find(((f1(p)-f1)>0 &(f2(p)-f2)>0) | ((f2(p)-f2)==0 &(f1(p)-f1)>0) | ((f1(p)-f1)==0 &(f2(p)-f2)>0)));
    end

    front(1).fr=find(n==0);
    % Creating subsequent fronts
    while (~isempty(front(rank).fr))

        front_indiv=front(rank).fr;
        n(front_indiv)=inf;
        chromosome(front_indiv,V+M+1)=rank;
        rank=rank+1;
        front(rank).fr=[];
        for i = 1:length(front_indiv)

            temp=struct(front_indiv(i)).sp;
            n(temp)=n(temp)-1;
        end 
        q=find(n==0);
        front(rank).fr=[front(rank).fr q];

    end
    chromosome_sorted=sortrows(chromosome,V+M+1);    % Ranked population  

    %Crowding distance Assignment
    rowsindex=1;
    for i = 1:length(front)-1
        l_f=length(front(i).fr);

        if l_f > 2

            sorted_indf1=[];
            sorted_indf2=[];
            sortedf1=[];
            sortedf2=[];
              % sorting based on f1 and f2;
            [sortedf1 sorted_indf1]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+1));
            [sortedf2 sorted_indf2]=sortrows(chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+2));

            f1min=chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+1);
            f1max=chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+1);

            chromosome_sorted(sorted_indf1(1)+rowsindex-1,V+M+2)=inf;
            chromosome_sorted(sorted_indf1(end)+rowsindex-1,V+M+2)=inf;

            f2min=chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+2);
            f2max=chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+2);

            chromosome_sorted(sorted_indf2(1)+rowsindex-1,V+M+3)=inf;
            chromosome_sorted(sorted_indf2(end)+rowsindex-1,V+M+3)=inf;

            for j = 2:length(front(i).fr)-1

                if  (f1max - f1min == 0) || (f2max - f2min == 0)

                    chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=inf;
                    chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=inf;
                else
                    chromosome_sorted(sorted_indf1(j)+rowsindex-1,V+M+2)=(chromosome_sorted(sorted_indf1(j+1)+rowsindex-1,V+1)-chromosome_sorted(sorted_indf1(j-1)+rowsindex-1,V+1))/(f1max-f1min);
                    chromosome_sorted(sorted_indf2(j)+rowsindex-1,V+M+3)=(chromosome_sorted(sorted_indf2(j+1)+rowsindex-1,V+2)-chromosome_sorted(sorted_indf2(j-1)+rowsindex-1,V+2))/(f2max-f2min);
                end
            end


        else

            chromosome_sorted(rowsindex:(rowsindex+l_f-1),V+M+2:V+M+3)=inf;
        end
        rowsindex = rowsindex + l_f;
    end
    chromosome_sorted(:,V+M+4) = sum(chromosome_sorted(:,V+M+2:V+M+3),2); 
    chromosome_NDS_CD1 = [chromosome_sorted(:,1:V+M) zeros(pop_size1,1) chromosome_sorted(:,V+M+1) chromosome_sorted(:,V+M+4)]; % Final Output Variable

end

%%% Handling infeasible solutions
if problem_type==1 | problem_type==0.5
    infpop=sortrows(infchromosome,V+M+1);
    infpop=[infpop(:,1:V+M+1) (rank:rank-1+size(infpop,1))' inf*(ones(size(infpop,1),1))];
    for kk = (size(front,2)):(size(front,2))+(length(infchromosome))-1
        front(kk).fr= pop_size1+1;
    end
end
chromosome_NDS_CD = [chromosome_NDS_CD1;infpop]; 
end

function child_offspring  = genetic_operator(parent_selected,V,xl,xu,etac,etam,pm)

%% Description
% 1. Crossover followed by mutation
% 2. Input is in 'parent_selected' matrix of size [pop_size,V].
% 3. Output is also of same size in 'child_offspring'. 

%% Reference 
% Deb & samir agrawal,"A Niched-Penalty Approach for Constraint Handling in Genetic Algorithms". 
%% SBX cross over operation incorporating boundary constraint
[N] = size(parent_selected,1);
xl1=xl';
xu1=xu';
rc=randperm(N);
for i=1:(N/2)
    parent1=parent_selected((rc(2*i-1)),:);
    parent2=parent_selected((rc(2*i)),:);
    if (isequal(parent1,parent2))==1 & rand(1)>0.5
        child1=parent1;
        child2=parent2;
    else 
        for j = 1: V  
            if parent1(j)<parent2(j)
               beta(j)= 1 + (2/(parent2(j)-parent1(j)))*(min((parent1(j)-xl1(j)),(xu1(j)-parent2(j))));
            else
               beta(j)= 1 + (2/(parent1(j)-parent2(j)))*(min((parent2(j)-xl1(j)),(xu1(j)-parent1(j))));
            end   
        end
         u=rand(1,V);
         alpha=2-beta.^-(etac+1);
         betaq=(u<=(1./alpha)).*(u.*alpha).^(1/(etac+1))+(u>(1./alpha)).*(1./(2 - u.*alpha)).^(1/(etac+1));
        child1=0.5*(((1 + betaq).*parent1) + (1 - betaq).*parent2);
        child2=0.5*(((1 - betaq).*parent1) + (1 + betaq).*parent2);
    end
    child_offspring((rc(2*i-1)),:)=poly_mutation(child1,V, xl, xu, etam, pm);           % polynomial mutation
    child_offspring((rc(2*i)),:)=poly_mutation(child2,V, xl, xu, etam, pm);             % polynomial mutation
end
end

function [vn1,vn2]=filterSingle(v1,v2)

    vn1=[v1(1)];
    vn2=[v2(1)];
    j=1;
    np=length(v1(:,1));
    for i=1:np-1
        if v1(i+1,1)~=v1(i,1)
            vn1=[vn1;v1(i+1,1)];
            vn2=[vn2;v2(i+1,1)];
        end
    end
end