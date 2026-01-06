function [distrRebarNew,listRebarDiamsNew]=distrRebarRecBeamCuts3DSec...
                            (nb3,db3,b,h,brec,hrec,vSep,nbAfterCut3)

dvs=10;

hp=h-2*hrec-2*dvs-db3(1);
bp=b-2*brec-2*dvs-max([db3]);

nbl1=nb3(1);
nbl2=nb3(2);
nbl3=nb3(3);

dbl1=db3(1);
dbl2=db3(2);
dbl3=db3(3);

nb=sum(nb3); % total number of rebars in tension
distrRebar=zeros(nb,2);

sepLowLay=(bp)/((nbl1)-1); % rebar separation - center to center

%% Conventional distribution
%% Layer 1

xl1=linspace(-0.5*bp,0.5*bp,nbl1);
distrRebar(1:nbl1,1)=xl1';
distrRebar(1:nbl1,2)=-0.5*hp+0.5*dbl1;

listRebarDiams(1:nbl1,1)=dbl1;

%% Layer 2

if nbl2>1
    xl2=linspace(-0.5*bp,0.5*bp,nbl2);
    distrRebar(nbl1+1:nbl1+nbl2,1)=xl2';
    distrRebar(nbl1+1:nbl1+nbl2,2)=-0.5*hp+dbl1+...
                                   vSep+0.5*dbl2;
elseif nbl2==1
    distrRebar(nbl1+1:nbl1+nbl2,1)=0;
    distrRebar(nbl1+1:nbl1+nbl2,2)=-0.5*hp+dbl1+...
                                    vSep+0.5*dbl2;
end
listRebarDiams(nbl1+1:nbl1+nbl2,1)=dbl2;

%% Layer 3

if nbl3>1
    xl2=linspace(-0.5*bp,0.5*bp,nbl3);
    distrRebar(nbl1+nbl2+1:nbl1+nbl2+nbl3,1)=xl2';
    distrRebar(nbl1+nbl2+1:nbl1+nbl2+nbl3,2)=-0.5*hp+dbl1+...
                                              vSep+dbl2+...
                                              vSep+0.5*dbl3;
elseif nbl3==1
    distrRebar(nbl1+nbl2+1:nbl1+nbl2+nbl3,1)=0;
    distrRebar(nbl1+nbl2+1:nbl1+nbl2+nbl3,2)=-0.5*hp+dbl1+...
                                              vSep+dbl2+...
                                              vSep+0.5*dbl3;
end
listRebarDiams(nbl1+nbl2+1:nbl1+nbl2+nbl3,1)=dbl3;

%% Redistribution
%% Layer 1
nb1new=nbAfterCut3(1);

nbcut1=nbl1-nb1new;
startRemove=(nb1new)/2+1; % if cuts are to be executed, the bars in the
                          % middle outwards are removed
finishRemove=startRemove-1+nbcut1;
indecesKeep=[];
for i=1:nbl1
    if i<startRemove || i>finishRemove
        indecesKeep=[indecesKeep;i];
    end
end

distrRebarNew1(:,1)=distrRebar(indecesKeep,1);
distrRebarNew1(:,2)=distrRebar(indecesKeep,2);

listRebarDiamsNew1=listRebarDiams(indecesKeep,1);

%% Layer 2
nb2new=nbAfterCut3(2);
nbcut2=nbl2-nb2new;
startRemove=(nb2new)/2+1;
finishRemove=startRemove-1+nbcut2;
indecesKeep=[];
for i=1:nbl2
    if i<startRemove || i>finishRemove
        indecesKeep=[indecesKeep;i];
    end
end
indecesKeep=indecesKeep+nbl1;
distrRebarNew2(:,1)=distrRebar(indecesKeep,1);
distrRebarNew2(:,2)=distrRebar(indecesKeep,2);

listRebarDiamsNew2=listRebarDiams(indecesKeep,1);

%% Layer 3
nb3new=nbAfterCut3(3);
nbcut3=nbl3-nb3new;
startRemove=(nb3new)/2+1;
finishRemove=startRemove-1+nbcut3;
indecesKeep=[];
for i=1:nbl3
    if i<startRemove || i>finishRemove
        indecesKeep=[indecesKeep;i];
    end
end
indecesKeep=indecesKeep+nbl1+nbl2;
distrRebarNew3(:,1)=distrRebar(indecesKeep,1);
distrRebarNew3(:,2)=distrRebar(indecesKeep,2);

listRebarDiamsNew3=listRebarDiams(indecesKeep,1);

%% Assembling of distributions
distrRebarNew=[distrRebarNew1;
               distrRebarNew2;
               distrRebarNew3];
           
listRebarDiamsNew=[listRebarDiamsNew1;
                  listRebarDiamsNew2;
                  listRebarDiamsNew3];
              