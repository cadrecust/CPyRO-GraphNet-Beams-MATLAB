function [LenRebar1,LenRebar2,LenRebar3]=lenRebarRecBeamCuts1DiamLayer(nb3,....
        nbAfterCut3,span,ldL3)

nbl1=nb3(1);
nbl2=nb3(2);
nbl3=nb3(3);

%% Redistribution
%% Layer 1

nb1new=nbAfterCut3(1);

nbcut1=nbl1-nb1new;
startRemove=(nb1new)/2+1; % if cuts are to be executed, the bars in the
                          % middle outwards are removed
finishRemove=startRemove-1+nbcut1;
indecesKeep=[];
for i=1:nbl2
    if i<startRemove || i>finishRemove
        indecesKeep=[indecesKeep;i];
    end
end

LenRebar1=zeros(nbl1,1)+ldL3(1);
LenRebar1(indecesKeep,1)=span;

%% Layer 2
% Ends
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

LenRebar2=zeros(nbl2,1)+ldL3(2);
LenRebar2(indecesKeep,1)=span;

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

LenRebar3=zeros(nbl3,1)+ldL3(3);
LenRebar3(indecesKeep,1)=span;

