function [volRebarSpans,LenRebarL,LenRebarM,LenRebarR,sepRebarSpans,db9Spans,...
  EffSpans,MrSpans,cSpans,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
  DistrRebarLeft,DistrRebarMid,DistrRebarRight,nbLMRspan,totalnbSpan,CFAspans]=...
  OptimMSFSBeamsRebar3DSec(b,h,span,brec,hrec,hagg,pmin,pmax,rebarAvailable,fcu,...
  load_conditions,fy,wac,cutxLoc,Wfac,pltOptConv,pltOptReb,span2Plt,nparticles,...
   niter)

nspans=length(load_conditions(:,1));

sepRebarSpans=[];
ListRebarDiamLeft=[];
ListRebarDiamMid=[];
ListRebarDiamRight=[];
DistrRebarLeft=[];
DistrRebarMid=[];
DistrRebarRight=[];
nbLMRspan=[];
nbcut3secSpan=[];
LenRebarL=[];
LenRebarM=[];
LenRebarR=[];

% Consequetive spans
bestnbtLMR=zeros(1,9);
bestdbc=zeros(1,9);

nblowRight=zeros(1,3);
dblowRight=zeros(1,3);
for i=1:nspans
    [bestPerformance,bestLenRebar,bestsepRebar,bestPosition,bestdbc,bestEfL,...
    bestEffMid,bestEffRight,bestMrLeft,bestMrMid,bestMrRight,cbestLeft,cBestMid,...
    cBestRight,bestListRebarDiamLeft,bestListRebarDiamMid,bestListRebarDiamRight,...
    bestDistrRebarLeft,bestRebarDistrMid,bestDistrRebarRight,bestnbtLMR,...
    bestnbcut3sec,nblowRight,dblowRight,CFA]=PSOBeamsRebar3DSec(b(i),h(i),...
    span(i),brec,hrec,hagg,pmin,pmax,rebarAvailable,fcu(i),load_conditions(i,:),...
    fy,wac,cutxLoc(i,:),bestdbc(1,7:9),bestnbtLMR(1,7:9),dblowRight,...
    nblowRight,Wfac,pltOptConv,i,nparticles,niter);
    
    nbcut3secSpan=[nbcut3secSpan;
                    bestnbcut3sec];
    db9Spans(i,:)=bestdbc;
    volRebarSpans(i,1)=bestPerformance;
    sepRebarSpans=[sepRebarSpans;
                   bestsepRebar];
    
    EffSpans(i,:)=[bestEfL,bestEffMid,bestEffRight];

    MrSpans(i,:)=[bestMrLeft,bestMrMid,bestMrRight];
    cSpans(i,:)=[cbestLeft,cBestMid,cBestRight];
	
    nbleft=length(bestListRebarDiamLeft);
    ListRebarDiamLeft=[ListRebarDiamLeft;
                      bestListRebarDiamLeft];
					  
    nbMid=length(bestListRebarDiamMid);
    ListRebarDiamMid=[ListRebarDiamMid;
                      bestListRebarDiamMid];
					  
    nbright=length(bestListRebarDiamRight);
    ListRebarDiamRight=[ListRebarDiamRight;
                        bestListRebarDiamRight];

    DistrRebarLeft=[DistrRebarLeft;
                    bestDistrRebarLeft];

    DistrRebarMid=[DistrRebarMid;
                   bestRebarDistrMid];

    DistrRebarRight=[DistrRebarRight;
                   bestDistrRebarRight];
    
    LenRebarL=[LenRebarL;
                bestLenRebar(1:sum(bestnbtLMR(1:3)))];
    LenRebarM=[LenRebarM;
                bestLenRebar(sum(bestnbtLMR(1:3))+1:sum(bestnbtLMR(1:6)))];
    LenRebarR=[LenRebarR;
                bestLenRebar(sum(bestnbtLMR(1:6))+1:sum(bestnbtLMR))];

    nbLMRspan=[nbLMRspan;
               bestnbtLMR];

    totalnbSpan(i,:)=[nbleft,nbMid,nbright];

    CFAspans(i,1)=CFA;
    if sum(bestdbc)==0
        disp(strcat('No solution was found for span No.',num2str(i)))
        return;
    end
end

nspanplt=length(span2Plt);
if pltOptReb==1
    for i=1:nspanplt
        %% Plot designs
        spanPlt=span2Plt(i);
        % Distribution of rebars - Left section

        nbleft=totalnbSpan(spanPlt,1);

        nbL=sum(totalnbSpan(1:spanPlt-1,1));
        Mul=load_conditions(spanPlt,2);
        
        if Mul<0
            DistrRebarLeftSec=-DistrRebarLeft(nbL+1:nbL+nbleft,:);
        else
            DistrRebarLeftSec=DistrRebarLeft(nbL+1:nbL+nbleft,:);
        end
        ListRebarDiamLeftSec=ListRebarDiamLeft(nbL+1:nbL+nbleft,1);

        % Distribution of rebars - Mid section

        nbmid=totalnbSpan(spanPlt,2);

        nbM=sum(totalnbSpan(1:spanPlt-1,2));

        DistrRebarMidSec=DistrRebarMid(nbM+1:nbM+nbmid,:);
        ListRebarDiamMidSec=ListRebarDiamMid(nbM+1:nbM+nbmid,1);

        % Distribution of rebars - Right section

        nbright=totalnbSpan(spanPlt,3);
        nbR=sum(totalnbSpan(1:spanPlt-1,3));

        Mur=load_conditions(spanPlt,4);

        if Mur<0
            DistrRebarRightSec=-DistrRebarRight(nbR+1:nbR+nbright,:);
        else
            DistrRebarRightSec=DistrRebarRight(nbR+1:nbR+nbright,:);
        end
        ListRebarDiamRightSec=ListRebarDiamRight(nbR+1:nbR+nbright,1);

        %% Plot
        plotBeamBar3sec(b(spanPlt),h(spanPlt),DistrRebarLeftSec,ListRebarDiamLeftSec,...
                DistrRebarMidSec,ListRebarDiamMidSec,DistrRebarRightSec,...
                ListRebarDiamRightSec,spanPlt);
    end
end
