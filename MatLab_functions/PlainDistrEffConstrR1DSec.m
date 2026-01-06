function [minArea,bestEff,bestMr,bestc,xBest,bestListDiam,bestRebarDistr,...
    isfeasible]=DistrEffConstrR1DSec(Mu,fcu,Es,h,b,hagg,brec,hrec,pmin,...
    pmax,sepMin,sepRebaright,nbAfterCut3,distrRebarComp,listRebarDiamComp,...
    nb3l,dbc,nbr)

[vSepMin,~]=sepMinMaxHK13(dbc(1),hagg,1);

fy=Es*0.00217;

amin=pmin*b*h;
amax=pmax*b*h;

dvs=10; % default diameter of shear links
bp=b-2*brec-2*dvs;
db1=dbc(1);

ab1=pi/4*db1^2;

nb1=nbr(1);
nb2=nbr(2);
nb3=nbr(3);

nbCombo3=[nb1,nb2,nb3];

% To determine if such possible solution is reliable
areaRebar=ab1*nb1+ab1*nb2+ab1*nb3;

%% Distribution of rebars
if sum(nb3l-nbAfterCut3)==0 % if there were cuts for the left section
    [distrRebarTen,listRebarDiamTen]=distrRebarRecBeam1DSec(nbCombo3,...
        dbc,b,h,brec,hrec,vSepMin);
else
    nb3l2=[nbCombo3(1),nbCombo3(1),nbCombo3(1)];
    [distrRebarTen,listRebarDiamTen]=distrRebarRecBeamCuts1DSec...
        (nb3l2,dbc,b,h,brec,hrec,vSepMin,nbCombo3);
end

if Mu<0
    distrRebar=[distrRebarTen;
                -distrRebarComp];
else
    distrRebar=[-distrRebarTen;
                distrRebarComp];
end
listRebarDiam=[listRebarDiamTen;
               listRebarDiamComp];
[Eff,Mr,c]=EfRecBeamBars(Mu,fcu,Es,fy,h,b,distrRebar,...
                             listRebarDiam,hrec);

[ccr]=rebarDistrConstr3LayerRecBeam1DiamLayer(bp,nbCombo3);

if all([Eff<1.0,areaRebar>=amin,areaRebar<=amax,sepRebaright(1)>=sepMin(1),ccr==1])
    
    minArea=areaRebar;
    xBest=nbCombo3;
    bestEff=Eff;
    bestMr=Mr;
    bestListDiam=listRebarDiam;
    bestRebarDistr=distrRebar;
    bestc=c;
    isfeasible=true;
else
    minArea=1e10;
    xBest=zeros(1,3);
    bestEff=0;
    bestMr=0;
    bestListDiam=0;
    bestRebarDistr=[];
    bestc=0;
    isfeasible=false;
end

