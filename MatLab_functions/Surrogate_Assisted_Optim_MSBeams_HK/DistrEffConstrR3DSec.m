function [minArea,bestEff,bestMr,bestc,xBest,bestListDiam,bestRebarDistr,...
    isfeasible]=PlainDistrEffConstrR1DLayer(Mu,fcu,Es,h,b,hagg,brec,hrec,pmin,...
    pmax,nbAfterCut3,distrRebarComp,listRebarDiamComp,...
    nb3l,dbc,nbr)

[vSepMin,vSepMax]=sepMinMaxHK13(max(dbc),hagg,1);
vSepMin=25;

fy=Es*0.00217;

amin=pmin*b*h;
amax=pmax*b*h;

dvs=10; % default diameter of shear links
bp=b-2*brec-2*dvs;
db1=dbc(1);
db2=dbc(2);
db3=dbc(3);

ab1=pi/4*db1^2;
ab2=pi/4*db2^2;
ab3=pi/4*db3^2;

nb1=nbr(1);
nb2=nbr(2);
nb3=nbr(3);

nbCombo3=[nb1,nb2,nb3];

% To determine if such possible solution is reliable
areaRebar=ab1*nb1+ab2*nb2+ab3*nb3;

%% Distribution of rebars
if sum(nb3l-nbAfterCut3)==0 % if there were no cuts for the left section
    [distrRebarTen,listRebarDiamTen]=distrRebarRecBeam1DiamLayer(nbCombo3,...
        dbc,b,h,brec,hrec,vSepMin);
else
    nb3l2=[nbCombo3(1),nbCombo3(1),nbCombo3(1)];
    if nbCombo3(1)<nb3l(1)
        [distrRebarTen,listRebarDiamTen]=distrRebarRecBeamCuts1DiamLayer...
            (nb3l,dbc,b,h,brec,hrec,vSepMin,nbCombo3);
    else
        [distrRebarTen,listRebarDiamTen]=distrRebarRecBeamCuts1DiamLayer...
            (nb3l2,dbc,b,h,brec,hrec,vSepMin,nbCombo3);
    end
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

if all([Eff<1.0,areaRebar>=amin,areaRebar<=amax])
    
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


