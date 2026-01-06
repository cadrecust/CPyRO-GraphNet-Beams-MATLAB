function [minimumFitness,bestLenRebar,bestsepRebar,bestNbCombo9,bestEffLeft,bestEffMid,...
    bestEffRight,bestMrLeft,bestMrMid,bestMrRight,bestcLeft,bestcMid,bestcRight,...
    bestListRebarDiamLeft,bestListRebarDiamMid,bestListRebarDiamRight,...
    bestDistrRebarLeft,bestRebarDistrMid,bestDistrRebarRight,...
    bestnbcut3sec,bestnblowLeft,bestdblowLeft,bestnbTopMid,bestdbTopMid,...
    bestnblowRight,bestdblowRight,bestCFA,constr]=SACutRedistrOptimRecBeam3DSec...
    (load_conditions,fcu,Es,h,b,span,dbc,hagg,brec,hrec,pmin,pmax,sepMin,...
    cutLoc,Ao3,Wfac,nblm)
            
fy=Es*0.00217;

amin=pmin*b*h;
amax=pmax*b*h;

dvs=10;
bpl=b-2*brec-2*dvs-max(dbc(1:3));
bpm=b-2*brec-2*dvs-max(dbc(4:6));
    
%% Rebar diameters
db1l=dbc(1);
db2l=dbc(2);
db3l=dbc(3);

db1m=dbc(4);
db2m=dbc(5);
db3m=dbc(6);

%% Rebar area

ab1l=pi/4*db1l^2;
ab2l=pi/4*db2l^2;
ab3l=pi/4*db3l^2;

ab1m=pi/4*db1m^2;
ab2m=pi/4*db2m^2;
ab3m=pi/4*db3m^2;

abl3=[ab1l,ab2l,ab3l];
abm3=[ab1m,ab2m,ab3m];

nbl1=nblm(1); nbl2=nblm(2); nbl3=nblm(3);
nbm1=nblm(4); nbm2=nblm(5); nbm3=nblm(6);

nb3l=nblm(1:3);
nb3m=nblm(4:6);

Abl=ab1l*nbl1+ab2l*nbl2+ab3l*nbl3;
Abm=ab1m*nbm1+nbm2*ab2m+nbm3*ab3m;

%% Rebar separation

% Left cross-section
sepRebarl1=(bpl)/((nbl1)-1);
if nbl2>1
    sepRebarl2=(bpl)/(nbl2-1);
else
    sepRebarl2=1e10;
end
if nbl3>1
    sepRebarl3=(bpl)/(nbl3-1);
else
    sepRebarl3=1e10;
end

% Mid cross-section
sepRebarm1=(bpm)/(nbm1-1);
if nbm2>1
    sepRebarm2=(bpm)/(nbm2-1);
else
    sepRebarm2=1e10;
end
if nbm3>1
    sepRebarm3=(bpm)/(nbm3-1);
else
    sepRebarm3=1e10;
end

sepRebarLeft=[sepRebarl1,sepRebarl2,sepRebarl3];
sepRebarMid=[sepRebarm1,sepRebarm2,sepRebarm3];

dbl3=[db1l,db2l,db3l];
dbm3=[db1m,db2m,db3m];

%% Cuts and redistribution Left-Mid span
[nbAfterCut3L,nbAfterCut3M,nbAfterCut3R,redistrRebarL2M,relistRebarDiamL2M,...
redistrRebarM2L,relistRebarM2L,redistrRebarM2R,relistRebarM2R]=...
cutRedistRebarSpanRecBeam3DSec(load_conditions,nb3l,nb3m,dbl3,dbm3,b,h,...
brec,hrec,hagg,fcu,fy);

%% Efficiency analysis
% Optimization of Right section
Mur=load_conditions(:,4);
[Abr,EffRight,MrRight,cRight,nbarsRight,sepRebarRight,ListDiamRight,...
 RebarDistrRight,isfeasibleRight]=SAPSOBeamsRebarR3DSec(Mur,fcu,h,b,...
 hagg,brec,hrec,pmin,pmax,sepRebarLeft,redistrRebarM2R,relistRebarM2R,...
 nbAfterCut3L,nb3l,dbl3,Ao3);

[sepMin1r,~]=sepMinMaxHK13([dbl3(1)],hagg,0);
[sepMin2r,~]=sepMinMaxHK13([dbl3(2)],hagg,0);
[sepMin3r,~]=sepMinMaxHK13([dbl3(3)],hagg,0);

sepMinr=[sepMin1r,sepMin2r,sepMin3r];
sepMin(7:9)=sepMinr;

abr3=pi/4*dbl3.^2;

nblowLeft=nbAfterCut3M;
dblowLeft=dbm3;

nbTopMid=nbAfterCut3L;
dbTopMid=dbl3;

nblowRight=nbAfterCut3R;
dblowRight=dbm3;

nbCombo6=[nbl1,nbl2,nbl3,...
            nbm1,nbm2,nbm3];

if sum(nbarsRight)==0
    nbCombo9=zeros(1,9);
else
    nbCombo9=[nbCombo6,nbarsRight];
end

% Left section
Mul=load_conditions(:,2);

[vsepminl,vsepmaxl]=sepMinMaxHK13(dbc(1:3),hagg,1);
if Mul<0
    distrRebarLeftComp=[-redistrRebarM2L];
    listRebarDiamLeftComp=[relistRebarM2L];

    [distrRebarLeftTen,listRebarDiamsLeftTen]=distrRebarRecBeam1DiamLayer...
        (nb3l,dbc(1:3),b,h,brec,hrec,vsepminl);

    distrRebarLeft=[distrRebarLeftTen;
                    distrRebarLeftComp];
    listRebarDiamLeft=[listRebarDiamsLeftTen;
                      listRebarDiamLeftComp];
else

    [distrRebarLeftComp,listRebarDiamLeftComp]=distrRebarRecBeam1DiamLayer...
        (nb3l,dbc(1:3),b,h,brec,hrec,vsepminl);

    distrRebarLeft=[-distrRebarLeftComp;
                    redistrRebarM2L];
    listRebarDiamLeft=[listRebarDiamLeftComp;
                        relistRebarM2L];
end               

[EffLeft,MrLeft,cLeft]=EfRecBeamBars(Mul,fcu,Es,fy,h,b,distrRebarLeft,...
                                     listRebarDiamLeft,hrec);         
% Mid section

[vminSepm,vmaxSepm]=sepMinMaxHK13(dbc(1:3),hagg,1);
Mum=load_conditions(:,3);
distrRebarMidComp=[-redistrRebarL2M];
listRebarDiamMidComp=[relistRebarDiamL2M];

[distrRebarMidTen,listRebarDiamsMidTen]=distrRebarRecBeam1DiamLayer...
    (nb3m,dbc(4:6),b,h,brec,hrec,vminSepm);

distrRebarMid=[distrRebarMidTen;
                distrRebarMidComp];
listRebarDiamMid=[listRebarDiamsMidTen;
                  listRebarDiamMidComp];

[EffMid,MrMid,cMid]=EfRecBeamBars(Mum,fcu,Es,fy,h,b,distrRebarMid,...
                             listRebarDiamMid,hrec);
                         
%% Quantifying cuts

nbcut3sec(1,:)=nb3l-nbAfterCut3L;
nbcut3sec(2,:)=nb3m-nbAfterCut3M;
nbcut3sec(3,:)=nbarsRight-nbAfterCut3L;

%% Anchorage lengths rebars to cut

ld1L=anchorLenBarTen(fcu,fy,h,hrec,db1l);
ld2L=anchorLenBarTen(fcu,fy,h,hrec,db2l);
ld3L=anchorLenBarTen(fcu,fy,h,hrec,db3l);

ldL3=[ld1L,ld2L,ld3L];

ld1m=anchorLenBarTen(fcu,fy,h,hrec,db1m);
ld2m=anchorLenBarTen(fcu,fy,h,hrec,db2m);
ld3m=anchorLenBarTen(fcu,fy,h,hrec,db3m);

ld3M=[ld1m,ld2m,ld3m];

ld1r=anchorLenBarTen(fcu,fy,h,hrec,dbl3(1));
ld2r=anchorLenBarTen(fcu,fy,h,hrec,dbl3(2));
ld3r=anchorLenBarTen(fcu,fy,h,hrec,dbl3(3));

ld3R=[ld1r,ld2r,ld3r];

%% Total rebar volume in beam
[volRebar,lenRebarL,lenRebarM,lenRebarR]=volRebarDesignBeamSpan1DiamLayer(nb3l,...
 nb3m,nbarsRight,abl3,abm3,abr3,ldL3,ld3M,ld3R,nbAfterCut3L,...
 nbAfterCut3M,cutLoc,span);

lenRebar=[lenRebarL;lenRebarM;lenRebarR];

%% Constructability
dbr3=dbl3;

Wunb=Wfac(1:2);
Wnd=Wfac(3);
Wnb=Wfac(4);
Wc=Wfac(5:6);

Wassem=Wfac(7);
Wcutbend=Wfac(8);

dbmin=10;
[sepMinDbmin,~]=sepMinMaxHK13(dbmin,hagg,0);
nbMaxLayer=fix((b-2*brec-2*dvs+sepMinDbmin)/(dbmin+sepMinDbmin));

[FCS1,FCS2,NB,UNB,UND,UC,CS]=CSRebarBeamsRec3DSec(nbMaxLayer,nb3l,...
nb3m,nbarsRight,dbl3,dbm3,dbr3,nbcut3sec(1,:),nbcut3sec(2,:),nbcut3sec(3,:),...
Wunb,Wnd,Wc,Wnb,Wassem,Wcutbend);

CFA=[CS,FCS1,FCS2];

%% Evaluate restrictions

sepRebar(1,:)=sepRebarLeft;
sepRebar(2,:)=sepRebarMid;
sepRebar(3,:)=sepRebarRight;

if all([EffMid<1.0,EffLeft<1.0,EffRight<1.0,Abr>=amin,Abr<=amax,Abm>=amin,...
        Abm<=amax,Abl<=amax,Abl>=amin,isfeasibleRight])
    
    minimumFitness=volRebar;
    bestLenRebar=lenRebar;

    bestnblowLeft=nblowLeft;
    bestdblowLeft=dblowLeft;

    bestnbTopMid=nbTopMid;
    bestdbTopMid=dbTopMid;

    bestnblowRight=nblowRight;
    bestdblowRight=dblowRight;

    bestNbCombo9=nbCombo9;
    bestsepRebar=sepRebar;
    bestnbcut3sec=nbcut3sec;

    bestEffLeft=EffLeft;
    bestEffMid=EffMid;
    bestEffRight=EffRight;

    bestMrLeft=MrLeft;
    bestMrMid=MrMid;
    bestMrRight=MrRight;

    bestRebarDistrMid=distrRebarMid;
    bestListRebarDiamMid=listRebarDiamMid;

    bestDistrRebarLeft=distrRebarLeft;
    bestListRebarDiamLeft=listRebarDiamLeft;

    bestDistrRebarRight=RebarDistrRight;
    bestListRebarDiamRight=ListDiamRight;

    bestcLeft=cLeft;
    bestcMid=cMid;
    bestcRight=cRight;

    bestCFA=CFA;
    constr=0;
else
    constr=0;
    if EffMid>1.0
        constr=constr+1;
    end
    if EffLeft>1.0
        constr=constr+1;
    end
    if EffRight>1.0
        constr=constr+1;
    end
    if Abr<amin
        constr=constr+1;
    end
    if Abr>amax
        constr=constr+1;
    end
    if Abm<amin
        constr=constr+1;
    end
    if Abm>amax
        constr=constr+1;
    end
    if Abl>amax
        constr=constr+1;
    end
    if Abl<amin
        constr=constr+1;
    end
    if isfeasibleRight==false
        constr=constr+1;
    end
    minimumFitness=1e10;
    bestLenRebar=0;
    
    bestNbCombo9=zeros(1,9);

    bestEffLeft=0;
    bestEffMid=0;
    bestEffRight=0;

    bestMrLeft=0;
    bestMrMid=0;
    bestMrRight=0;

    bestsepRebar=zeros(3,3);

    bestRebarDistrMid=[];
    bestListRebarDiamMid=0;
    bestDistrRebarLeft=[];

    bestListRebarDiamLeft=0;
    bestDistrRebarRight=[];
    bestListRebarDiamRight=0;

    bestcLeft=0;
    bestcMid=0;
    bestcRight=0;
    
    bestnblowLeft=zeros(1,3);
    bestdblowLeft=zeros(1,3);

    bestnbTopMid=zeros(1,3);
    bestdbTopMid=zeros(1,3);

    bestnblowRight=zeros(1,3);
    bestdblowRight=zeros(1,3);
    
    bestnbcut3sec=zeros(3,3);
    
    bestCFA=zeros(1,3);
end