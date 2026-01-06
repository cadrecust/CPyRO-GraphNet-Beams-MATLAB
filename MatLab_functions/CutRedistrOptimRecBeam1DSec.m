function [minimumFitness,bestLenRebar,bestsepRebar,bestNbCombo9,bestEffLeft,bestEffMid,...
    bestEffRight,bestMrLeft,bestMrMid,bestMrRight,bestcLeft,bestcMid,bestcRight,...
    bestListRebarDiamLeft,bestListRebarDiamMid,bestListRebarDiamRight,...
    bestDistrRebarLeft,bestRebarDistrMid,bestDistrRebarRight,...
    bestnbcut3sec,bestnblowLeft,bestdblowLeft,bestnbTopMid,bestdbTopMid,...
    bestnblowRight,bestdblowRight,bestCFA,constr]=CutRedistrOptimRecBeam1DSec...
    (load_conditions,fcu,Es,h,b,span,dbc,hagg,brec,hrec,pmin,pmax,sepMin,...
    cutLoc,Wfac,nblm,rebarAvailable)
            
Wunb=Wfac(1:2);
Wnd=Wfac(3);
Wcut=Wfac(4);
Wnb=Wfac(5);
Wcs1=Wfac(6);
Wcs2=Wfac(7);

fy=Es*0.00217;

amin=pmin*b*h;
amax=pmax*b*h;

dvs=10;
bpl=b-2*brec-2*dvs-max(dbc(1));
bpm=b-2*brec-2*dvs-max(dbc(2));

%% Rebar diameters
db1l=dbc(1);
db1m=dbc(2);

%% Rebar area

ab1l=pi/4*db1l^2;
ab1m=pi/4*db1m^2;

nbl1=nblm(1); nbl2=nblm(2); nbl3=nblm(3);
nbm1=nblm(4); nbm2=nblm(5); nbm3=nblm(6);

nb3l=nblm(1:3);
nb3m=nblm(4:6);

Abl=ab1l*nbl1+ab1l*nbl2+ab1l*nbl3;
Abm=ab1m*nbm1+nbm2*ab1m+nbm3*ab1m;

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

%% Cuts and redistribution Left-Mid span
[nbAfterCut3L,nbAfterCut3M,nbAfterCut3R,redistrRebarL2M,relistRebarDiamL2M,...
redistrRebarM2L,relistRebarM2L,redistrRebarM2R,relistRebarM2R]=...
cutRedistRebarSpanRecBeam1DSec(load_conditions,nb3l,nb3m,db1l,db1m,b,h,...
brec,hrec,hagg,fcu,fy);

%% Efficiency analysis
% Optimization of Right section
Mur=load_conditions(:,4);
[Abr,EffRight,MrRight,cRight,nb3r,sepRebarRight,ListDiamRight,...
 RebarDistrRight]=PSOBeamsRebarR1DSec(Mur,fcu,h,b,hagg,brec,hrec,pmin,...
 pmax,sepRebarLeft,redistrRebarM2R,relistRebarM2R,nbAfterCut3L,nb3l,db1l,...
 rebarAvailable);

[sepMin1r,~]=sepMinMaxHK13([db1l],hagg,0);

sepMinr=[sepMin1r];
sepMin(3)=sepMinr;

abr=pi/4*db1l.^2;

nblowLeft=nbAfterCut3M;
dblowLeft=db1m;

nbTopMid=nbAfterCut3L;
dbTopMid=db1l;

nblowRight=nbAfterCut3R;
dblowRight=db1m;

nbCombo6=[nbl1,nbl2,nbl3,...
            nbm1,nbm2,nbm3];

if sum(nb3r)==0
    nbCombo9=zeros(1,9);
else
    nbCombo9=[nbCombo6,nb3r];
end

% Left section
Mul=load_conditions(:,2);

[vsepminl,vsepmaxl]=sepMinMaxHK13(dbc(1),hagg,1);
if Mul<0
    distrRebarLeftComp=[-redistrRebarM2L];
    listRebarDiamLeftComp=[relistRebarM2L];

    [distrRebarLeftTen,listRebarDiamsLeftTen]=distrRebarRecBeam1DSec...
        (nb3l,dbc(1),b,h,brec,hrec,vsepminl);

    distrRebarLeft=[distrRebarLeftTen;
                    distrRebarLeftComp];
    listRebarDiamLeft=[listRebarDiamsLeftTen;
                      listRebarDiamLeftComp];
else

    [distrRebarLeftComp,listRebarDiamLeftComp]=distrRebarRecBeam1DSec...
        (nb3l,dbc(1),b,h,brec,hrec,vsepminl);

    distrRebarLeft=[-distrRebarLeftComp;
                    redistrRebarM2L];
    listRebarDiamLeft=[listRebarDiamLeftComp;
                        relistRebarM2L];
end               

[EffLeft,MrLeft,cLeft]=EfRecBeamBars(Mul,fcu,Es,fy,h,b,distrRebarLeft,...
                                     listRebarDiamLeft,hrec);         
% Mid section

[vminSepm,vmaxSepm]=sepMinMaxHK13(dbc(2),hagg,1);
Mum=load_conditions(:,3);
distrRebarMidComp=[-redistrRebarL2M];
listRebarDiamMidComp=[relistRebarDiamL2M];

[distrRebarMidTen,listRebarDiamsMidTen]=distrRebarRecBeam1DSec...
    (nb3m,dbc(2),b,h,brec,hrec,vminSepm);

distrRebarMid=[distrRebarMidTen;
                distrRebarMidComp];
listRebarDiamMid=[listRebarDiamsMidTen;
                  listRebarDiamMidComp];

[EffMid,MrMid,cMid]=EfRecBeamBars(Mum,fcu,Es,fy,h,b,distrRebarMid,...
                             listRebarDiamMid,hrec);
                         
%% Quantifying cuts

nbcut3sec(1,:)=nb3l-nbAfterCut3L;
nbcut3sec(2,:)=nb3m-nbAfterCut3M;
nbcut3sec(3,:)=nb3r-nbAfterCut3L;

%% Anchorage lengths rebars to cut

ld1L=anchorLenBarTen(fcu,fy,h,hrec,db1l);
ld1M=anchorLenBarTen(fcu,fy,h,hrec,db1m);
ld1R=anchorLenBarTen(fcu,fy,h,hrec,db1l);

%% Total rebar volume in beam
[volRebar,lenRebarL,lenRebarM,lenRebarR]=volRebarDesignBeamSpan1DSec(nb3l,...
 nb3m,nb3r,ab1l,ab1m,abr,ld1L,ld1M,ld1R,nbAfterCut3L,...
 nbAfterCut3M,cutLoc,span);

lenRebar=[lenRebarL;lenRebarM;lenRebarR];

%% Constructability
[FCS1,FCS2,NB,UNB,UND,UC,CS]=CSRebarBeamsRec1DSec(nb3l,nb3m,nb3r,...
    dbl,dbm,dbr,nbcut3l,nbcut3m,nbcut3r,Wunb,Wnd,Wcut,Wnb,Wcs1,Wcs2)
[UNBS,UNDS,UCS,BSS,CFAS,BS,CS]=CFABeamsRec1DSec(nb3l,nb3m,nb3r,db1l,...
db1m,db1l,nbcut3sec(1,:),nbcut3sec(2,:),nbcut3sec(3,:),Wunb,Wnd,Wcut);

%% Rebar distribution restriction
% Left section
[ccl]=rebarDistrConstr3LayerRecBeam1DiamLayer(bpl,nb3l);

% Mid section
[ccm]=rebarDistrConstr3LayerRecBeam1DiamLayer(bpm,nb3m);

% Right section
[ccr]=rebarDistrConstr3LayerRecBeam1DiamLayer(bpm,nbarsRight);

%% Evaluate restrictions

sepRebar(1,:)=sepRebarLeft;
sepRebar(2,:)=sepRebarMid;
sepRebar(3,:)=sepRebarRight;

if all([EffMid<1.0,EffLeft<1.0,EffRight<1.0,Abr>=amin,Abr<=amax,Abm>=amin,...
   Abm<=amax,Abl<=amax,Abl>=amin,ccl==1,ccm==1,ccr==1])
    
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

    bestCFA=CS;
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
    if ccl==0
        constr=constr+1;
    end
    if ccm==0
        constr=constr+1;
    end
    if ccr==0
        constr=constr+1;
    end
    minimumFitness=1e10;
    bestLenRebar=0;
    
    bestNbCombo9=nbCombo9;

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
    bestdblowLeft=0;

    bestnbTopMid=zeros(1,3);
    bestdbTopMid=0;

    bestnblowRight=zeros(1,3);
    bestdblowRight=0;
    
    bestnbcut3sec=nbcut3sec;
    
    bestCFA=0;
end