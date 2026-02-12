function [bestPerformance,bestEffR,bestMrR,bestcRight,bestnb,bestsepRebar,...
    bestListDiam,bestRebarDistr,isfeasibleRight]=SAPSOBeamsRebarR3DSec...
    (Mur,fc,h,b,hagg,brec,hrec,pmin,pmax,sepRebarleft,distrRebarComp,...
    listRebarDiamComp,nbAfterCut3,nb3l,dbc,Ao3)

%------------------------------------------------------------------------
% Syntax:
% [c_best,bestMr,bestEf,best_area,tbest,h]=PSO3layerBeamsRebar(b,h,duct,...
%    b_rec,h_rec,vSep,fc,Mu,fy,graphConvergencePlot)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: SI - (Kg,cm)
%                  US - (lb,in)
%------------------------------------------------------------------------
% PURPOSE: To determine an optimal reinforcement area for a given beam 
% cross-section with specified initially dimensions (b,h) through the SGD 
% method.
% 
% OUTPUT: c_best,bestMr,bestEf: The neutral axis depth for the optimal 
%                               design, the resistant bending moment for 
%                               the optimal design,
%
%         best_area,tbest:      The optimal reinforcement area, the optimal 
%                               t width of the ISR
%
%         h:                    The final cross-section height in case it 
%                               is modified from the given initial proposal 
%                               value
%
% INPUT:  load_conditions:      vector as [nload,Mu] size: nloads x 2
%
%         factor_fc:            is determined by de applicable design code. 
%                               The ACI 318-19 specifies it as 0.85
%
%         duct:                 is the ductility demand parameter, with 
%                               possible values of 1,2 or 3, for low 
%                               ductility, medium ductility or high 
%                               ductility respectively
%
%         h_rec,b_rec:          is the concrete cover along the height dimension
%                               and the width cross-section dimension, respectively
%                               (cm)
%
%         h,b:                  cross-section dimensions (cm)
%
%         E:                    is the Elasticity Modulus of reinforcing steel
%                               (Kg/mm2)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-07-03
% Copyright (c)  School of Engineering
%                HKUST
%------------------------------------------------------------------------

Es=200e3;
dvs=10;
bp=b-2*brec-2*dvs-max(dbc);

bestPerformance=1e10;
bestEffR=0;

%% Number of rebars
       
[sepMin1,sepMax1]=sepMinMaxHK13(dbc(1),hagg,0);
[sepMin2,sepMax2]=sepMinMaxHK13(dbc(2),hagg,0);
[sepMin3,sepMax3]=sepMinMaxHK13(dbc(3),hagg,0);

sepMin=[sepMin1,sepMin2,sepMin3];
abc=dbc.^2*pi/4;

rr1=b-2*brec-2*dvs;
rr2=b-2*brec-2*dvs;
rr3=b-2*brec-2*dvs;

rr21=sepMin(1);
rr22=sepMin(2);
rr23=sepMin(3);

nbmax1r=fix((rr1)/(rr21));
nbmax2r=fix((rr2)/(rr22));
nbmax3r=fix((rr3)/(rr23));
nbmaxr=[nbmax1r,nbmax2r,nbmax3r];

[nb1,nb2,nb3,isfeasible]=nb3DSec2(nbmaxr,Ao3(3),abc,nb3l,nbAfterCut3) ;

nbr=[nb1,nb2,nb3];
        
%% Rebar separation

sepRebar1=(bp)/(nbr(1)-1);
if nbr(2)>1
    sepRebar2=(bp)/(nbr(2)-1);
else
    sepRebar2=1e5;
end
if nbr(3)>1
    sepRebar3=(bp)/(nbr(3)-1);
else
    sepRebar3=1e5;
end
sepRebarRight=[sepRebar1,sepRebar2,sepRebar3];

for ii=1:3
    if nb3l(ii)-nbAfterCut3(ii)==0 && nbAfterCut3(ii)>2
        sepRebarRight(ii)=sepRebarleft(ii);
    end
end

[Abr,EffR,MrR,cRight,xBest,ListDiam,RebarDistr,isfeasibleRight]=...
DistrEffConstrR3DSec(Mur,fc,Es,h,b,hagg,brec,hrec,pmin,pmax,...
nbAfterCut3,distrRebarComp,listRebarDiamComp,nb3l,dbc,nbr);

if (isfeasibleRight)

    bestPerformance=Abr;
    bestEffR=EffR;
    bestMrR=MrR;
    bestcRight=cRight;
    bestnb=nbr;
    bestsepRebar=sepRebarRight;
    bestListDiam=ListDiam;
    bestRebarDistr=RebarDistr;
end
        
if bestEffR==0
    bestPerformance=1e10;
    bestEffR=0;
    bestMrR=0;
    bestcRight=0;
    bestnb=zeros(1,3);
    bestsepRebar=zeros(1,3);
    bestListDiam=0;
    bestRebarDistr=[];
end
end
