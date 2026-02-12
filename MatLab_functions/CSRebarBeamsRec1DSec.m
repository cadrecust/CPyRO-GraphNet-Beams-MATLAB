function [FCS1,FCS2,NB,UNB,UND,CSDCut,CS]=CSRebarBeamsRec1DSec(nbmaxLay,nb3l,...
    nb3m,nb3r,dbl,dbm,dbr,nbcut3l,nbcut3m,nbcut3r,Wunb,Wnd,Wcut,Wnb,Wcs1,Wcs2)

    %------------------------------------------------------------------------
    % Syntax:
    % [FCS1,FCS2,NB,UNB,UND,UC,CS]=CSRebarBeamsRec1DSec(nb3l,nb3m,nb3r,...
    %  dbl,dbm,dbr,nbcut3l,nbcut3m,nbcut3r,Wunb,Wnd,Wcut,Wcs1,Wcs2)
    %-------------------------------------------------------------------------
    % SYSTEM OF UNITS: Any.
    %
    %------------------------------------------------------------------------
    % PURPOSE: To compute the constructability score of a rebar design in a
    % concrete beam, consisting of a max of 1 rebar diameter size per 
    % cross-section
    % 
    % OUTPUT: CS:       Constructability Score
    %
    % INPUT:  Wunb:     is the weight factor for the Uniformity of Number
    %                   of Rebars (rebar distribution)
    %
    %         Wnd:      is the weight factor for the diversity of number of
    %                   rebar diameter sizes
    %
    %         Wcut:     is the weight factor for the number of cuts
    %
    %         Wnb:      is the weight factor for the number of rebars
    %
    %         Wcs1:     is the weight factor for the constructability score
    %                   of rebar assembly
    %
    %         Wcs2:     is the weight factor for the constructability score
    %                   of rebar cutting and bending
    %
    %------------------------------------------------------------------------
    % LAST MODIFIED: L.F.Veduzco    2025-02-05
    %                School of Engineering
    %                The Hong Kong University of Science and Technology (HKUST)
    %------------------------------------------------------------------------
        
    %% Number of rebars
    NB1 = 1-min(sum(nb3l)/(3*nbmaxLay),1);
    NB2 = 1-min(sum(nb3m)/(3*nbmaxLay),1);
    NB3 = 1-min(sum(nb3r)/(3*nbmaxLay),1);
    
    NB = (1/3*(NB1 + NB2 + NB3)).^Wnb(1) ;
    
    %% Uniformity of Number of Rebars along the beam's span
    % Right section
    nb3sec=[nb3l;
            nb3m;
            nb3r];

    nbcuts3sec=[nbcut3l;
                nbcut3m;
                nbcut3r];
    db3sec=[dbl;
            dbm;
            dbr];

    for i=1:3
        nb3s=nb3sec(i,:);
        nlays=sum((nb3s~=0));
        nlays=ceil(nlays);
        UNBL=zeros(nlays-1,1);
        for j=1:nlays-1
            if nb3s(j+1)<=(nb3s(1))
                UNBL(j)=(nb3s(j+1)/(nb3s(1)))^Wunb(2);
            else
                UNBL(j)=1;
            end
        end
        if nlays>1
            UNBS(i)=1/(nlays-1)^Wunb(1)*sum(UNBL);
        else
            UNBS(i)=1;
        end
    end

    UNB=(sum(UNBS)/3)^Wunb(2);
    
    %% Number of diameters sizes

    dbsec=db3sec(:,1);
    nb3s=nb3sec(:,1);

    % Quantifying the number of different diameter sizes per section
    NDS = difDiamSizesLayers(dbsec',nb3s');

    UND=(1/NDS)^Wnd(1);
    
    %% Cutting and bending
    % Number of cuts and number of layers
    for i=1:3
        nb3s=nb3sec(i,:);
        nbcuts=sum(nbcuts3sec(i,:));
        nlays=sum((nb3s~=0));
        nlays=ceil(nlays);
        NC=(nbcuts/sum(nb3s))^Wcut(2);
        
        if nbcuts>0
            UCS(i)=(1/nlays)^Wcut(1)*NC;
        else
            UCS(i)=1;
        end
    end
    NC=sum(UCS)/3;
    
    % Irregularity of cut lengths
    
    CSDCut=(1/NDS)^0.5;
    
    %% Constructability of assembly and placing
    FCS1 = ((UNB + UND + NB)/3 ) ;
    
    %% Constructability of cutting and bending
    FCS2 = ( NC + CSDCut )/2 ;
    
    %% Constructability Score
    CS = 1/2*(FCS1^Wcs1 + FCS2^Wcs2);
end

%% Function appendix
function nDiams = difDiamSizesLayers(vector,refDifZero)
    nitemsVec=length(vector);
    dbsecdif0=[];
    for ii=1:nitemsVec
        if refDifZero(ii)~=0
            dbsecdif0=[dbsecdif0,vector(1,ii)];
        end
    end
    [dbsort,ind]=sort(dbsecdif0);
    ndbdif0=length(dbsort);
    nbj=1; 
    jj=1;
    a=[];
    for ii=1:ndbdif0-1
        if dbsort(ii)~=dbsort(ii+1)
            a=[a,nbj]; 
            nbj=1; 
            jj=jj+1;
        else
            nbj=nbj+1;
        end
    end
    a=[a,nbj];
    nDiams=length(a);
end