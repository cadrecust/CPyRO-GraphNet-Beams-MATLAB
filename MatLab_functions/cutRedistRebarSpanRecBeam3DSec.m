function [nbnewl3,nbnewm3,nbnewr3,redistrRebarL2M,relistRebarDiamL2M,...
        redistrRebarM2L,relistRebarM2L,redistrRebarM2R,relistRebarM2R]=...
        cutRedistRebarSpanRecBeam3DSec(load_conditions,nb3l,nb3m,db3l,db3m,...
        b,h,brec,hrec,hagg,fc,fy)

    %------------------------------------------------------------------------
    % Syntax:
    % [nbnewl6,nbnewm6,nbnewr6,redistrRebarL2M,relistRebarDiamL2M,...
    %   redistrRebarM2L,relistRebarM2L,redistrRebarM2R,relistRebarM2R]=...
    %   cutRedistRebarSpanRecBeam(load_conditions,nb6l,nb6m,db6l,db6m,...
    %   b,h,brec,hrec,hagg,fc,fy)
    %
    %-------------------------------------------------------------------------
    % SYSTEM OF UNITS: Any.
    %
    %------------------------------------------------------------------------
    % PURPOSE: To compute the redistribution of cross-section rebar design 
    % coordinates, according to a original rebar distribution and the number
    % cuts to execute. The rebar design must consist of a max of 12 diameter
    % sizes, along the whole length of the beam.
    % 
    % OUTPUT: redistrRebarM2R:    vector (nbars x 2) that contains the new
    %                             rebar cross-section coordinates on the 
    %                             right beam cross-section from the mid
    %                             cross-section, after rebar cutting
    %
    % INPUT:  nb3l:               is the vector containing the number of 
    %                             rebars for each rebar layer component 
    %                             (1 x 3) for the left cross-section of a 
    %                             beam
    %
    %------------------------------------------------------------------------
    % LAST MODIFIED: L.F.Veduzco    2025-02-05
    %                School of Engineering
    %                The Hong Kong University of Science and Technology (HKUST)
    %------------------------------------------------------------------------
    
    [vsepl,~]=sepMinMaxHK13(db3l,hagg,1);
    Mleft=load_conditions(1,2);
    Mmid=load_conditions(1,3);
    Mright=load_conditions(1,4);
    if Mleft*Mmid<0
        % Cuts Left-section
        
        Mtl=0.5*Mleft;
        [nbnew1l,nbnew2l,nbnew3l,Abcutl]=cutRebarSecBeamMr3DSec...
            (nb3l,db3l,Mtl,b,h,hrec,vsepl,fc,fy);
    
        % Re-distribution after cuts of rebar in tension
        % Mid section
        nbnewl3=[nbnew1l,nbnew2l,nbnew3l];
        
        [redistrRebarL2M,relistRebarDiamL2M]=distrRebarRecBeamCuts3DSec...
            (nb3l,db3l,b,h,brec,hrec,vsepl,nbnewl3);
    else
        nbnewl3=nb3l;
        [redistrRebarL2M,relistRebarDiamL2M]=distrRebarRecBeamCuts3DSec...
            (nb3l,db3l,b,h,brec,hrec,vsepl,nbnewl3);
    end
    
    [vsepm,vsepmaxm]=sepMinMaxHK13(db3m,hagg,1);
    if Mmid*Mleft<0
        % Cuts mid-section
        Mtm=0.5*Mmid;
        [nbnew1m,nbnew2m,nbnew3m,Abcutm]=cutRebarSecBeamMr3DSec(nb3m,...
                                     db3m,Mtm,b,h,hrec,vsepm,fc,fy);
        % Re-distribution after cuts of rebar in tension
        % Left section
        nbnewm3=[nbnew1m,nbnew2m,nbnew3m];
        [redistrRebarM2L,relistRebarM2L]=distrRebarRecBeamCuts3DSec...
            (nb3m,db3m,b,h,brec,hrec,vsepm,nbnewm3);
    else
        if abs(Mleft)>0.5*abs(Mmid)
            nbnewm3=nb3m;
    
            [redistrRebarM2L,relistRebarM2L]=distrRebarRecBeamCuts3DSec...
                (nb3m,db3m,b,h,brec,hrec,vsepm,nbnewm3);
        elseif abs(Mleft)<=0.5*abs(Mmid)
            Mtm=0.5*Mmid;
            [nbnew1m,nbnew2m,nbnew3m,Abcutm]=...
                cutRebarSecBeamMr3DSec(nb3m,db3m,Mtm,b,h,hrec,vsepm,fc,fy);
            
            % Re-distribution after cuts of rebar in tension
            % Left section
            nbnewm3=[nbnew1m,nbnew2m,nbnew3m];
            [redistrRebarM2L,relistRebarM2L]=distrRebarRecBeamCuts3DSec...
                (nb3m,db3m,b,h,brec,hrec,vsepm,nbnewm3);
        end
    end
    vsepr=vsepm;
    % Cuts for right-section
    if Mmid*Mright<0
        
        % Re-distribution after cuts of rebar in tension
        % Left section
        nbnewr3=nbnewm3;
        redistrRebarM2R=redistrRebarM2L;
        relistRebarM2R=relistRebarM2L;
    else
        
        if abs(Mright)>0.5*abs(Mmid)
            nbnewr3=nb3m;
    
            [redistrRebarM2R,relistRebarM2R]=distrRebarRecBeamCuts3DSec...
                (nb3m,db3m,b,h,brec,hrec,vsepr,nbnewr3);
        elseif abs(Mright)<=0.5*abs(Mmid)
            Mtm=0.5*Mmid;
            [nbnew1r,nbnew2r,nbnew3r,Abcutr]=cutRebarSecBeamMr3DSec(nb3m,...
                                     db3m,Mtm,b,h,hrec,vsepr,fc,fy);
                                 
            % Re-distribution after cuts of rebar in tension
            % Left section
            nbnewr3=[nbnew1r,nbnew2r,nbnew3r];
            [redistrRebarM2R,relistRebarM2R]=distrRebarRecBeamCuts3DSec...
                (nb3m,db3m,b,h,brec,hrec,vsepr,nbnewr3);
        end
    end
end