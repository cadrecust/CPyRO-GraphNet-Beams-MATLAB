function cutxLoc=cutLocationSSRecBeam(M,dL)

    %------------------------------------------------------------------------
    % Syntax:
    % cutxLoc=cutLocationSSRecBeam(M,dL)
    %
    %-------------------------------------------------------------------------
    % SYSTEM OF UNITS: Any.
    %
    %------------------------------------------------------------------------
    % PURPOSE: To compute the location of rebar cuts along the length of a
    % beam span
    % 
    % OUTPUT: cutxLoc:  is the vector containing the rebar cut locations as
    %                   longitudinal coordinates [X1,X2,X3,X4]:
    %
    %                   X1 corresponds to the cut location of tension rebar
    %                   for the left cross-section
    %
    %                   X2, X3 corresponds to the cut location of tension rebar
    %                   for the mid cross-section on the left side and
    %                   right side, respectively
    %
    %                   X4 corresponds to the cut location of tension rebar
    %                   for the right cross-section
    %
    % INPUT:  M:        is a vector containing the bending moment
    %                   quantities for each discretization point along the 
    %                   beam, according to dL (described below)
    %
    %         dL:       is the discretization length of the beam for the
    %                   analysis of bending moment distribution along the 
    %                   span of a beam
    %
    %------------------------------------------------------------------------
    % LAST MODIFIED: L.F.Veduzco    2025-02-05
    %                School of Engineering
    %                The Hong Kong University of Science and Technology (HKUST)
    %------------------------------------------------------------------------
    
    ne=length(M(1,:));
    Mleft=M(1,1);
    Mright=M(1,ne);
    [Mmid,mp]=max(M(1,:));
    
    %% Cut location
    if Mleft*Mmid<0 && Mright*Mmid<0 % Standard case
        cutxLoc=[];
        for n=1:ne
            Mref=Mleft;
            if abs(M(1,n))<=abs(0.5*Mref)
                n1=n;
                cutxLoc=[cutxLoc,dL*n];
                break;
            end
        end
    
        Mref=Mmid;
        for n=mp:-1:1
            if abs(M(1,n))<=abs(0.5*Mref)
                cutxLoc=[cutxLoc,dL*n];
                break;
            end
        end
    
        Mref=Mmid;
        for n=mp:1:ne
            if abs(M(1,n))<=abs(0.5*Mref)
                cutxLoc=[cutxLoc,dL*n];
                break;
            end
        end
        
        Mref=Mright;
        for n=ne:-1:1
            if abs(M(1,n))<=abs(0.5*Mref)
                cutxLoc=[cutxLoc,dL*n];
                break;
            end
        end
    elseif Mleft*Mmid<0 && Mright*Mmid>0
        if abs(Mmid)*0.5>abs(Mright) % Case 2.1.a
            cutxLoc=[];
            for n=1:ne
                Mref=Mleft;
                if abs(M(1,n))<=abs(0.5*Mref)
                    n1=n;
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
    
            Mref=Mmid;
            for n=mp:-1:1
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            Mref=Mmid;
            for n=mp:1:ne
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            
            cutxLoc=[cutxLoc,dL*ne];
        elseif abs(Mmid)*0.5<=abs(Mright) % Case 2.2.a
            cutxLoc=[];
            for n=1:ne
                Mref=Mleft;
                if abs(M(1,n))<=abs(0.5*Mref)
                    n1=n;
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
    
            Mref=Mmid;
            for n=mp:-1:1
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            cutxLoc=[cutxLoc,dL*ne];
            cutxLoc=[cutxLoc,dL*ne];
        end
    elseif Mleft*Mmid>0 && Mright*Mmid<0
        if abs(Mmid)*0.5>abs(Mleft) % Case 2.1.b
            cutxLoc=[];
            cutxLoc=[cutxLoc,0];
            Mref=Mmid;
            for n=mp:-1:1
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            
            Mref=Mmid;
            for n=mp:1:ne
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            
            
            Mref=Mright;
            for n=ne:-1:1
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            
        elseif abs(Mmid)*0.5<=abs(Mleft) % Case 2.2.b
            cutxLoc=[]; 
            cutxLoc=[cutxLoc,0];
            cutxLoc=[cutxLoc,0];
            
            Mref=Mmid;
            for n=mp:1:ne
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
            
            
            Mref=Mright;
            for n=ne:-1:1
                if abs(M(1,n))<=abs(0.5*Mref)
                    cutxLoc=[cutxLoc,dL*n];
                    break;
                end
            end
        end
    end
end