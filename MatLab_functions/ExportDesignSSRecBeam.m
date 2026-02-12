function ExportDesignSSRecBeam(directionData,DimBeamsCollec,db9Spans,...
    LenRebarL,LenRebarM,LenRebarR,DistrRebarLeft,DistrRebarMid,DistrRebarRight,...
    totnbSpan,tenbLMRspan,ListRebarDiamLeft,ListRebarDiamMid,ListRebarDiamRight,...
    diamListdSdb,distrSb,beamNSb,ShearBeamDesignCollec)

%------------------------------------------------------------------------
% Syntax:
% ExportResultsBeam(directionData,dim_beams_collection,coordEndBeams,...
%   disposition_rebar_beams3sec,nbarbeamsCollection,arrangemetbarbeams)
%
%-------------------------------------------------------------------------
% SYSTEM OF UNITS: Any.
%
%------------------------------------------------------------------------
% PURPOSE: Computes the exportation of the design results of a beam element
%          into a .txt file on a prescribed folder route.
% 
% INPUT:  directionData:            is the folder disc location to save the
%                                   results
%
%         dim_beams_collection:     is the array containing the cross-section
%                                   dimensions data of the beam element
%
%         coordEndBeams:            is the array containing the centroid 
%                                   coordinates of the initial cross-
%                                   section's end of the beam
%
%         disposition_rebar_beams3sec: is the array containing the local 
%                                      rebar coordinates of each of the 
%                                      three critical design cross-sections
%                                      of the beam
%
%         nbarbeamsCollection:      is the total number of rebars of each 
%                                   of the three design cross-sections, both 
%                                   in tension and compression. 
%                                   Size = [1,6] in format:
%
%               [nbarsLeft_{ten},nbarsLeft_{comp},nbarsCenter_{ten},...
%               nbarsCenter_{comp},nbarsRight_{ten},nbarsRight_{comp}]
%
%         arrangemetbarbeams:       is the list of all the rebar types used 
%                                   in the element
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2025-02-05
% Copyright (c)  School of Engineering
%                HKUST
%------------------------------------------------------------------------

nameFile1='nBarBeamsCollec.csv';
fileid01 = [directionData, nameFile1];
writematrix(tenbLMRspan, fileid01);

nameFile2='dimBeamsCollec.csv';
fileid02 = [directionData, nameFile2];
writematrix(DimBeamsCollec, fileid02);

nameFile3='TendbSec.csv';
fileid03 = [directionData, nameFile3];
writematrix(db9Spans, fileid03);

DistrRebarLMRBeam=[DistrRebarLeft;
                   DistrRebarMid;
                   DistrRebarRight];

nameFile5='DistrBarBeams3sec.csv';
fileid05 = [directionData, nameFile5];
writematrix(DistrRebarLMRBeam, fileid05);

LenRebarLMRBeam=[LenRebarL;
                   LenRebarM;
                   LenRebarR];

nameFile12='LenRebarLMRBeams.csv';
fileid12 = [directionData, nameFile12];
writematrix(LenRebarLMRBeam, fileid12);

nameFile11='NbarsTot3sec.csv';
fileid11 = [directionData, nameFile11];
writematrix(totnbSpan, fileid11);

ListDiamLMRBeam=[ListRebarDiamLeft;
                 ListRebarDiamMid;
                 ListRebarDiamRight];

nameFile6='ListDiamLMRBeams.csv';
fileid06 = [directionData, nameFile6];
writematrix(ListDiamLMRBeam, fileid06);

nameFile7='DistrSb.csv';
fileid07 = [directionData, nameFile7];
writematrix(distrSb, fileid07);

nameFile8='beamNSb.csv';
fileid08 = [directionData, nameFile8];
writematrix(beamNSb, fileid08);

nameFile9='DiamListSb.csv';
fileid09 = [directionData, nameFile9];
writematrix(diamListdSdb, fileid09);

% Exporting shear design data (if required)
nameFile10='ShearBeamDesign.csv';
fileid10 = [directionData, nameFile10];
writematrix(ShearBeamDesignCollec, fileid10);
