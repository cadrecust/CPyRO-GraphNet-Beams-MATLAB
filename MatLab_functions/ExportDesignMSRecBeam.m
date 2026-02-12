function ExportDesignMSRecBeam(directionData,nMSbeams,eleMSBeams,DimMSBeams,...
    dbMSBeams,cutLocMSBeam,DistrRebarLeft,DistrRebarMid,DistrRebarRight,...
    NTotalMSBeams,tenbLMRMSBeam,ListRebarDiamLeft,ListRebarDiamMid,...
    ListRebarDiamRight,LenRebarLeftMSBeams,LenRebarMidMSBeams,...
    LenRebarRightMSBeams,diamListdSdb,distrSb,beamNSb,ShearBeamDesignCollec)

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

%% Max number of spans in all Multi-span beams

for i=1:nMSbeams
    nbeamsMS(i)=length(eleMSBeams{i});
    nbarsTenBeam(i)=sum(tenbLMRMSBeam{i});
    nbTotBeam(i)=sum(NTotalMSBeams{i});
end

nbeamsMax=max(nbeamsMS);
elemMSBeams=zeros(nMSbeams,nbeamsMax);
for i=1:nMSbeams
    elemMSBeams(i,1:nbeamsMS(i))=eleMSBeams{i};
end
nameFile0='eleMSBeams.csv';
fileid00 = [directionData, nameFile0];
writematrix(elemMSBeams, fileid00);

maxBarsTenMSBeam=max(nbarsTenBeam);
maxBarsTotalMSBeam=max(nbTotBeam);
nsbMaxBeam=max(beamNSb);

%% Number of rebars per span, per multi-span
tenbLMRspan=zeros(nMSbeams,9*nbeamsMax);
for i=1:nMSbeams
    tenbLMRspan(i,1:9*nbeamsMS(i))=tenbLMRMSBeam{i};
end

nameFile1='nBarBeamsCollec.csv';
fileid01 = [directionData, nameFile1];
writematrix(tenbLMRspan, fileid01);

%% Dimensions
DimBeams=zeros(nMSbeams,nbeamsMax*3+2);
for i=1:nMSbeams
    DimBeams(i,1:nbeamsMS(i)*3+2)=DimMSBeams{i};
end

nameFile2='dimBeamsCollec.csv';
fileid02 = [directionData, nameFile2];
writematrix(DimBeams, fileid02);

%% Rebar diameters per cross-section, per beam, per multi-span beam
dbLMRMSBeams=zeros(nMSbeams,9*nbeamsMax);
for i=1:nMSbeams
    dbLMRMSBeams(i,1:9*nbeamsMS(i))=dbMSBeams{i};
end

nameFile3='DiamBarBeamsCollec3sec.csv';
fileid03 = [directionData, nameFile3];
writematrix(dbLMRMSBeams, fileid03);

%% Locations of rebar cuts
cutLocMSBeams=zeros(nMSbeams,4*nbeamsMax);
for i=1:nMSbeams
    cutLocMSBeams(i,1:4*nbeamsMS(i))=cutLocMSBeam{i};
end
nameFile4='CutLocationRebarMSBeams.csv';
fileid04 = [directionData, nameFile4];
writematrix(cutLocMSBeams, fileid04);

%% List of rebar diameter sizes

ListDiamLMRMSBeam=zeros(maxBarsTotalMSBeam,nMSbeams);
for i=1:nMSbeams
    ListDiamLMRMSBeam(1:nbTotBeam(i),i)=[ListRebarDiamLeft{i};
                                        ListRebarDiamMid{i};
                                        ListRebarDiamRight{i}];
end

nameFile5='ListDiamLMRMSBeams.csv';
fileid05 = [directionData, nameFile5];
writematrix(ListDiamLMRMSBeam, fileid05);

%% Distribution of rebars per cross-section per multi-span
DistrRebarLMRMSBeam=zeros(maxBarsTotalMSBeam,2*nMSbeams);
for i=1:nMSbeams
    DistrRebarLMRMSBeam(1:nbTotBeam(i),(i-1)*2+1:i*2)=[-DistrRebarLeft{i};
                                                        DistrRebarMid{i};
                                                       -DistrRebarRight{i}];
end

nameFile6='DistrBarBeams3sec.csv';
fileid06 = [directionData, nameFile6];
writematrix(DistrRebarLMRMSBeam, fileid06);

%% Side rebars - Distribution MS Beams
DistrSbMSBeam=zeros(nsbMaxBeam,2*nMSbeams);
for i=1:nMSbeams
    DistrSbMSBeam(1:beamNSb(i),(i-1)*2+1:i*2)=distrSb{i};
end

nameFile7='DistrSb.csv';
fileid07 = [directionData, nameFile7];
writematrix(DistrSbMSBeam, fileid07);

%% Number of side bars per MS - Beam
nameFile8='beamNSb.csv';
fileid08 = [directionData, nameFile8];
writematrix(beamNSb, fileid08);

%% Side rebars - List of diameter sizes MS Beams
ListDiamSbMSBeam=zeros(nsbMaxBeam,nMSbeams);
for i=1:nMSbeams
    ListDiamSbMSBeam(1:beamNSb(i),i)=diamListdSdb{i};
end

nameFile9='DiamListSb.csv';
fileid09 = [directionData, nameFile9];
writematrix(ListDiamSbMSBeam, fileid09);

%% Shear reinforcement - MS Beams
shearMSBeams=zeros(nMSbeams,6*nbeamsMax);
for i=1:nMSbeams
    shearMSBeams(i,1:6*nbeamsMS(i))=ShearBeamDesignCollec{i};
end

nameFile10='ShearBeamDesign.csv';
fileid10 = [directionData, nameFile10];
writematrix(shearMSBeams, fileid10);

%% Number of beams per MS Beam
nameFile11='nBeamsMSBeams.csv';
fileid11 = [directionData, nameFile11];
writematrix(nbeamsMS, fileid11);

%% Total number of rebar per beam
NTotalSpansMSBeams=zeros(nMSbeams,nbeamsMax*3);
for i=1:nMSbeams
    nspans=nbeamsMS(i);
    NTotalSpansMSBeams(i,1:nspans*3)=NTotalMSBeams{i};
end

nameFile11='totalNbMSBeams.csv';
fileid11 = [directionData, nameFile11];
writematrix(NTotalSpansMSBeams, fileid11);

%% Length of longitudinal rebar
LenRebarLMRMSBeam=zeros(maxBarsTenMSBeam,nMSbeams);
for i=1:nMSbeams
    LenRebarLMRMSBeam(1:nbarsTenBeam(i),i)=[LenRebarLeftMSBeams{i};
                                         LenRebarMidMSBeams{i};
                                         LenRebarRightMSBeams{i}];
end
nameFile12='LenRebarLMRBeams.csv';
fileid12 = [directionData, nameFile12];
writematrix(LenRebarLMRMSBeam, fileid12);

% ------------------------------------------------------------------------