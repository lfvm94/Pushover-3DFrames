clc
clear all

% Pushover_3D_Frame_Ex02
%----------------------------------------------------------------
% PURPOSE 
%    To compute the non-linear static Pushover analysis for a
%    reinforced concrete 3D frame by importing the topology from 
%    SAP2000 model.
%
%----------------------------------------------------------------
%    Notes: 
%           To download the SM Toolbox visit: 
%           https://github.com/RJ-Soft/SM-Toolbox/releases/tag/7.0.2
%----------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

%% Topology and node coordinates

% OOPTION 2: Exported from a SAP2000 model (the SM Toolbox is required)
APIDLLPath ='C:\Program Files\Computers and Structures\SAP2000 22\SAP2000v1.dll';
ProgramPath ='C:\Program Files\Computers and Structures\SAP2000 22\SAP2000.exe';

ModelName = 'Frame_Ex01.sdb';
ModelPath = fullfile('C:\Users\luizv\OneDrive\Pushover_3DFrames\FrameSAP2000_Ex02',ModelName);

[coordxyz,NiNf]=ExtractTopologySAP2000(ProgramPath,APIDLLPath,...
                                            ModelPath);
                                        
coordxyz=coordxyz*2.54; % to change the location coordinates from in to cm

% OPTION 1: Coordinates and connectivity manually given
%{
NiNf=[1,2;2,3;3,4;3,5;5,6;7,8;8,9;9,10;9,11;11,12;8,2;9,3;11,5;2,13;13,14;14,3;14,15;15,5;8,16;16,17;17,9;17,18;18,11;16,13;17,14;18,15;13,19;19,20;20,14;20,21;21,15;16,22;22,23;23,17;23,24;24,18;22,19;23,20;24,21];
coordxyz=[0,0,0;0,0,118.110236220472;137.795275590551,0,118.110236220472;137.795275590551,0,0;393.700787401575,0,118.110236220472;393.700787401575,0,0;0,157.480314960630,0;0,157.480314960630,118.110236220472;137.795275590551,157.480314960630,118.110236220472;137.795275590551,157.480314960630,0;393.700787401575,157.480314960630,118.110236220472;393.700787401575,157.480314960630,0;0,0,236.220472440945;137.795275590551,0,236.220472440945;393.700787401575,0,236.220472440945;0,157.480314960630,236.220472440945;137.795275590551,157.480314960630,236.220472440945;393.700787401575,157.480314960630,236.220472440945;0,0,354.330708661417;137.795275590551,0,354.330708661417;393.700787401575,0,354.330708661417;0,157.480314960630,354.330708661417;137.795275590551,157.480314960630,354.330708661417;393.700787401575,157.480314960630,354.330708661417]*2.54;
%}

nbars=39;

%% Materials
% f'c of each element
fc1=300;
fpc=zeros(nbars,1)+fc1;

% Elasticity modulus of each element in function of f'c
E=14000.*sqrt(fpc);

v1=0.2;
v=zeros(nbars,1)+v1; % Poisson modulus

G=E./(2.*(1+v)); % Shear modulus

%% Geometry/Topology
% cross-section dimensions of each element (rectangular geometry)
dimensions=[30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40];

% cross-section area of each element
A=dimensions(:,1).*dimensions(:,2);

% Cross-section inertia
Iy=1/12.*dimensions(:,1).*dimensions(:,2).^3;
Iz=1/12.*dimensions(:,2).*dimensions(:,1).^3;

adim=dimensions(:,2).*0.5;
bdim=dimensions(:,1).*0.5;

% Saint Venant constant (polar inertia - Torsion)
J=adim.*bdim.^3.*(16/3-3.36.*bdim./adim.*(1-bdim.^4./(12.*adim.^4)));
      
%% Prescribed boudnary conditions [dof, displacement]
bc=[1 0;
    2 0;
    3 0;
    4 0;
    5 0;
    6 0;
    19 0;
    20 0;
    21 0;
    22 0;
    23 0;
    24 0;
    31 0;
    32 0;
    33 0;
    34 0;
    35 0;
    36 0;
    37 0;
    38 0;
    39 0;
    40 0;
    41 0;
    42 0;
    55 0;
    56 0;
    57 0;
    58 0;
    59 0;
    60 0;
    67 0;
    68 0;
    69 0;
    70 0;
    71 0;
    72 0];
       
%% Additional data (optional)
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Col";
           6 "Col";
           7 "Beam";
           8 "Col";
           9 "Beam";
           10 "Col";
           11 "Beam";
           12 "Beam";
           13 "Beam";
           14 "Col";
           15 "Beam";
           16 "Col";
           17 "Beam";
           18 "Col";
           19 "Col";
           20 "Beam";
           21 "Col";
           22 "Beam";
           23 "Col";
           24 "Beam";
           25 "Beam";
           26 "Beam";
           27 "Col";
           28 "Beam";
           29 "Col";
           30 "Beam";
           31 "Col";
           32 "Col";
           33 "Beam";
           34 "Col";
           35 "Beam";
           36 "Col";
           37 "Beam";
           38 "Beam";
           39 "Beam";];
       
elemcols=[];
elembeams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elembeams=[elembeams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elemcols=[elemcols,j];
    end
end

%% Elements' end conditions
supports=[1 "Fixed" "Fixed";
           2 "Fixed" "Fixed";
           3 "Fixed" "Fixed";
           4 "Fixed" "Fixed";
           5 "Fixed" "Fixed";
           6 "Fixed" "Fixed";
           7 "Fixed" "Fixed";
           8 "Fixed" "Fixed";
           9 "Fixed" "Fixed";
           10 "Fixed" "Fixed";
           11 "Fixed" "Fixed";
           12 "Fixed" "Fixed";
           13 "Fixed" "Fixed";
           14 "Fixed" "Fixed";
           15 "Fixed" "Fixed";
           16 "Fixed" "Fixed";
           17 "Fixed" "Fixed";
           18 "Fixed" "Fixed";
           19 "Fixed" "Fixed";
           20 "Fixed" "Fixed";
           21 "Fixed" "Fixed";
           22 "Fixed" "Fixed";
           23 "Fixed" "Fixed";
           24 "Fixed" "Fixed";
           25 "Fixed" "Fixed";
           26 "Fixed" "Fixed";
           27 "Fixed" "Fixed";
           28 "Fixed" "Fixed";
           29 "Fixed" "Fixed";
           30 "Fixed" "Fixed";
           31 "Fixed" "Fixed";
           32 "Fixed" "Fixed";
           33 "Fixed" "Fixed";
           34 "Fixed" "Fixed";
           35 "Fixed" "Fixed";
           36 "Fixed" "Fixed";
           37 "Fixed" "Fixed";
           38 "Fixed" "Fixed";
           39 "Fixed" "Fixed"];
       
%% Local z axis of each element
for i=1:nbars
    if type_elem(i,2)=="Col"
        eobars(i,:)=[0 1 0];
    else
        eobars(i,:)=[0 0 1];
    end
end
    
%% Loads       
beams_LL=-70; % Uniformly distributed loads over the beams

% Assignation of distributed loads on beams
qbarz=zeros(nbars,2);
qbarz(elembeams,2)=beams_LL;

% Lateral equivalent seismic forces from a modal analysis
nbays=2;
seismicForces=[1500; 
               1500;
               2000;
               2000;
               2500;
               2500];

% Degrees of freedom over which each seismic force is applied (one for
% each seismic force)
nodeForces=[2 8 13 16 19 22];
dofSeismicForces=nodeForces*6-5; % x direction

%% Plastic moments of each element's ends
Mp=[12080000 12080000;
    6490000 6490000;
    13063000 13076940;
    5490000 5490000;
    12080000 12080000;
    11863000 1876940;
    7363000 7976940;
    13090000 13090000;
    8680000 8680000;
    12090000 12090000;
    9363000 9976940;
    5490000 5490000;
    8680000 8680000; 
    12080000 12080000;
    6490000 6490000;
    13063000 13076940;
    5490000 5490000;
    12080000 12080000;
    11863000 1876940;
    7363000 7976940;
    13090000 13090000;
    8680000 8680000;
    12090000 12090000;
    9363000 9976940;
    5490000 5490000;
    8680000 8680000;
    12080000 12080000;
    6490000 6490000;
    13063000 13076940;
    5490000 5490000;
    12080000 12080000;
    11863000 1876940;
    7363000 7976940;
    13090000 13090000;
    8680000 8680000;
    12090000 12090000;
    9363000 9976940;
    5490000 5490000;
    8680000 8680000]; % Kg-cm

% Height of each floor
hfloor=[300;300;300];  
nfloors=length(hfloor);

%% PUSHOVER IN POSITIVE DIRECTION OF FORCES

[lambdaRight,pdriftDIRight,driftDIRight,defBasedDIRight,maxDispRight,...
 barPlasNodeRight]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
 NiNf(:,1),NiNf(:,2),eobars,supports,bc,seismicForces,hfloor,...
 dofSeismicForces,nbays,0.01,0.2,3);

%% PUSHOVER IN NEGATIVE DIRECTION OF FORCES

seismicForces=-seismicForces;
    
[lambdaLeft,pdriftDILeft,driftDILeft,defBasedDILeft,maxDispLeft,...
 barPlasNodeLeft]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
 NiNf(:,1),NiNf(:,2),eobars,supports,bc,seismicForces,hfloor,...
 dofSeismicForces,nbays,0.01,0.2,3);

%% Final results
SafetyFac=min([max(lambdaRight), max(lambdaLeft)])
pdriftDI=min([sum(pdriftDIRight)/nfloors,sum(pdriftDILeft)/nfloors])
driftDI=min([sum(driftDIRight)/nfloors,sum(driftDILeft)/nfloors])
dbDI=min([sum(defBasedDIRight)/nfloors,sum(defBasedDILeft)/nfloors])

Max_Displacement=max(max(maxDispLeft),max(maxDispRight))
