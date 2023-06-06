clc
clear all

% Pushover_3D_Frame_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the non-linear static Pushover analysis for a
%    reinforced concrete 3D frame
%
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

nnodes=8;
nbars=8;

%% Materials
% f'c of each element
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element in function of f'c
E=14000.*sqrt(fpc);

v=[0.2; % Poisson modulus
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2];

G=E./(2.*(1+v)); % Shear modulus

%% Geometry/Topology
% cross-section dimensions of each element (rectangular geometry)
dimensions=[40 40;
            25 50;
            40 40;
            25 50;
            25 50;
            40 40;
            25 50;
            40 40];
        
% cross-section area of each element
A=dimensions(:,1).*dimensions(:,2);
        
% Cross-section inertia
Iy=1/12.*dimensions(:,1).*dimensions(:,2).^3;
Iz=1/12.*dimensions(:,2).*dimensions(:,1).^3;

adim=dimensions(:,2).*0.5;
bdim=dimensions(:,1).*0.5;

% Saint Venant constant (polar inertia - Torsion)
J=adim.*bdim.^3.*(16/3-3.36.*bdim./adim.*(1-bdim.^4./(12.*adim.^4)));
      
%% Topology and node coordinates

% OPTION 1: Manually given
% Coordinates of each node
coordxyz=[0 0 0;
          0 0 300;
          0 500 300;
          0 500 0;
          400 0 0;
          400 0 300;
          400 500 300;
          400 500 0];
      
% Connectivity
NiNf=[1 2;
    2 3;
    3 4;
    3 7;
    2 6;
    5 6;
    6 7;
    7 8];

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
    25 0;
    26 0;
    27 0;
    28 0;
    29 0;
    30 0;
    43 0;
    44 0;
    45 0;
    46 0;
    47 0;
    48 0];
       
%% Additional data (optional)
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Beam";
           6 "Col";
           7 "Beam";
           8 "Col"];
       
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
           8 "Fixed" "Fixed"];
       
%% Local z axis of each element
eobars=[0 1 0;
        0 0 1;
        0 1 0;
        0 0 1;
        0 0 1;
        0 1 0;
        0 0 1;
        0 1 0];
    
%% Loads       
beams_LL=[1 -40; % Uniformly distributed loads over the beams
          2 -40;
          3 -40;
          4 -40];

% Assignation of distributed loads on beams
qbarz=zeros(nbars,2);
for i=1:beams
    qbarz(elembeams(i),2)=beams_LL(i,2);
end
   
% Lateral equivalent seismic forces from a modal analysis
seismicForces=[1000; 
               1000];
            
% Degrees of freedom over which each seismic force is applied (one for
% each seismic force)
dofSeismicForces=[7 13];

%% Plastic moments of each element's ends
Mp=[11680000 11680000;
    4490000 4490000;
    12363000 12363000;
    14900000 14900000;
    16800000 16800000;
    12976940 12976940;
    3976940 3976940;
    11490000 11490000]; % Kg-cm

% Height of each floor
hfloor=[300];  
nfloors=length(hfloor);

nbays=2;

%% PUSHOVER IN POSITIVE DIRECTION OF FORCES

[lambdaRight,pdriftDIRight,driftDIRight,defBasedDIRight,maxDispRight,...
 barPlasNodeRight]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
 NiNf(:,1),NiNf(:,2),eobars,supports,bc,seismicForces,hfloor,...
 dofSeismicForces,nbays,0.01,0.2,4);

%% PUSHOVER IN NEGATIVE DIRECTION OF FORCES

seismicForces=-seismicForces;
    
[lambdaLeft,pdriftDILeft,driftDILeft,defBasedDILeft,maxDispLeft,...
 barPlasNodeLeft]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
 NiNf(:,1),NiNf(:,2),eobars,supports,bc,seismicForces,hfloor,...
 dofSeismicForces,nbays,0.01,0.2,4);

%% Final results
SafetyFac=min([max(lambdaRight), max(lambdaLeft)])
pdriftDI=min([sum(pdriftDIRight)/nfloors,sum(pdriftDILeft)/nfloors])
driftDI=min([sum(driftDIRight)/nfloors,sum(driftDILeft)/nfloors])
dbDI=min([sum(defBasedDIRight)/nfloors,sum(defBasedDILeft)/nfloors])

Max_Displacement=max(max(maxDispLeft),max(maxDispRight))
