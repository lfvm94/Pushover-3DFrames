function [historyIncLoad,pdriftDI,driftDI,defBasedDI,maxDisplacement,...
         barPlasNode]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
         ni,nf,eobars,support,bc,seismicforces,Hfloor,dofForces,nbays,...
         dload,dydu,NmaxPlasteps)

%------------------------------------------------------------------------
% [historyIncLoad,pdriftDI,driftDI,defBasedDI,maxDisplacement,...
%  barPlasNode]=Pushover3DFrames(qbarz,A,Mp,E,G,Iy,Iz,J,coordxyz,...
%  ni,nf,eobars,support,bc,seismicforces,Hfloor,dofForces,nbays,...
%  dload,dydu,NmaxPlasteps)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute a static non-linear pushover analysis of a 3D frame
%  
% 
% INPUT:  A = [area_bar;
%               ...]                 area of all elements
%
%         Mp = [Mpi Mpj;             Plastic Moment for each member 
%               ... ]                (i) initial node, (j) final node
%
%         E = [e_bar;                Elasticity modulus of each element
%               ...]                    
%
%         G = [g_bar;                Shear modulus of elasticity of each
%               ...]                    
%
%         v = [v-bar;                Poisson modulus for each bar
%             ...]
%
%         Iy = [inertia_bar;         momentum of inertia for all elements'
%                       ...]         cross-section with respect to their
%                                    local Y axis (see doc.)
%         
%         Iz = [inertia_bar;         momentum of inertia for all elements'
%                       ...]         cross-section with respect to their
%                                    local Z axis (see doc.)
%
%         J = [polar-inertia_bar;    polar momentum of inertia for all
%                       ...]         elements' cross-sections
%
%         coordxyz = [x,y,z;         node coordinates for all nodes
%                       ...];
%
%         ni                         list of initial nodes of all bars,
%         nf                         list of final nodes of all bars:
%                                         size = [nbars,1]
% 
%         eobars                     local Z axis of each element in the
%                                    global system of reference
%
%         qbarz = [bar, load;
%                   ..    .. ]       uniform distributed loads on elements
%
%         support = [i, j]           support at each bar's end
%                                    options: "Art" or "Fixed"
%                                    (i) initial node, (j) final node
%
%         bc                         restricted dof
%
%         seismicForces = [f(1);]    lateral forces per floor:
%                          f(n);]    size = [nfloors,1]
%
%         Hfloor = [h(1);            Height of each floor from bottom
%                    h(n)]           to top: size = [nfloors,1]
%
%         dofForces = [dof-f(1),     dof at which the lateral forces are
%                       dof-f(n)]    applied (from bottom to top) - global
%
%
%         dydu:                      max ratio between the elastic floor
%                                    deformation and the plastic floor
%                                    deformation -> du/dy -> (dmax-dy)/dy 
%
%         NmaxPlasteps:              Max number of plastic formation steps
%
% OUTPUT: historyIncLoad             history of incremental load factors at
%                                    at which plastic moments are reached
%
%         pdriftDI                   Plastic inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         driftDI                    Inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         defBasedDI                 Deformation based Damage Index
%                                    per floor size = [nfloors,1]
%
%         maxDisplacement            Max absoloute lateral displacement
%                                    for each floor: size = [nfloors,1]
% 
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nbars=length(E);
nnodes=length(coordxyz(:,1));

% Topology matrix 
Edof=zeros(nbars,13);
for i=1:nbars
    Edof(i,1)=i;
    
    Edof(i,2)=ni(i)*6-5;
    Edof(i,3)=ni(i)*6-4;
    Edof(i,4)=ni(i)*6-3;
    Edof(i,5)=ni(i)*6-2;
    Edof(i,6)=ni(i)*6-1;
    Edof(i,7)=ni(i)*6;
    
    Edof(i,8)=nf(i)*6-5;
    Edof(i,9)=nf(i)*6-4;
    Edof(i,10)=nf(i)*6-3;
    Edof(i,11)=nf(i)*6-2;
    Edof(i,12)=nf(i)*6-1;
    Edof(i,13)=nf(i)*6;
    
end

[ndof,edof]=nonRestrcDof(nnodes,bc,6);

mbar=zeros(nbars,2); % to save the plastic moments at each articulation
                     % of each bar as a plastification occurs
nfloors=length(Hfloor); % will be used to compute the relative floor
                         % displacements and damage indices of each floor

dispHistFloorLeft=[];
forceFloorLeft=[];

dispHistFloorRight=[];
forceFloorRight=[];

zero_disp=zeros(nfloors,1);
dispHistFloorLeft=[dispHistFloorLeft,zero_disp];
forceFloorLeft=[forceFloorLeft,zero_disp];

dispHistFloorRight=[dispHistFloorRight,zero_disp];
forceFloorRight=[forceFloorRight,zero_disp];

plast_bars=zeros(2,nbars);
barPlasNode=[];
iter_collection=[];
historyIncLoad=[];
looping=0;
incLoad=1.0;
iteration=0;
while looping==0
    Kglobal=zeros(6*nnodes);
    fglobal=zeros(6*nnodes,1);
    fglobal(dofForces)=seismicforces*incLoad;
        
    elmmat=zeros(12*nbars,12);

    for i=1:nbars      
        
        ex=[coordxyz(ni(i),1) coordxyz(nf(i),1)];
        ey=[coordxyz(ni(i),2) coordxyz(nf(i),2)];
        ez=[coordxyz(ni(i),3) coordxyz(nf(i),3)];
     
        Ex(i,:)=ex;
        Ey(i,:)=ey;
        Ez(i,:)=ez;
        
        eo=eobars(i,:);
        ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
         
        eq=[0 0 qbarz(i,2) 0];

        [Kebar,febar]=beam3e(ex,ey,ez,eo,ep,eq); % This is a CALFEM
                                                 % function
                                                 % Download at: 
                                                 % https://www.byggmek.lth.se/english/calfem/

        if support(i,2)=="Fixed" && support(i,3)=="Art" % DOF: M3
            Mpl=mbar(i,2);
            eq=[0 0 qbarz(i,2) 0 Mpl];
            
            [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,1);
             
         elseif support(i,2)=="Art" && support(i,3)=="Fixed" % DOF: M3

             Mpl=mbar(i,1);
             eq=[0 0 qbarz(i,2) 0 Mpl];
             [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,2);
         elseif support(i,2)=="Art" && support(i,3)=="Art" % DOF: M3

             Mpl=mbar(i,1);
             Mp2=mbar(i,2);
             eq=[0 0 qbarz(i,2) 0 [Mpl,Mp2]];
             [Kebar,febar]=beamArt3e(ex,ey,ez,ep,eq,eo,3);
         end
         elmmat((i-1)*12+1:12*i,:)=Kebar; % storing Kebar
         fbe(:,i)=febar; % storing element forces for further use
         
         % Assembling global stiffness matrix
         [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);

    end 
    globalKreduced=Kglobal(edof,edof');
    if iteration==0
        det_Kred=det(globalKreduced);
    else
        det_Kred_post=det(globalKreduced);
    end
    
    [Uglobal,Reactions]=solveq(Kglobal,fglobal,bc);
    
    % --- computation of mechanic elements at the ends of bars --- %
    for i=1:nbars
        ue=Uglobal(Edof(i,2:13));
        ke=elmmat((i-1)*12+1:12*i,:);

        fe=-(ke*ue-fbe(:,i));

        reac_bars(:,i)=fe;
    end
    
    iteration=iteration+1;
    current_plas=0; % to register if there is a plastification in the
                    % current load step

    plastified_bars=zeros(nbars,1);
    for i=1:nbars

        % Detect if any end of this bar (i) has been plastified
        bar_plas_check=0;
        if plast_bars(1,i)~=0 || plast_bars(2,i)~=0
            bar_plas_check=bar_plas_check+1;
        end
        if bar_plas_check==1
            % Detect if the other end has been also plastified
            if abs(reac_bars(5,i))>=Mp(i,1) && ...
               abs(reac_bars(11,i))>=Mp(i,2) % if both ends are plastified
           
                if plast_bars(1,i)==1 && plast_bars(2,i)==0 
                    % The bar is currently Art-Fixed and will be Art-Art
                    current_plas=1;

                    plast_bars(2,i)=1;

                    % Change condition Fixed-Art to Art-Art
                    support(i,3)="Art";
                    
                    % Storing plastic moment / registration
                    mplas=reac_bars(11,i);
                    mbar(i,2)=mplas;
                    
                    plastified_bars(i,1)=2;
                
                elseif plast_bars(1,i)==0 && plast_bars(2,i)==1
                    % The bar is currently Fixed-Art and will be Art-Art
                    current_plas=1;
                    plast_bars(1,i)=1;

                    % Change condition Fixed-Art to Art-Art
                    support(i,2)="Art";
                    
                    % Storing plastic moment / registration
                    mplas=reac_bars(5,i);
                    mbar(i,1)=mplas;
                    
                    plastified_bars(i,1)=1;
                                
                end
            end
        elseif bar_plas_check==0
            if abs(reac_bars(5,i))>=Mp(i,1) && ...
                    abs(reac_bars(11,i))<Mp(i,2)

                current_plas=1;
                mplas=reac_bars(5,i);
                plast_bars(1,i)=1;

                % change condition to Art
                support(i,2)="Art";
                
                % % Storing plastic moment / registration
                mbar(i,1)=mplas;
                    
                plastified_bars(i,1)=1;

            elseif abs(reac_bars(11,i))>=Mp(i,2) && ...
                    abs(reac_bars(5,i))<Mp(i,1)
                current_plas=1;
                plast_bars(2,i)=1;
                mplas=reac_bars(11,i);

                % change condition to Fixed-Art
                support(i,3)="Art";
                
                % Equivalent forces
                mbar(i,2)=mplas;
                
                plastified_bars(i,1)=2;
                
            elseif abs(reac_bars(11,i))>=Mp(i,2) && ...
                    abs(reac_bars(5,i))>=Mp(i,1)
                current_plas=1;
                plast_bars(2,i)=1;
                plast_bars(1,i)=1;
                
                mplas1=reac_bars(5,i);
                mplas2=reac_bars(11,i);

                % change condition to Art-Art
                support(i,2)="Art";
                support(i,3)="Art";
                
                % Storing plastic moment / registration
                mbar(i,1)=mplas1;
                mbar(i,2)=mplas2;
                
                plastified_bars(i,1)=3;

            end
        end
    end
    
    if current_plas==0
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
        continue;
    else
        historyIncLoad=[historyIncLoad,incLoad];
        iter_collection=[iter_collection,iteration];
        
        barPlasNode=[barPlasNode,plastified_bars];
        if sum(seismicforces)<0
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(nbays*i)))];
                else
                    % Relative displacement
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(nbays*i)))-...
                        abs(Uglobal(dofForces(nbays*i-nbays)))];
                    
                end
                
                force_iter=[force_iter;
                            abs(seismicforces(nbays*i))*incLoad]; 
                     
            end
            dispHistFloorLeft=[dispHistFloorLeft,disp_iter];                
            forceFloorLeft=[forceFloorLeft,force_iter];
            
        else
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(nbays*i)))];
                else
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(nbays*i)))-...
                        abs(Uglobal(dofForces(nbays*i-nbays)))];
                        
                end
                force_iter=[force_iter;
                            abs(seismicforces(nbays*i))*incLoad]; 
                
            end
            dispHistFloorRight=[dispHistFloorRight,disp_iter];                
            forceFloorRight=[forceFloorRight,force_iter];
            
        end
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
        
        if sum(seismicforces)<0
            nplas=length(dispHistFloorLeft(1,:))-1;
            du=max(dispHistFloorLeft(1,:))-dispHistFloorLeft(1,2);
            dy=dispHistFloorLeft(1,2);
        else
            nplas=length(dispHistFloorRight(1,:))-1;
            du=max(dispHistFloorRight(1,:))-dispHistFloorRight(1,2);
            dy=forceFloorRight(1,2);
        end
        if du/dy>dydu || nplas>NmaxPlasteps
            break;
        
        else
            continue;
        end
    end

end
            
maxDisplacement=zeros(nfloors,1);
k_direction=zeros(1,nfloors);
pdriftDI=zeros(1,nfloors);
defBasedDI=zeros(1,nfloors);
driftDI=zeros(1,nfloors);
if sum(seismicforces)<0
    
    nd=length(dispHistFloorLeft(1,:));
    for i=1:nfloors
        k_direction(i)=forceFloorLeft(i,2)/dispHistFloorLeft(i,2);
        
        % Plastic drift damage index
        pdriftDI(i)=(max(dispHistFloorLeft(i,:))-dispHistFloorLeft(i,2))/...
                    (Hfloor(i))*100;
    
        % Interstory drift damage index
        driftDI(i)=max(dispHistFloorLeft(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(dispHistFloorLeft(i,:));
        
        delta_u=0.04*Hfloor(i);
        defBasedDI(i)=(maxDisplacement(i)-dispHistFloorLeft(i,2))/...
                      (delta_u-dispHistFloorLeft(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(1)
        if i==1
            plot(dispHistFloorLeft(i,:),...
                forceFloorLeft(i,:),'b -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(dispHistFloorLeft(i,:),...
                forceFloorLeft(i,:),'-',...
                'LineWidth',1.8,'DisplayName',floorText(i,:))
            hold on
        end   
    end
    xlabel('Lateral relative displacement')
    ylabel('Load')
    title('Load-Displacement Historial per floor (-)')
    hold on
    
    Ed=extract(Edof,Uglobal);
    
    figure(2)
    %-----Undeformed mesh-----%
    xlabel('x')
    ylabel('y')
    zlabel('z')
    plotpar=[1,2,1];
    elnum=Edof(:,1);
    eldraw3(Ex,Ey,Ez,plotpar,elnum)
    %----Deformed mesh---------%
    magnfac=100;
    title(strcat('Deformed-undeformed structure - Negative Forces; ',...)
    ' Scale: 1- ', num2str(magnfac)));
    plotpar=[1,3,1];
    [magnfac]=eldisp3(Ex,Ey,Ez,Ed,plotpar,magnfac);
    
else
    
    nd=length(dispHistFloorRight(1,:));
    for i=1:nfloors
        k_direction(i)=forceFloorRight(i,2)/dispHistFloorRight(i,2);
        
        % Plastic drift damage index
        pdriftDI(i)=(max(dispHistFloorRight(i,:))-dispHistFloorRight(i,2))/...
                    (Hfloor(i))*100;
        
        % Interstory drift damage index
        driftDI(i)=max(dispHistFloorRight(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(dispHistFloorRight(i,:));
    
        delta_u=0.04*Hfloor(i);
        
        defBasedDI(i)=(maxDisplacement(i)-dispHistFloorRight(i,2))/...
                      (delta_u-dispHistFloorRight(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(3)
        if i==1
            plot(dispHistFloorRight(i,:),...
                forceFloorRight(i,:),'k -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(dispHistFloorRight(i,:),...
                forceFloorRight(i,:),'-','LineWidth',1.8,...
                'DisplayName',floorText(i,:))
            hold on
        end
    end
    xlabel('Lateral relative displacement')
    ylabel('Load')
    title('Load-Displacement Historial per floor (+)')
    hold on
    
    Ed=extract(Edof,Uglobal);
    
    figure(4)
    %-----Undeformed mesh-----%
    magnfac=100;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(strcat('Deformed-undeformed structure - Positive Forces; ',...)
    ' Scale: 1- ', num2str(magnfac)));
    plotpar=[1,2,1];
    elnum=Edof(:,1);
    eldraw3(Ex,Ey,Ez,plotpar,elnum)
    
    %----Deformed mesh---------%
    plotpar=[1,3,1];
    [magnfac]=eldisp3(Ex,Ey,Ez,Ed,plotpar,magnfac);

end

%---------------------------------end----------------------------------