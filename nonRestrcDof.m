function [ndof,edof]=nonRestrcDof(nnodes,bc,ndofnode)
%------------------------------------------------------------------------
% Syntax:
% [ndof,edof]=nonRestrcDof(nnodes,bc,ndofnode)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the DOF vector containing the DOF that are not restricted
%  based on the pre-stablished restricted DOF vector (bc). Through this
%  vector the reduced stiffness matrix could be computed in order to
%  compare the determinants of each structural system at every stage of
%  analysis.
%  
% 
% INPUT:  nnodes:                total number of nodes of the structure
%
%         bc:                    boundary condition vector, containing the
%                                DOF that are restricted in motion or
%                                displacement. Size: [restric-DOF,2] in
%                                format [DOF,displacement]
%
%         ndofnode:              number of dof per node
%
% OUTPUT: edof:                  is vector containing the non-restricted
%                                DOF. Size: [Non-restric-DOF x 1]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

edof_bc=zeros(nnodes*ndofnode,1);
edof=[];
for i=1:length(bc(:,1))
    
    edof_bc(bc(i,1))=bc(i,1);
end

for i=1:nnodes*ndofnode
    if edof_bc(i)==0
        edof=[edof;i];
    end
end

ndof=length(edof);