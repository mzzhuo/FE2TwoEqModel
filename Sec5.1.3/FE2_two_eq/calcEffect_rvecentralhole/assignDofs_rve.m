%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 
% mesh parameters
ncoord = 2;     % 2D problem
if strcmp(elemType,'Tri3')
    nelnodes = 3;
elseif strcmp(elemType,'Quad4') 
    nelnodes = 4;   % triangle element
end

nfacenodes = 2; % line element
ndof = 1;       % each node has 1 dofs: concentration/temp.
% 
%% Assign dofArray: An array of dof indices: phi c lambda
% 
% nDofLambda = size(addLambda,1); % each pair one dof
%%
dofArray = zeros( nnode,ndof); 

% dofs for conc/temp
dofArray(1:nnode,1) = 1:nnode; 

% no. of matched node pair/lagrange multiplier for 
% periodic bcs, not including 4 corners nodes
nlamb = length(bound.master);
% 
% all boundry nodes including corner nodes
bounodes = [bound.master; bound.slave; bound.corner];
bounodes = unique(bounodes);
% 
if ( bc == 'linear' )
    ndofs = ndof*nnode + ndof*nfix + ndof;
elseif ( bc == 'periodic' )
    ndofs = ndof*nnode + ndof*nlamb + 4*ndof; % 3 corners nodes (1 for ref.) + 2 conservation
else
    error('wrong boundary condition type!');
end
% 