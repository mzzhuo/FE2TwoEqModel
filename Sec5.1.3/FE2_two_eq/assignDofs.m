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
ndof = 2;       % each node has 2 dofs: macro temp. for matrix and inclusion
% 
%% Assign dofArray: An array of dof indices: phi c lambda
% 
% nDofLambda = size(addLambda,1); % each pair one dof
%%
dofArray = zeros( nnode,ndof ); 
% dofs for conc/temp
dofArray(1:nnode,1) = 1:nnode; 
dofArray(1:nnode,2) = nnode+1:2*nnode; 
%%
% 
ndofs = ndof*nnode;
% 
% conc dofs   
% dofsC = dofArray( 1:nnode,1 );
