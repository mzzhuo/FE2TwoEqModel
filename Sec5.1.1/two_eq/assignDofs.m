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
ndof = 1;       % each node has 2 dofs: potential and concentration
% 
%% Assign dofArray: An array of dof indices: phi c lambda
% matched interface node number
nnodepair = size(interface.nodes,1); 
%%
dofArray = zeros( nnode+nnodepair,ndof ); 
% dofs for conc/temp
dofArray(1:nnode+nnodepair,1) = 1:nnode+nnodepair; 
%%
% 
ndofs = ndof*nnode + nnodepair;
% 
% conc dofs   
% dofsC = dofArray( 1:nnode,1 );
