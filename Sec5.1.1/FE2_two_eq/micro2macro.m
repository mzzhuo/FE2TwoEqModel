%% flux with interface nodes
% volume = vol_mat + vol_inc;
% node_int = interface.nodes;
% 
% bounodes_mat = [bounodes; node_int(:,2)];
% F_bi = F_bond + F_infc;
% 
% cols_mat = dofArray(bounodes_mat,1);
% D_b = coords(bounodes_mat,:);
% flux_mat = - D_b' * F_bi(cols_mat)/volume;
% 
% % consistent matrix
% delta = eye(ndofs,ndofs);
% delta = delta([nnode+1:nnode+nlamb+3 nnode+nlamb+5+1:end],:);
% rhs = [zeros(ndof*nnode,4); -C];
% solu = K \ rhs;
% % stiff = K1'*delta*solu;
% stiff = [K_peri; K_corn; K_inte]'*delta*solu;
% % 
% dFlux_mat = D_b'*stiff(cols_mat,:)/volume; 

%%
% bounodes_inc = node_int(:,3);
% cols_inc = dofArray(bounodes_inc,1);
% D_b = coords(bounodes_inc,:);
% flux_inc = - D_b' * F_bi(cols_inc)/volume;
% 
% delta = eye(ndofs,ndofs);
% delta = delta([nnode+1:nnode+nlamb+3 nnode+nlamb+5+1:end],:);
% rhs = [zeros(ndof*nnode,4); -C];
% solu = K \ rhs;
% stiff = [K_peri; K_corn; K_inte]'*delta*solu;
% % 
% dFlux_inc = D_b'*stiff(cols_inc,:)/volume; 
%% flux no interfacial part
% calculate homegenized 
volume = vol_mat + vol_inc;

% bounodes = [ bounodes; node_int(:,2) ];
cols = dofArray(bounodes,1);
D_b = coords(bounodes,:);
flux_mat = - D_b' * F_bond(cols)/volume;
% flux = - coords'*F_int/volume;

rhs = [zeros(ndof*nnode,4); -C];
solu = K \ rhs;
% stiff = K1'*delta*solu;
solu = solu(nnode+1:nnode+nlamb+3,:);
stiff = [K_peri; K_corn]'*solu;
% 
D_b  = coords(bounodes,:);
cols = dofArray(bounodes,1);
dFlux_mat = D_b'*stiff(cols,:)/volume; 

%% consistent matrix
node_int = interface.nodes;
rhs = [zeros(ndof*nnode,4); -C];
solu = K \ rhs;
% stiff = K1'*delta*solu;
solu = solu(nnode+nlamb+5+1:end,:);
stiff = K_inte'*solu;
% 
% node_int = interface.nodes;
% 
% interfacial flux outward positive, inward flux is positive source
sour_mat = -sum(F_infc(node_int(:,2),:))/volume; 
dSou_mat = -sum(stiff(node_int(:,2),:))/volume; 
% 
sour_inc = -sum(F_infc(node_int(:,3),:))/volume; 
dSou_inc = -sum(stiff(node_int(:,3),:))/volume; 

%%
% %% check internal nodes 
% allnodes = unique(connect);
% allnodes(bounodes) = [];
% internodes = allnodes;
% internodes = [internodes; internodes+nnode];
% Stiff(internodes,:);
% check = F_int(internodes,:);
% end