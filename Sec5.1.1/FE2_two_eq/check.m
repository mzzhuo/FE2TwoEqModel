
node_inc = unique(connect.inc);

node_mat = unique(connect.mat);


%%
cols = dofArray(bounodes,1); 
flux_bound = sum(F_bond(cols))/volume;
%%
% cons_mat = sum(F_cons(node_mat))/volume;
% cons_inc = sum(F_cons(node_inc))/volume;

cons_mat = sum(F_cons_mat)/volume;
cons_inc = sum(F_cons_inc)/volume;


%%
u_mat * sum( K_consA' / coefs.hc0 / vol_mat  * u(nnode+nlamb+4) )
u_inc * sum( K_consB' / coefs.hc1 / vol_inc  * u(nnode+nlamb+5) ) 


%%
u_mat *  cons_mat
u_inc *  cons_inc

%%
sum(u(1:nnode) .* F_cons_mat) /volume
sum(u(1:nnode) .* F_cons_inc) /volume
%%

sum(u(node_int(:,2)) .* F_infc(node_int(:,2),:)) / volume

sum(u(node_int(:,3)) .* F_infc(node_int(:,3),:)) / volume


%%
K_consA * u(1:nnode) / ( coefs.hc0 * vol_mat )


%% 

x_mat = coords(node_int(:,2),:);
x_mat' * F_infc(node_int(:,2),:) / volume

x_inc = coords(node_int(:,2),:);
x_inc' * F_infc(node_int(:,3),:) / volume


%%

coords' * F_cons_mat / volume


coords' * F_cons_inc / volume

%%
cencoord = [(max(coords(:,1))+min(coords(:,1)))/2 (max(coords(:,2))+min(coords(:,2)))/2];

cencoord * sum(F_cons_mat) / volume
