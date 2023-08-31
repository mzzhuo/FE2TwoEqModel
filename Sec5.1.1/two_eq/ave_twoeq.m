% analyze fully resolved solution

% parameters
% L = 1/140; %0.00357143; 
% H = 1/140; %0.00357143; 
% nx = 1;  % circle number in x direction
% ny = 1;   % circle number in y direction
% r_ca = L/nx * 0.3250;


T_mat = zeros(timesteps,2); % first colum time, sec column temp.
% 
node_mat = unique(connect.mat);

for j = 1:timesteps
    count = 0;
    for i = 1:length(node_mat)
%         xcoord = coords(node_mat(i),1);
%         if (xcoord >= 0) && (xcoord <= L/nx)
            count = count + 1;
            T_mat(j,2) = T_mat(j,2) + u(dofArray(node_mat(i),1),j); % 1 means temp 
%         end
    end
    T_mat(j,2) = T_mat(j,2)./count; 
end
T_mat(:,1) = (0:n_time)*dt;


T_inc = zeros(timesteps,2); % first colum time, sec column temp.
% 
node_inc = unique(connect.inc);

for j = 1:timesteps
    count = 0;
    for i = 1:length(node_inc)
%         xcoord = coords(node_inc(i),1);
%         if (xcoord >= 0) && (xcoord <= L/nx)
            count = count + 1;
            T_inc(j,2) = T_inc(j,2) + u(dofArray(node_inc(i),1),j); % 1 means temp 
%         end
    end
    T_inc(j,2) = T_inc(j,2)./count; 
end
T_inc(:,1) = (0:n_time)*dt;

%%
T_rve = zeros(timesteps,2); % first colum time, sec column temp.
% 
node = unique([connect.mat; connect.inc]);

for j = 1:timesteps
    count = 0;
    for i = 1:length(node)
%         xcoord = coords(node_mat(i),1);
%         if (xcoord >= 0) && (xcoord <= L/nx)
            count = count + 1;
            T_rve(j,2) = T_rve(j,2) + u(dofArray(node(i),1),j); % 1 means temp 
%         end
    end
    T_rve(j,2) = T_rve(j,2)./count; 
end
T_rve(:,1) = (0:n_time)*dt;
%%
data = [T_rve(:,1), T_rve(:,2)];
fileID = fopen('eg1_aveTemp_singlescale.txt','w');
fprintf(fileID,'%12.6e, %12.6e \n',data');
fclose(fileID);
%%
figure;
hold on;
plot(T_mat(:,1),T_mat(:,2),'-','MarkerSize',5)
plot(T_inc(:,1),T_inc(:,2),'-','MarkerSize',5)

%%
data = [T_mat(:,1), T_mat(:,2), T_inc(:,2)];
fileID = fopen('eg1_aveTemp_singlescale.txt','w');
fprintf(fileID,'%12.6e, %12.6e, %12.6e \n',data');
fclose(fileID);
%%
data = [T_mat(:,1), [0;sour_mat]];
fileID = fopen('eg1_interflux_singlescale.txt','w');
fprintf(fileID,'%12.6e, %12.6e \n',data');
fclose(fileID);

%%
% homogenized conc/potential in matrix
%
%
% T_mat = zeros(nx,2); % first colum count of node, sec column potential; third column conc 
% % 
% node_mat = unique(connect.mat);
% 
% for i = 1:length(node_mat)
%     xcoord = coords(node_mat(i),1);
%     for j = 1:nx
%         if (xcoord >= (j-1)*L/nx) && (xcoord <= j*L/nx)
%             T_mat(j,1) = T_mat(j,1) + 1; 
%             T_mat(j,2) = T_mat(j,2) + u(dofArray(node_mat(i),1),end); % 1 means temp 
%         end
%     end
% end
% T_mat(:,2) = T_mat(:,2)./T_mat(:,1); 
% T_mat(:,1) = L/nx/2 + (0:nx-1)*L/nx;  

%%
% figure;
% hold on;
% plot([0;T_mat(:,1)],[1;T_mat(:,2)],'-s','Linewidth',3,'MarkerSize',3)

%%
% homogenized conc/potential in inclusion
%
%
% T_inc = zeros(nx,2); % first colum count of node, sec column potential; third column conc 
% % 
% node_inc = unique(connect.inc);
% 
% for i = 1:length(node_inc)
%     xcoord = coords(node_inc(i),1);
%     for j = 1:nx
%         if (xcoord >= (j-1)*L/nx) && (xcoord <= j*L/nx)
%             T_inc(j,1) = T_inc(j,1) + 1; 
%             T_inc(j,2) = T_inc(j,2) + u(dofArray(node_inc(i),1),end); % 1 means temp 
%         end
%     end
% end
% T_inc(:,2) = T_inc(:,2)./T_inc(:,1); 
% T_inc(:,1) = L/nx/2 + (0:nx-1)*L/nx;  
%%
% figure;
% hold on;
% plot(T_inc(:,1),T_inc(:,2),'-o','Linewidth',3,'MarkerSize',3)

%%
% fileID = fopen('fullscale16_15unit.tikz','w');
% fprintf(fileID,'%12.6f, %12.6f, %12.6f \n',T_mat');
% fclose(fileID);
% 
%% paraview data
% data = [fullscale4xaxis(:,4), fullscale4xaxis(:,1), fullscale4xaxis(:,2)];
% fileID = fopen('fullscale4_xaxis.tikz','w');
% fprintf(fileID,'%12.6f, %12.6f, %12.6f \n',data');
% fclose(fileID);



%%
% F_int = - K(nnode+1:end,1:nnode)' * u(nnode+1:end,end);
% 
% node_int = interface.nodes;
% 
% % source_mat = sum(F_int(node_int([1:40 401:440],2)))/(1/20 * 0.3478 *2*pi);
% % 
% % source_inc = sum(F_int(node_int([1:40 401:440],3)))/(1/20 * 0.3478 *2*pi);
% 
% source_mat = sum(F_int(node_int([1:40 401:440],2)))/(1/10)^2;
% 
% source_inc = sum(F_int(node_int([1:40 401:440],3)))/(1/10)^2;

%
% source_mat = sum(F_int(node_int(1:20,2)))/(1/20)^2;
% source_mat2 = sum(F_int(node_int(401:420,2)))/(1/20)^2;
% 
% source_mat3 = sum(F_int(node_int(21:40,2)))/(1/20)^2;
% source_mat4 = sum(F_int(node_int(421:440,2)))/(1/20)^2;
% source_inc = sum(F_int(node_int(1:20,3)))/(1/20)^2;

%%

% F_int = - K(nnode+1:end,1:nnode)' * u(nnode+1:end,end);
% 
% node_int = interface.nodes;
% 
% % source_mat = sum(F_int(node_int(1:20,2)))/(1/10 * 0.3478 *2*pi);
% % 
% % source_inc = sum(F_int(node_int(1:20,3)))/(1/10 * 0.3478 *2*pi);
% 
% source_mat = sum(F_int(node_int(1:20,2)))/(1/10)^2;
% 
% source_inc = sum(F_int(node_int(1:20,3)))/(1/10)^2;
