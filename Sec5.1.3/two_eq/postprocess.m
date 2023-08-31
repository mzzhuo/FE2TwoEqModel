% output
figure;
hold on;
plot(T_mat(:,1),T_mat(:,2),'o','MarkerSize',5)
plot(T_inc(:,1),T_inc(:,2),'s','MarkerSize',5)
%% plot the concentration at last step
figure
hold on;
plot(coords(:,1),u(1:nnode,end),'o','MarkerSize',3);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 0.004 6 24]);
box on
%%
incnodes = unique(connect.inc);
plot(coords(incnodes,1),u(incnodes,end),'ro','MarkerSize',3);

%%
incnodes = unique(connect.inc);
plot(coords(incnodes,1),u(incnodes,end),'ro','MarkerSize',3);
% ---------------------------------------------------------------------------
%%
hold on;
plot(temp0(:,1),temp0(:,2),'r-','Linewidth',3);

%%
x = coords ((coords(:,2) == 0.0),1); 
y = u_macro(coords(:,2) == 0.0);

data = [x y];
[~,idx] = unique(data(:,1));
data = data(idx,:);

% figure;
% plot(x,y,'o')

% ---------------------------------------------------------------------------
%%
data = T_mat;
fileID = fopen('macroTemp_mat.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);
%%
% ( coefs.hc0 * (T_mat(2)*vol_mat + T_inc(2)*vol_inc) ) / vol_inc / 2;
% heat goes to the matrix
coefs.hc0 * T_mat(end,2) * (volume*0.6682) / (dt * n_time)
%%
% heat flux by interface
sum(F_int(node_int(:,2),:))/volume
%%
leftBC = bound.left;
nodesleft = unique(leftBC);
rightBC = bound.right;
nodesright = unique(rightBC);
lowerBC = bound.lower;
nodeslower = unique(lowerBC);
upperBC = bound.right;
nodesupper = unique(upperBC);
nodesbound = [nodesleft;nodesright;nodeslower;nodesupper];
sum(F_int(nodesbound))
    
    