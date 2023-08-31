% output
bc = bound.right;
nodesright = unique(bc);
conc_right = mean(u( dofArray(nodesright,2), :));
%%
% output
bc = bound.right;
nodesright = unique(bc);
poten_right = mean(u( dofArray(nodesright,1), :));

%% plot the temp. evolution
% figure
hold on;
t = (0:n_time)*dt;
plot(t, mean(u(dofArray(1:nnode,1),:)),'o','MarkerSize',8);
plot(t, mean(u(dofArray(1:nnode,2),:)),'o','MarkerSize',8);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 4 0 1]);
box on
%%
% figure
hold on;
t = (0:n_time)*dt;
plot(t,[0; interfluxstore(:,1)]/9e7,'r-s','Linewidth',2);
% ---------------------------------------------------------------------------
%% plot the temp. at last step --- rve
figure
hold on;
plot(coords(:,1),u(1:nnode,end),'o','MarkerSize',5);
plot(coords(:,1),u(nnode+1:end,end),'s','MarkerSize',5);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 4 0 1]);
box on

%%
x = coords ((coords(:,2) == 0.0),1); 
y = u_macro(coords(:,2) == 0.0);

data = [x y];
[~,idx] = unique(data(:,1));
data = data(idx,:);

figure;
plot(x,y,'o')

% ---------------------------------------------------------------------------
%%
data = [t; mean(u(dofArray(1:nnode,1),:)); mean(u(dofArray(1:nnode,2),:))]';
fileID = fopen('eg1_macroTemp_fe2.txt','w');
fprintf(fileID,'%12.6e, %12.6e, %12.6e \n',data');
fclose(fileID);
%%
t = (0:n_time)*dt;
data = [t', [0; interfluxstore(:,1)]];
fileID = fopen('eg1_interflux_fe2.txt','w');
fprintf(fileID,'%12.6e, %12.6e \n',data');
fclose(fileID);

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