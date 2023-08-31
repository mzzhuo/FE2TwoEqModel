%% plot the temp. evolution
figure
hold on;
t = (0:n_time)*dt;
plot(t, mean(u(dofArray(1:nnode,1),:)),'o','MarkerSize',8);
% plot(t, mean(u(dofArray(1:nnode,2),:)),'s','MarkerSize',8);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 4 0 1]);
box on

%% plot the concentration at last step
% figure
hold on;
plot(coords(:,1),u(:,end),'o','MarkerSize',3);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 4 0 1]);
box on
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

figure;
plot(x,y,'o')

% ---------------------------------------------------------------------------
%%
t = (0:n_time)*dt;
data = [t; mean(u(dofArray(1:nnode,1),:))]';
fileID = fopen('macroTemp_oneq.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);
