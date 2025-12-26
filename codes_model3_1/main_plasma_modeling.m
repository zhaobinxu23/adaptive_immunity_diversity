%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)


clc
clear

x0(1) = 0; %% A
x0(2) = 500; %% B
x0(3) = 0; %% C1
x0(4) = 0; %% C2
x0(5) = 10;%% V
x0(6) = 0; %% ASC cell


para(1) = 1e-7; % k1
para(2) = 1e-14; % k-1
para(3) = 2;% k2
para(4) = 1e-3;% alpha
para(5) = 0.005;% d1
para(6) = 0.02;% d2
para(7) = 0.5;% d3
para(8) = 0.5;% d4
para(9) = 0.5;% k3
para(10) = 1e-7;
para(11) = 0.01; %% decay ratio of ASC cell
para(12) = 1e5;
para(13) = 0.5;


[t,y]=ode15s(@pathway_model_plasma_model,[0 150],x0,[],para);

% plot(t,y(:,1),'linewidth',2);
% hold on
% plot(t,y(:,2),'linewidth',2);
% hold on
% plot(t,y(:,3),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t,y(:,5),'linewidth',2);
% hold on

% x0(1) = interp1(t,y(:,1),500);
% x0(2) = interp1(t,y(:,2),500);
% x0(3) = interp1(t,y(:,3),500);
% x0(4) = interp1(t,y(:,4),500);
% x0(5) = 1e10;
% x0(6) = interp1(t,y(:,6),500); %% ASC cell
% 
% [t_new,y_new]=ode15s(@pathway_model_plasma_model,[500 1000],x0,[],para);
% 
plot(t,log10(y(:,5)),'linewidth',2);
hold on
% plot(t,y(:,2),'linewidth',2);
% hold on
% plot(t,log10(y(:,5)),'linewidth',2);
% hold on

para(9) = 0.55;% k3

[t_new,y_new]=ode15s(@pathway_model_plasma_model,[0 150],x0,[],para);

plot(t_new,log10(y_new(:,5)),'linewidth',2);
hold on
% plot(t,y(:,2),'linewidth',2);
% hold on
% plot(t,log10(y(:,5)),'linewidth',2);
% hold on


% plot(t,y(:,1),'linewidth',2);
% hold on
% plot(t,y(:,3),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t,y(:,5),'linewidth',2);
% hold on
% 

% plot(t_new,y_new(:,1),'linewidth',2);
% % hold on
% plot(t_new,y_new(:,1),'linewidth',2);
% hold on
% plot(t,y(:,3),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t_new,y_new(:,5),'linewidth',2);
% hold on
% 
