%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)
clc
clear




x0(1) = 0; %% Tu  466.70152158166456501930952072144
x0(2) = 0; %% V  9999.9948277755174785852432250977
x0(3) = 0; %% C  1933403.043163299560546875  1009094.06948363780975341796875
x0(4) = 1e5; %% A  1933404043.163330078125  5000000000.0
x0(5) = 1e3;%% Tc  264321.40618126433498247251918656  1009137.72138082981109619140625


para(1) = 1.0; % k1
para(2) = 0.05; % d1
para(3) = 10;% β1
para(4) = 1e7;% km
para(5) = 1e-4;% α1
para(6) = 1e-7;% k2
para(7) = 0;% k-2
para(8) = 1e4;% ρ1
para(9) = 1;% k3
para(10) = 2;% k4
para(11) = 1e-3;% d2
para(12) = 1e2;% Π2
para(13) = 0.1;% ρ2
para(14) = 2.0e-4;% k5
para(15) = 1;% Π3
para(16) = 1e-3;% d3
para(17) = 100;






[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);

title('Dynamics of tumor cells');
xlabel('Time');
ylabel('Concentration');
%% Modeling the increase of cancer cell replication on overall tumor dynamics
para(1) = 1.1; % k1
[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);

title('Dynamics of tumor cells');
xlabel('Time');
ylabel('Concentration');

%% Modeling the decrease of antibody regernation constant on overall tumor dynamics
para(1) = 1.0; % k1
para(10) = 1.9;% k4 from 2 to 1.9
[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);

title('Dynamics of tumor cells');
xlabel('Time');
ylabel('Concentration');

%% Modeling the further decrease of antibody regernation constant on overall tumor dynamics
para(1) = 1.0; % k1
para(10) = 1.8;% k4 from 2 to 1.8
[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);

title('Dynamics of tumor cells');
xlabel('Time');
ylabel('Concentration');


% x0(1) =  interp1(t,y(:,1),20);
% x0(2) =  interp1(t,y(:,2),20);
% x0(3) = interp1(t,y(:,3),20);
% x0(4) = interp1(t,y(:,4),20);
% x0(5) = interp1(t,y(:,5),20);
% 
% para(10) = 3;% k4
% 
% 
% [z,y_2]=ode15s(@pathway_model_tumor_5,[20 100],x0,[],para);
% 
% plot(z,y_2(:,1),'linewidth',2);
% hold on
% para(10) = 1.5;% k4
% para(14) = 1.5e-4;% k5
% 
% [t_new,y_new]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);
% % 
% plot(t_new,y_new(:,1),'linewidth',2);
% hold on
% 
% para(10) = 1.78; 
% para(14) = 1.5e-4;% k5
% para(1) = 1.0; % k1
% 
% 
% [z,y_2]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);
% % 
% plot(z,y_2(:,1),'linewidth',2);
% hold on
% plot(z,y_2(:,4),'linewidth',2);
% hold on
% plot(z,y_2(:,5),'linewidth',2);
% hold on
% x0(1) =  0.1*interp1(z,y_2(:,1),20);
% x0(2) =  interp1(z,y_2(:,2),20);
% x0(3) = interp1(z,y_2(:,3),20);
% x0(4) = interp1(z,y_2(:,4),20);
% x0(5) = interp1(z,y_2(:,5),20);
% 
% para(10) = 2;% k4
% para(14) = 2e-4;% k5
% 
% [z,y_2]=ode15s(@pathway_model_tumor_5,[20 100],x0,[],para);
% 
% % plot(z,y_2(:,1),'linewidth',2);
% % hold on
% plot(z,y_2(:,4),'linewidth',2);
% hold on
% plot(z,y_2(:,5),'linewidth',2);
% hold on
% x0(6) = interp1(t,y(:,6),10000);
% 
% x0(7) = interp1(t,y(:,7),10000);
% 
% 
% para(10) = 1;% k4
% para(14) = 1e-4;% k5
% 
% [t_new,y_new]=ode15s(@pathway_model_tumor_7,[0 1000],x0,[],para);
% 
% plot(t_new,y_new(:,1),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t,y(:,5),'linewidth',2);
% hold on
% %%
% x0(1) = interp1(t,y(:,1),10);
% x0(2) = interp1(t,y(:,2),10);
% x0(3) = interp1(t,y(:,3),10);
% x0(4) = interp1(t,y(:,4),10);
% x0(5) = interp1(t,y(:,5),10);
% x0(6) = interp1(t,y(:,6),10);
% x0(7) = interp1(t,y(:,7),10);
% 
% para(10) = 2.0;% k4
% para(14) = 1e-3;% k5
% para(13) = 0.1;% ρ2
% 
% [t_new,y_new]=ode15s(@pathway_model_tumor_7,[10 100],x0,[],para);
% 
% plot(t_new,y_new(:,1),'linewidth',2);
% hold on
% plot(t_new,y_new(:,4),'linewidth',2);
% hold on




