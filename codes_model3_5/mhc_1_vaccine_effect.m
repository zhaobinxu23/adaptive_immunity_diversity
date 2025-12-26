%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)
clc
clear



x0(1) = 0; %% Tu  
x0(2) = 0; %% V  
x0(3) = 0; %% C  
x0(4) = 1e5; %% A  
x0(5) = 1e3;%% Tc  



para(1) = 1.1; % k1
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
para(18) = 0.05; % decay ratio of monoclonal antiobdy






[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);
plot(t,y(:,1),'linewidth',5);
hold on

para(14) = 4.0e-4;% k5

[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);
plot(t,y(:,1),'linewidth',5);
hold on

para(14) = 6.0e-4;% k5

[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);
plot(t,y(:,1),'linewidth',5);
hold on

xlabel('Time');
ylabel('Concentration');







