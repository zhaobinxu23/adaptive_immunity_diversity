%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)
clc
clear

figure('Position', [100, 100, 1200, 800]);

for i = 1:9
    % 创建子图 - 2行5列排列
    subplot(3, 3, i);
    



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






[t,y]=ode15s(@pathway_model_tumor_5,[0 5],x0,[],para);
plot(t,y(:,1),'linewidth',5);
hold on

x0(1) =  interp1(t,y(:,1),5);
x0(2) =  interp1(t,y(:,2),5)+(i-1)*10^8;
x0(3) = interp1(t,y(:,3),5);
x0(4) = interp1(t,y(:,4),5);
x0(5) = interp1(t,y(:,5),5);


[t,y]=ode15s(@pathway_model_tumor_5,[5 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);
hold on

title(['Dynamics of tumor cells with neoantigen dosage =  ', num2str(i-1),'*10^8 ']);
xlabel('Time');
ylabel('Concentration');
hold off
end






