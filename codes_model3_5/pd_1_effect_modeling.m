%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)
clc
clear

figure('Position', [100, 100, 1200, 800]);

for i = 1:6
    % 创建子图 - 2行5列排列
    subplot(3, 2, i);
    



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
para(13) = 10^(i-6);% ρ2
para(14) = 2.0e-4;% k5
para(15) = 1;% Π3
para(16) = 1e-3;% d3
para(17) = 100;






[t,y]=ode15s(@pathway_model_tumor_5,[0 100],x0,[],para);

plot(t,y(:,1),'linewidth',5);
hold on

title(['Dynamics of tumor cells when ρ2 = 10^ ', num2str(i-6)]);
xlabel('Time');
ylabel('Concentration');
hold off
end






