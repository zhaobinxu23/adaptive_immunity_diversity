%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)

clc
clear

% 创建一个大图
figure('Position', [100, 100, 1200, 800]);

for i = 1:9
 
    % 创建子图 - 2行5列排列
subplot(3, 3, i);


x0(1) = 0; %% A
x0(2) = 500; %% B
x0(3) = 0; %% C1
x0(4) = 0; %% C2
x0(5) = 1;%% V
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
para(13) = 0.05*i;


[t,y]=ode15s(@pathway_model_plasma_model,[0 50],x0,[],para);







plot(t,y(:,5),'linewidth',2);
hold on


legend('Dynamics of Virus level');
title(['Iteration ', num2str(i)]);
xlabel('Time');
ylabel('Concentration');
grid on;
    
hold off;
end

