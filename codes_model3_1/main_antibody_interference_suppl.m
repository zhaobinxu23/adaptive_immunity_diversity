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
x0(5) = 1e8;%% V
x0(6) = 0; %% ASC cell


para(1) = 1e-7; % k1
para(2) = 1e-14; % k-1
para(3) = 2;% k2
para(4) = 1e-3;% alpha
para(5) = 0.005;% d1
para(6) = 0.02;% d2
para(7) = 0.5;% d3
para(8) = 0.5;% d4
para(9) = -0.02;% k3
para(10) = 1e-7;
para(11) = 0.01; %% decay ratio of ASC cell
para(12) = 1e5;


[t,y]=ode15s(@pathway_model_plasma_model,[0 100],x0,[],para);


x0(1) = interp1(t,y(:,1),100);
x0(2) = interp1(t,y(:,2),100);
x0(3) = interp1(t,y(:,3),100);
x0(4) = interp1(t,y(:,4),100);
x0(5) = 10^(0.5*(i-1))*1e6;
x0(6) = interp1(t,y(:,6),100); %% ASC cell

[t_new,y_new]=ode15s(@pathway_model_plasma_model,[100 1000],x0,[],para);




plot(t,y(:,1),'linewidth',2);
hold on


plot(t_new,y_new(:,1),'linewidth',2);
hold on

legend('Dynamics of antibody levels after primary vaccination',  'Dynamics of antibody levels after secondary vaccination');
title(['Iteration ', num2str(i)]);
xlabel('Time');
ylabel('Concentration');
grid on;
    
hold off;
end

