%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)

clc
clear

% 创建一个大图
figure('Position', [100, 100, 1200, 800]);


for i = 1:12
    % 创建子图 - 2行5列排列
subplot(4, 3, i);

x0(1) = 0; % drug
x0(2) = 1; % virus
x0(3) = 0; % drug-virus complex
x0(4) = 1000; % antibody
x0(5) = 0; % virus-antibody complex
x0(6) = 0; % drug-virus-antibody complex

para(1) = 1e-7; % k1 drug binding constant
para(2) = 1e-14; % k-1  drug dissociation constant
para(3) = 1e-2;% k2 % drug decay constant
para(4) = 0.5;% k3  complex decay constant
para(5) = 1;% k4 virus replication constant
para(6) = 1e-7;% k5 antibody binding constant
para(7) = 1e-14;% k-5  antibody dissociation constant
para(8) = 10;% k6 antibody replenish constant
para(9) = 1e-2;% k7 antibody decay constant
para(10) = 2;% k8  feedback constant
para(11) = 0.1;



[t,y]=ode15s(@pathway_model_with_drug_inhib,[0 20],x0,[],para);

% plot(t,y(:,1),'linewidth',2);
% hold on
plot(t,y(:,2),'linewidth',2);
hold on


%%
x0(1) = i*1e8;
x0(2) = interp1(t,y(:,2),20);
x0(3) = 0;
x0(4) = interp1(t,y(:,4),20);
x0(5) = interp1(t,y(:,5),20);
x0(6) = 0;

[t_new,yy]=ode15s(@pathway_model_with_drug_inhib,[20 100],x0,[],para);


plot(t_new,yy(:,2),'linewidth',2);
hold on

legend('Virus dynamics before therapy',  'Virus dynamics after therapy');
title(['Iteration ', num2str(i)]);
xlabel('Time');
ylabel('Concentration');
grid on;
    
hold off;
end


