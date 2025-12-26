% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
clc
clear

x0 = zeros(1402,1);
%%  IgM distribution
mu = 9; % 均值
sigma = 0.8; % 标准差


% prob_A(1) = 0.001;
% prob_A(2:9) = 0;
prob_A(1) = normcdf(5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end

for i = 1:10
    for j = 1:10
        AA(10*(i-1)+j) = prob_A(i)*prob_A(j);
    end
end

M_1 = 1e15;
M_2 = 4e18;
G_1 = 1e15;
G_2 = 4e19;
E = 1e18;
C_1 = 1e13;
C_2 = 4e16;
C_3 = 1e13;
C_4 = 4e17;
P_M = 5e7;
P_G = 1e8;


for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = M_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+100) = M_2*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+200) = G_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+300) = G_2*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+400) = C_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+500) = C_2*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+600) = C_3*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+700) = C_4*prob_A(i)*prob_A(j);
    end
end

for i = 801:1200
   
 x0(i) = 0;
  
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1200) = P_M*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1300) = P_G*prob_A(i)*prob_A(j);
    end
end


 x0(1401) = E;
 x0(1402) = 10;





para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 2.2e17+1e13; % replenish constant pi 1 
para(4) = 0.5e13;% replenish constant pi 2
para(5) = 0.01; % decay constant of BCR IgM
para(6) = 0.005;% decay constant of BCR IgG
para(7) = 1; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 5e-7;% k2' feedback constant on PC cell regeneration
para(9) = 0.1; % decay constant of plasma Cell IgM
para(10) = 0.05;%% decay constant of plasma Cell IgG
para(11) = 4.4e9;% production constant of IgM
para(12) = 1.2e10;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.02;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 1.2; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;





[t, y] = ode15s(@pathway_model_many_antibody_immune_include_plasma,[0 10000],x0,[],para,para_new,para_new_1,AA);


% [t, y]=ode15s(@pathway_model_many_antibody_immune_include_T,[0 10],x0,[],para,para_new,prob_A,options);


% plot(t,y(:,1402),'linewidth',2); %% generate figure 1


IgG_primary_3D = plot_figure_IgG(t,y); %  generate figure 3


%%
strong_antibody = interp1(t,y(:,251),10000)+interp1(t,y(:,261),10000)+interp1(t,y(:,262),10000)+interp1(t,y(:,271),10000)+interp1(t,y(:,272),10000)+interp1(t,y(:,273),10000)+interp1(t,y(:,281),10000)+interp1(t,y(:,282),10000)+interp1(t,y(:,283),10000)+interp1(t,y(:,284),10000)...
    +interp1(t,y(:,291),10000)+interp1(t,y(:,292),10000)+interp1(t,y(:,293),10000)+interp1(t,y(:,294),10000)+interp1(t,y(:,295),10000);

very_strong_antibody = interp1(t,y(:,261),10000)+interp1(t,y(:,271),10000)+interp1(t,y(:,272),10000)+interp1(t,y(:,281),10000)+interp1(t,y(:,282),10000)+interp1(t,y(:,283),10000)...
    +interp1(t,y(:,291),10000)+interp1(t,y(:,292),10000)+interp1(t,y(:,293),10000)+interp1(t,y(:,294),10000);

super_strong_antibody = interp1(t,y(:,271),10000)+interp1(t,y(:,281),10000)+interp1(t,y(:,282),10000)...
    +interp1(t,y(:,291),10000)+interp1(t,y(:,292),10000)+interp1(t,y(:,293),10000);

super_super_strong_antibody = interp1(t,y(:,281),10000)...
    +interp1(t,y(:,291),10000)+interp1(t,y(:,292),10000);

super_super_super_strong_antibody = interp1(t,y(:,291),10000);

for iii = 5
    iii
x0 = zeros(1402,1);
%%  IgM distribution
mu = 9; % 均值
sigma = 0.8; % 标准差


% prob_A(1) = 0.001;
% prob_A(2:9) = 0;
prob_A(1) = normcdf(5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end

for i = 1:10
    for j = 1:10
        AA(10*(i-1)+j) = prob_A(i)*prob_A(j);
    end
end

M_1 = 1e15;
M_2 = 4e18;
G_1 = 1e15;
G_2 = 4e19;
E = 1e18;
C_1 = 1e13;
C_2 = 4e16;
C_3 = 1e13;
C_4 = 4e17;
P_M = 5e7;
P_G = 1e8;



for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = M_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+100) = M_2*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+200) = G_1*prob_A(i)*prob_A(j)*(1-10^(-iii))+10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000);
        zz0(10*(i-1)+j+200) = 10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000)+strong_antibody*prob_A(i)*prob_A(j);
        zz1(10*(i-1)+j+200) = 10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000)+very_strong_antibody*prob_A(i)*prob_A(j);
        zz2(10*(i-1)+j+200) = 10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000)+super_strong_antibody*prob_A(i)*prob_A(j);
        zz3(10*(i-1)+j+200) = 10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000)+super_super_strong_antibody*prob_A(i)*prob_A(j);
        zz4(10*(i-1)+j+200) = 10^(-iii)*interp1(t,y(:,10*(i-1)+j+200),10000)+super_super_super_strong_antibody*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+300) = G_2*prob_A(i)*prob_A(j)*(1-10^(-iii))+10^(-iii)*interp1(t,y(:,10*(i-1)+j+300),10000);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+400) = C_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+500) = C_2*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+600) = C_3*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+700) = C_4*prob_A(i)*prob_A(j);
    end
end

for i = 801:1200
   
 x0(i) = 0;
  
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1200) = P_M*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1300) = P_G*prob_A(i)*prob_A(j)*(1-10^(-iii))+10^(-iii)*interp1(t,y(:,10*(i-1)+j+1300),10000);
    end
end


 x0(1401) = E;
 x0(1402) = 10;





para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 2.2e17+1e13; % replenish constant pi 1 
para(4) = 0.5e13;% replenish constant pi 2
para(5) = 0.01; % decay constant of BCR IgM
para(6) = 0.005;% decay constant of BCR IgG
para(7) = 1; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 5e-7;% k2' feedback constant on PC cell regeneration
para(9) = 0.1; % decay constant of plasma Cell IgM
para(10) = 0.05;%% decay constant of plasma Cell IgG
para(11) = 4.4e9;% production constant of IgM
para(12) = 1.2e10;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.02;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 1.2; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;





[t_new, z] = ode15s(@pathway_model_many_antibody_immune_include_plasma,[0 1000],x0,[],para,para_new,para_new_1,AA);




for i = 1:1000

strong_antibody_new(i) = interp1(t_new,z(:,251),i)+interp1(t_new,z(:,261),i)+interp1(t_new,z(:,262),i)+interp1(t_new,z(:,271),i)+interp1(t_new,z(:,272),+i)...
    +interp1(t_new,z(:,273),i)+interp1(t_new,z(:,281),i)+interp1(t_new,z(:,282),i)+interp1(t_new,z(:,283),i)+interp1(t_new,z(:,284),i)+interp1(t_new,z(:,291),i)...
    +interp1(t_new,z(:,292),i)+interp1(t_new,z(:,293),i)+interp1(t_new,z(:,294),i)+interp1(t_new,z(:,295),i);

% very_strong_antibody_new(i) = interp1(t_new,z(:,261),1000+i)+interp1(t_new,z(:,271),1000+i)+interp1(t_new,z(:,272),1000+i)...
%     +interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,282),1000+i)+interp1(t_new,z(:,283),1000+i)+interp1(t_new,z(:,291),1000+i)...
%     +interp1(t_new,z(:,292),1000+i)+interp1(t_new,z(:,293),1000+i)+interp1(t_new,z(:,294),1000+i);
% 
% super_strong_antibody_new(i) = interp1(t_new,z(:,271),1000+i)...
%     +interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,282),1000+i)+interp1(t_new,z(:,291),1000+i)...
%     +interp1(t_new,z(:,292),1000+i)+interp1(t_new,z(:,293),1000+i);
% 
% 
% super_super_strong_antibody_new(i) = interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,291),1000+i)...
%     +interp1(t_new,z(:,292),1000+i);
% 
% super_super_super_strong_antibody_new(i) = interp1(t_new,z(:,291),1000+i);
end

for i = 1:1000
    cross_over_1(iii,i) = (interp1(t_new,z(:,251),i)*(zz0(251))/x0(251) + interp1(t_new,z(:,261),i)*(zz0(261))/x0(261)+interp1(t_new,z(:,262),i)*(zz0(262))/x0(262)...
        +interp1(t_new,z(:,271),i)*(zz0(271))/x0(271)+interp1(t_new,z(:,272),i)*(zz0(272))/x0(272)+interp1(t_new,z(:,273),i)*(zz0(273))/x0(273)...
        +interp1(t_new,z(:,281),i)*(zz0(281))/x0(281)+interp1(t_new,z(:,282),i)*(zz0(282))/x0(282)+interp1(t_new,z(:,283),i)*(zz0(283))/x0(283)...
        +interp1(t_new,z(:,284),i)*(zz0(284))/x0(284)+interp1(t_new,z(:,291),i)*(zz0(291))/x0(291)+interp1(t_new,z(:,292),i)*(zz0(292))/x0(292)...
        +interp1(t_new,z(:,293),i)*(zz0(293))/x0(293)+interp1(t_new,z(:,294),i)*(zz0(294))/x0(294)+interp1(t_new,z(:,295),i)*(zz0(295))/x0(295))/strong_antibody_new(i);
    
    % cross_over_2(iii,i) = (interp1(t_new,z(:,261),1000+i)*(zz1(261))/x0(261)...
    %     +interp1(t_new,z(:,271),1000+i)*(zz1(271))/x0(271)+interp1(t_new,z(:,272),1000+i)*(zz1(272))/x0(272)...
    %     +interp1(t_new,z(:,281),1000+i)*(zz1(281))/x0(281)+interp1(t_new,z(:,282),1000+i)*(zz1(282))/x0(282)+interp1(t_new,z(:,283),1000+i)*(zz1(283))/x0(283)...
    %     +interp1(t_new,z(:,291),1000+i)*(zz1(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz1(292))/x0(292)...
    %     +interp1(t_new,z(:,293),1000+i)*(zz1(293))/x0(293)+interp1(t_new,z(:,294),1000+i)*(zz1(294))/x0(294))/very_strong_antibody_new(i);
    % 
    % cross_over_3(iii,i) = (interp1(t_new,z(:,271),1000+i)*(zz2(271))/x0(271)...
    %     +interp1(t_new,z(:,281),1000+i)*(zz2(281))/x0(281)+interp1(t_new,z(:,282),1000+i)*(zz2(282))/x0(282)...
    %     +interp1(t_new,z(:,291),1000+i)*(zz2(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz2(292))/x0(292)...
    %     +interp1(t_new,z(:,293),1000+i)*(zz2(293))/x0(293))/super_strong_antibody_new(i);
    % 
    % cross_over_4(iii,i) = (interp1(t_new,z(:,281),1000+i)*(zz3(281))/x0(281)...
    %     +interp1(t_new,z(:,291),1000+i)*(zz3(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz3(292))/x0(292))/super_super_strong_antibody_new(i);
    % 
    % 
    % cross_over_5(iii,i) = (interp1(t_new,z(:,291),1000+i)*(zz4(291))/x0(291))/super_super_super_strong_antibody_new(i);

    
    

%%
    cross_over_new_1(iii,i) = cross_over_1(iii,i)*strong_antibody_new(i)/(strong_antibody*(1-10^(-iii))+cross_over_1(iii,i)*strong_antibody_new(i));

    % cross_over_new_2(iii,i) = cross_over_2(iii,i)*very_strong_antibody_new(i)/(very_strong_antibody*(1-10^(-iii))+cross_over_2(iii,i)*very_strong_antibody_new(i));
    % 
    % cross_over_new_3(iii,i) = cross_over_3(iii,i)*super_strong_antibody_new(i)/(super_strong_antibody*(1-10^(-iii))+cross_over_3(iii,i)*super_strong_antibody_new(i));
    % 
    % cross_over_new_4(iii,i) = cross_over_4(iii,i)*super_super_strong_antibody_new(i)/(super_super_strong_antibody*(1-10^(-iii))+cross_over_4(iii,i)*super_super_strong_antibody_new(i));
    % 
    % cross_over_new_5(iii,i) = cross_over_5(iii,i)*super_super_super_strong_antibody_new(i)/(super_super_super_strong_antibody*(1-10^(-iii))+cross_over_5(iii,i)*super_super_super_strong_antibody_new(i));
end

ddd = interp1(t_new,z(:,1402),(0:1:1000));
max_virus(iii) = max(ddd);
end

% plot(t_new,z(:,1402),'linewidth',2); %% generate figure 2
   
IgG_secondary_3D = plot_figure_IgG(t_new,z); %  generate figure 4

time_interp = 0:5:1000; % 统一插值时间轴
groups = [-31:1:-13];

virus_primary   = interp1(t, y(:, 1402), time_interp, 'linear', 0);
virus_secondary = interp1(t_new, z(:, 1402), time_interp, 'linear', 0);

ratio_cross = cross_over_1(5,:);

%% 

%% 1. 创建 Nature 风格组合图
% 设置画布: 宽1000px, 高800px
figure('Units', 'pixels', 'Position', [100, 100, 1000, 800], 'Color', 'w');
std_font = 'Arial';

% =========================================================================
% Panel A: 病毒动力学对比 (Primary vs Secondary)
% =========================================================================
ax1 = subplot(2, 2, 1);
hold on;
plot(time_interp, virus_primary, 'b-', 'LineWidth', 2);
plot(time_interp, virus_secondary, 'r-', 'LineWidth', 2);
hold off;

% 美化
title('a', 'Units','normalized','Position',[-0.15, 1.05], 'FontWeight','bold', 'FontSize', 16, 'FontName', std_font);
text(0.5, 1.05, 'Viral Dynamics', 'Units','normalized','HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 12, 'FontName', std_font);
xlabel('Time (days)', 'FontName', std_font, 'FontSize', 10);
ylabel('Virus Load (Conc.)', 'FontName', std_font, 'FontSize', 10);
legend({'Primary Infection', 'Secondary (Variant)'}, 'Box', 'off', 'Location', 'Best', 'FontName', std_font);
grid on; set(gca, 'GridLineStyle', ':', 'LineWidth', 1.2, 'FontSize', 10, 'FontName', std_font, 'TickDir', 'out', 'Box', 'off');
xlim([0 200]); % 根据需要调整X轴范围展示波峰

% =========================================================================
% Panel B: 交叉反应/记忆占比随时间变化 (Figure 5)
% =========================================================================
ax2 = subplot(2, 2, 2);
plot(ratio_cross, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.5); % 紫色

% 美化
title('b', 'Units','normalized','Position',[-0.15, 1.05], 'FontWeight','bold', 'FontSize', 16, 'FontName', std_font);
text(0.5, 1.05, 'Cross-reactivity Ratio', 'Units','normalized','HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 12, 'FontName', std_font);
xlabel('Time (days)', 'FontName', std_font, 'FontSize', 10);
ylabel('Ratio (Memory / Total)', 'FontName', std_font, 'FontSize', 10);
ylim([0 1.05]); grid on; 
set(gca, 'GridLineStyle', ':', 'LineWidth', 1.2, 'FontSize', 10, 'FontName', std_font, 'TickDir', 'out', 'Box', 'off');

% =========================================================================
% Panel C: 初次感染 IgG 3D Landscape (Figure 3)
% =========================================================================
ax3 = subplot(2, 2, 3);
[X, Y] = meshgrid(time_interp,groups); % X=Time, Y=Group
surf(X, Y, IgG_primary_3D, 'EdgeColor', 'none', 'FaceColor', 'interp');

% 3D 美化
colormap(ax3, parula); shading interp;
light('Position',[0 -10 10], 'Style', 'local'); lighting gouraud; material dull;
title('c', 'Units','normalized','Position',[-0.15, 1.05], 'FontWeight','bold', 'FontSize', 16, 'FontName', std_font);
text(0.5, 1.05, 'Primary IgG Landscape', 'Units','normalized','HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 12, 'FontName', std_font);
xlabel('Time (days)', 'FontName', std_font);
ylabel('Affinity Group', 'FontName', std_font);
zlabel('Concentration', 'FontName', std_font);
view(45, 30); axis tight; xlim([0 200]); % 聚焦早期
set(gca, 'LineWidth', 1, 'FontSize', 9, 'FontName', std_font, 'TickDir', 'out', 'Box', 'off');
clb1 = colorbar; clb1.Label.String = 'Conc.';

% =========================================================================
% Panel D: 二次感染 IgG 3D Landscape (Figure 4 - Imprinting)
% =========================================================================
ax4 = subplot(2, 2, 4);
surf(X, Y, IgG_secondary_3D, 'EdgeColor', 'none', 'FaceColor', 'interp');

% 3D 美化 - 保持与 Panel C 一致的视觉风格
colormap(ax4, turbo); % 使用不同色系或保持 parula 均可，turbo对比度更高
shading interp;
light('Position',[0 -10 10], 'Style', 'local'); lighting gouraud; material dull;

title('d', 'Units','normalized','Position',[-0.15, 1.05], 'FontWeight','bold', 'FontSize', 16, 'FontName', std_font);
text(0.5, 1.05, 'Secondary IgG (Imprinting)', 'Units','normalized','HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 12, 'FontName', std_font);
xlabel('Time (days)', 'FontName', std_font);
ylabel('Affinity Group', 'FontName', std_font);
zlabel('Concentration', 'FontName', std_font);

% 统一 Z 轴坐标以便对比 (Optional: 取两图最大值)
max_z = max(max(IgG_primary_3D(:)), max(IgG_secondary_3D(:)));
zlim([0 max_z]); 

view(45, 30); axis tight;
set(gca, 'LineWidth', 1, 'FontSize', 9, 'FontName', std_font, 'TickDir', 'out', 'Box', 'off');
clb2 = colorbar; clb2.Label.String = 'Conc.';