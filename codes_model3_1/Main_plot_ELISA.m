% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
% clc
% clear
%% 下面这张图是产生初次感染的elisa结果的代码
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





[t, y] = ode15s(@pathway_model_many_antibody_immune_include_plasma,[0 200],x0,[],para,para_new,para_new_1,AA);


%% 转化生成结合系数的frequencies_2和frequencies数据

for i = 1:100
data_new(i,:) = interp1(t,y(:,i),(0:1:200));%% represent IgM-BCR 
data_new_IgG(i,:) = interp1(t,y(:,i+200),(0:1:200));%% represent IgG-BCR 
end

data_k_on = zeros(10, 201);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on(row_index, :) = data_k_on(row_index, :) + data_new(i, :);
end

data_k_on_IgG = zeros(10, 201);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on_IgG(row_index, :) = data_k_on_IgG(row_index, :) + data_new_IgG(i, :);
end



data_k_off = zeros(10, 201);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off(row_index, :) = data_k_off(row_index, :) + data_new(i, :);
end

data_k_off_IgG = zeros(10, 201);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off_IgG(row_index, :) = data_k_off_IgG(row_index, :) + data_new_IgG(i, :);
end


data_kd = zeros(19, 201);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:201
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end

data_kd_IgG = zeros(19, 201);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:201
        data_kd_IgG(row_index, j) = data_kd_IgG(row_index, j) + data_new_IgG(i, j);
    end
end

values = [-31:1:-13];
time_points = [0:1:200];
frequencies = data_kd;
[TimePoints, Values] = meshgrid(time_points,values); % 注意这里的顺序

frequencies_2 = data_kd_IgG;

Virus_Vector = interp1(t,y(:,1402),(0:1:200));


for i = 1:100
data_new_IgM(i,:) = interp1(t,y(:,i+100),(0:1:200));
data_new_IgG(i,:) = interp1(t,y(:,i+300),(0:1:200));
end

clear x0




for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 101:200
    x0(i) = data_new_IgM(i-100,k)/1e4;
end

for i = 301:400
    x0(i) = data_new_IgG(i-300,k)/1e4;
end

x0(1402) = 1e16;%% virus



para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
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


[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_elisa(k) = 1e16-interp1(t,zz(:,1402),2);

end


for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 101:200
    x0(i) = data_new_IgM(i-100,k)/1e4;
end



x0(1402) = 1e16;%% virus






para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
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



[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_IgM(k) = 1e16-interp1(t,zz(:,1402),2);

end


%% 
for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 301:400
    x0(i) = data_new_IgG(i-300,k)/1e4;
end



x0(1402) = 1e16;%% virus






para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
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



[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_IgG(k) = 1e16-interp1(t,zz(:,1402),2);

end



t_interp = 0:200; 
len = length(t_interp);

% 1. 提取病毒载量 (第1402列) 并插值
Virus_Vector = interp1(t, y(:, 1402), t_interp);

% 2. 提取 IgM 和 IgG 的原始数据 (用于 ELISA 计算)
% 注意：请根据您的模型确认 IgM 和 IgG 在 y 中的列索引
% 根据附件推测：IgM 可能在 101-200 或 1201-1300，此处以 1201-1300 为例
% 如果您的索引是 101-200，请手动修改下面的 range
idx_IgM_start = 101; 
idx_IgG_start = 301;

raw_IgM = zeros(100, len);
raw_IgG = zeros(100, len);

for i = 1:100
    % 提取并插值，防止维度不匹配
    col_m = idx_IgM_start + i - 1;
    col_g = idx_IgG_start + i - 1;
    raw_IgM(i, :) = interp1(t, y(:, col_m), t_interp);
    raw_IgG(i, :) = interp1(t, y(:, col_g), t_interp);
end

% 3. 模拟计算 ELISA 数值
% (注：真实的 ELISA 计算在您的文档中似乎需要再次调用 ode15s。
% 为了快速作图，这里使用 近似公式 或 请您替换为您已经算好的 out_IgM_elisa 变量)

% 如果您手头已经有 out_IgM_elisa, out_IgG_elisa, out_overall_elisa
% 请直接跳过下面一段。

% --- 模拟 ELISA 数据 (仅作演示，如果您有真实数据请替换) ---
% 假设 ELISA 数值与抗体浓度成正比，并受病毒中和影响

% ---------------------------------------------------------

% 此时您应该有以下变量准备好传入绘图函数：
% 1. t_interp
% 2. Virus_Vector
% 3. IgM_ELISA_Vector
% 4. IgG_ELISA_Vector
% 5. Total_ELISA_Vector
% 6. frequencies_2 (19x201 矩阵)

% 运行绘图函数
plot_final_nature_figure(t_interp, Virus_Vector, ...
                         out_overall_IgM, out_overall_IgG, out_overall_elisa, ...
                         frequencies_2);


