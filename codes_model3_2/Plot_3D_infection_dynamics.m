%% --- [步骤0] 数据准备 ---
% 1. 提取感染细胞矩阵 I_matrix
% y 的第 2 列到 p.N+1 列是感染细胞 I 随时间(行)和感染龄(列)的变化




% for i = 1:400
% I_matrix(i,:) = interp1(t,y(:,i+1),(0:1:2000));
% 
% end
 I_matrix= y_new(:,2:p.N+1);



% 2. 准备网格
% 时间向量 t 已经有了
% 感染龄向量 p.a_vec 已经有了
% 生成网格以便绘制 3D图

[AgeGrid, TimeGrid] = meshgrid(p.a_vec, t_new);

%%
%% --- [方案一] 3D 演化曲面图 ---
figure('Name', '3D Dynamics of Infected Cells', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.6]);

% 绘制曲面
s = surf(AgeGrid, TimeGrid, I_matrix);

%这几行是为了让图更好看
s.EdgeColor = 'none'; % 去掉网格线，只留颜色
shading interp;       % 平滑颜色过渡
colormap(jet);        % 使用高对比度色图
colorbar;             % 显示颜色条

% 坐标轴标签
xlabel('Infection Age (a)', 'FontWeight', 'bold');
ylabel('Time (hours)', 'FontWeight', 'bold');
zlabel('Infected Cell Density', 'FontWeight', 'bold');
title('Dynamics of Infected Cells Distribution I(t, a)', 'FontSize', 14);

% 调整视角 (可以手动旋转图窗，也可以用代码固定)
view(-45, 30); 
grid on;
axis tight;

%%
% %% --- [方案二] 感染动态热力图 ---
figure('Name', 'Heatmap of Infection (Corrected)', 'Color', 'w');

% 1. 使用与 surf 完全相同的网格数据 (复用前面的变量)
% 前提是你已经运行了: [AgeGrid, TimeGrid] = meshgrid(p.a_vec, t);
% I_matrix 也是之前的 [时间 x 感染龄] 矩阵

% 2. 绘制伪彩色图 (Pseudo-color Plot)
% 注意参数顺序：我想让 X轴是时间，Y轴是感染龄
% 所以第一个参数放 TimeGrid，第二个放 AgeGrid
h = pcolor(TimeGrid, AgeGrid, I_matrix);

% 3. 关键设置：去掉网格线并平滑颜色
set(h, 'EdgeColor', 'none'); % 去掉黑色网格线
shading interp;              % 使颜色平滑过渡（和surf效果一致）

% 4. 视觉美化
colormap(jet);      % 颜色风格
colorbar;           % 显示色条
ylabel(colorbar, 'Infected Cell Density');

% 5. 坐标轴标签
xlabel('Time (hours) [X-Axis]', 'FontWeight', 'bold');
ylabel('Infection Age (a) [Y-Axis]', 'FontWeight', 'bold');
title('Dynamics Heatmap: Age Distribution over Time', 'FontSize', 12);

% 6. 确保轴范围紧凑
axis tight;

%%
%% --- [方案三] 关键变量的 3D 相轨迹 ---
figure('Name', '3D Phase Portrait', 'Color', 'w');

% 提取需要的变量
T_vals = y_new(:, 1);                % 易感细胞
Total_I = sum(I_matrix, 2);      % 总感染细胞
V_vals = y_new(:, p.N+2);            % 病毒量

% 绘制 3D 轨迹曲线
plot3(T_vals, Total_I, V_vals, 'LineWidth', 2, 'Color', 'b');
hold on;

% 标记起点（绿色圆圈）和终点（红色叉）
plot3(T_vals(1), Total_I(1), V_vals(1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(T_vals(end), Total_I(end), V_vals(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

grid on;
xlabel('Susceptible Cells (T)');
ylabel('Total Infected Cells (I)');
zlabel('Virus Load (V)');
title('Trajectory of Infection Dynamics');

% 由于病毒通常是跨数量级变化，如果图形是被压扁的，可以开启对数坐标
set(gca, 'ZScale', 'log'); 
% set(gca, 'XScale', 'log', 'YScale', 'log'); % 视情况开启
view(40, 30);