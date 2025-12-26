%% --- [步骤0] 数据准备 ---
% 1. 提取感染细胞矩阵 I_matrix
% y 的第 2 列到 p.N+1 列是感染细胞 I 随时间(行)和感染龄(列)的变化
I_matrix = y_new(:, 2:p.N+1); 

% 2. 准备网格
% 时间向量 t 已经有了
% 感染龄向量 p.a_vec 已经有了
% 生成网格以便绘制 3D图
[AgeGrid, TimeGrid] = meshgrid(p.a_vec, t_new);

%%


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

