%% 1. 数据准备与排序设置 (Data Preparation & Sorting)

% 假设你的原始数据矩阵如下 (100行 x 时间点列)
% 请确保这些变量在工作区中已加载

for i = 1:100
data_IgM(i,:) = interp1(t_new,y_new(:,i+503),(0:1:2000));
data_IgG(i,:) = interp1(t_new,y_new(:,i+703),(0:1:2000));

data_IgM_V(i,:) = interp1(t_new,y_new(:,i+403+900),(0:1:2000));
data_IgG_V(i,:) = interp1(t_new,y_new(:,i+403+1100),(0:1:2000));
end


%% 1. 数据准备 (Data Preparation)
% 假设工作区已有数据: data_IgM_V, data_IgG_V, data_IgM, data_IgG (大小 100 x Time)
% 如果没有，请先加载你的数据。

% --- 【用户填写区域】定义分组映射关系 ---
% 格式：{ 分组名称(数字), [包含的抗体ID列表] }
% 这一步决定了热图从上到下的排列顺序
group_mapping = {
    -13, [91];          % Group -13
    -14, [92, 81];      % Group -14
    -15, [93,82,71];            
    -16, [94,83,72,61];     
    -17, [95,84,73,62,51];     
    -18, [96,85,74,63,52,41];     
    -19, [97,86,75,64,53,42,31];     
    -20, [98,87,76,65,54,43,32,21];     
    -21, [99,88,77,66,55,44,33,22,11];     
    -22, [100,89,78,67,56,45,34,23,12,1];     
    -23, [90,79,68,57,46,35,24,13,2];     
    -24, [80,69,58,47,36,25,14,3];     
    -25, [70,59,48,37,26,15,4];     
    -26, [60,49,38,27,16,5];     
    -27, [50,39,28,17,6];     
    -28, [40,29,18,7];     
    -29, [30,19,8]; 
    -30, [20,9]; 
    -31, [10];          % Group -31
};
% --- 自动生成排序索引 ---
new_order = [];       
group_struct = [];    
current_row = 1;

for i = 1:size(group_mapping, 1)
    g_name = group_mapping{i, 1};
    g_ids = group_mapping{i, 2};
    
    if ~isempty(g_ids)
        new_order = [new_order, g_ids]; 
        % 记录绘图所需的Y轴位置信息
        s.name = num2str(g_name);
        s.start_row = current_row;      
        s.end_row = current_row + length(g_ids) - 1;
        group_struct = [group_struct, s];
        current_row = current_row + length(g_ids);
    end
end

% 补全剩余ID，防止报错
if length(new_order) ~= 100
    remaining = setdiff(1:100, new_order);
    new_order = [new_order, remaining];
end

% --- 根据新顺序重排数据并取对数 ---
plot_data = { ...
    log10(data_IgM_V(new_order, :)), 'IgM-Virus Complex', 'parula';
    log10(data_IgG_V(new_order, :)), 'IgG-Virus Complex', 'parula';
    log10(data_IgG(new_order, :)),   'Total IgG',         'parula';
    log10(data_IgM(new_order, :)),   'Total IgM',         'parula'
};

%% 2. 绘图 (Plotting - Nature Style)

% 设置画布：Nature 全幅宽度约为 18cm，我们设大一点以便看清，但比例要对
fig = figure('Units', 'centimeters', 'Position', [1, 1, 30, 16], 'Color', 'w');

% 布局：左侧 Padding 设为 loose 以容纳 Affinity 标签
% TileSpacing 改为 'none' 让图紧挨着，这是热图拼接的常见做法
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'loose'); 

% 调整 TileLayout 的外边距 (Left, Bottom, Right, Top)
% 这里的 'LabelOffset' 是逻辑上的，实际利用 Padding
% 我们稍后会在 ax1 上手动处理边距

axes_handles = [];
max_time = size(plot_data{1,1}, 2);

% --- 循环绘制 4 个子图 ---
for i = 1:4
    ax = nexttile;
    axes_handles = [axes_handles, ax];
    
    % 绘制热图
    imagesc(ax, plot_data{i,1});
    
    % 基础样式
    apply_nature_style(ax, plot_data{i,2}, plot_data{i,3}, max_time);
    
    % 仅在第1个图显示左侧 Y 轴标签
    if i == 1
        % 设置 Y 轴显示重排后的抗体 ID
        set(ax, 'YTick', 1:100, 'YTickLabel', new_order);
        ylabel(ax, 'Antibody ID', 'FontWeight', 'bold', 'FontSize', 9);
        
        % 【关键】绘制左侧的分组括号
        draw_left_brackets(ax, group_struct, max_time);
    else
        set(ax, 'YTickLabel', []); % 隐藏其他图的 Y 轴文字
        ylabel(ax, '');
    end
end

% --- 绘制贯穿所有子图的白色分割线 ---
line_positions = [group_struct.end_row] + 0.5; 
line_positions(end) = []; % 去掉底部的线

for i = 1:length(line_positions)
    y_pos = line_positions(i);
    for k = 1:4
        % 绘制白色细线，透明度设为 0.6，视觉上分割 Group
        yline(axes_handles(k), y_pos, '-', 'Color', [1 1 1], 'LineWidth', 0.8, 'Alpha', 0.7);
    end
end

% 添加总标题 (可选)
title(t, 'Immune Response Dynamics after therapeutic vaccine injection (dosage = 1*10^1^4)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

%% 3. 辅助函数 (Helper Functions)

function apply_nature_style(ax, title_str, cmap_name, max_x)
    % 配色
    colormap(ax, cmap_name);
    h = colorbar(ax);
    h.Label.String = 'Log_{10}(Conc.)';
    h.Label.FontSize = 8;
    h.Label.FontName = 'Arial';
    % 调整 colorbar 宽度，使其更细致
    ax.ColorScale = 'log'; % 只要数据已经是log了就不需要这个，但保持线性映射
    
    % 标题与轴标签
    title(ax, title_str, 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Arial');
    xlabel(ax, 'Time', 'FontSize', 8, 'FontName', 'Arial');
    
    % 坐标轴美化
    set(ax, 'FontName', 'Arial', 'FontSize', 7); % 坐标轴字体 7pt
    set(ax, 'LineWidth', 0.75);                  % 线宽适中
    set(ax, 'TickDir', 'out');                   % 刻度向外
    set(ax, 'Box', 'off');                       % 去掉多余边框
    set(ax, 'YDir', 'reverse');                  % 确保第1行在最上面
    
    xlim([0.5, max_x+0.5]);
    ylim([0.5, 100.5]);
    
    clim([0 13]); % 统一色阶范围
end

function draw_left_brackets(ax, g_struct, max_x)
    % 允许在坐标轴外绘图
    set(ax, 'Clipping', 'off'); 
    
    % --- 位置参数调整 (解决看不清的问题) ---
    % 定义偏移量 (单位：X轴数据单位)
    % 假设 X 轴是 0-400。
    % bracket_x_offset: 括号竖线距离 Y 轴的距离
    % text_x_offset: 文字距离 Y 轴的距离
    
    bracket_x_offset = -0.10 * max_x;  % 稍微拉近一点，原先可能是 -0.15
    text_x_offset    = -0.12 * max_x;  % 文字紧跟括号
    tick_len         = 0.02 * max_x;   % 括号短横线的长度
    
    for i = 1:length(g_struct)
        gs = g_struct(i);
        y_top = gs.start_row - 0.5;
        y_bottom = gs.end_row + 0.5;
        y_center = (y_top + y_bottom) / 2;
        
        % 1. 绘制括号线条 (黑色，能够看清)
        % 形状：[横线-竖线-横线]
        line_x = [bracket_x_offset+tick_len, bracket_x_offset, bracket_x_offset, bracket_x_offset+tick_len];
        line_y = [y_top, y_top, y_bottom, y_bottom];
        
        % 如果只有一行，画一个小短线即可
        if gs.start_row == gs.end_row
             line(ax, [bracket_x_offset, bracket_x_offset+tick_len], [y_center, y_center], ...
                 'Color', 'k', 'LineWidth', 1, 'Clipping', 'off');
        else
             line(ax, line_x, line_y, 'Color', 'k', 'LineWidth', 1, 'Clipping', 'off');
        end
        
        % 2. 绘制 Affinity Group 数值
        text(ax, text_x_offset, y_center, gs.name, ...
            'HorizontalAlignment', 'right', ... % 右对齐，让文字末尾对其括号
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontSize', 8, ...       % 字体不需要太大，清晰即可
            'FontWeight', 'bold', ...
            'Color', 'k');
    end
    
    % 在顶部添加列名 "Affinity\nGroup" (\n换行省空间)
    text(ax, text_x_offset, -2, sprintf('Affinity\nGroup'), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 9, 'FontWeight', 'bold', 'FontName', 'Arial');
end