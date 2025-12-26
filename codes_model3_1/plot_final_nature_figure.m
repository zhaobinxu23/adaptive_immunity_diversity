function plot_final_nature_figure(time_vec, virus_data, ...
                                  igm_elisa, igg_elisa, total_elisa, ...
                                  affinity_matrix)
% PLOT_FINAL_NATURE_FIGURE 生成包含病毒、ELISA和Kd热图的顶级图表
%
% 输入:
%   time_vec: 时间向量 (e.g., 0:200)
%   virus_data: 病毒载量向量
%   igm_elisa: IgM ELISA 读数向量
%   igg_elisa: IgG ELISA 读数向量
%   total_elisa: 总抗体 ELISA 读数向量
%   affinity_matrix: 19 x N 的亲和力分布矩阵 (frequencies_2)
%                    Row 1 = Kd -31, Row 19 = Kd -13

    % 创建画布
    figure('Position', [100, 50, 850, 900], 'Color', 'w');
    layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'normal');

    % ==========================================================
    % Panel A: 病毒载量 (左轴) & ELISA (右轴)
    % ==========================================================
    ax1 = nexttile;
    
    % --- 左轴：病毒载量 (红色系) ---
    yyaxis left
    % 使用半对数坐标展示病毒的大范围变化
    h_virus = semilogy(time_vec, virus_data, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
    ylabel('Viral Load (copies/mL)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.8 0.2 0.2]);
    set(gca, 'YColor', [0.8 0.2 0.2], 'FontSize', 11, 'LineWidth', 1.2);
    % 调整刻度标签格式
    ax1.YAxis(1).Exponent = 0; 
    
    % --- 右轴：ELISA (蓝色/黑色系) ---
    yyaxis right
    hold on;
    % IgM (青色虚线)
    h_igm = plot(time_vec, igm_elisa, '--', 'Color', [0.2 0.7 0.8], 'LineWidth', 2);
    % IgG (蓝色实线)
    h_igg = plot(time_vec, igg_elisa, '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 2.5);
    % Total (黑色点线，作为总量参考)
    h_total = plot(time_vec, total_elisa, ':', 'Color', 'k', 'LineWidth', 2);
    
    ylabel('Antibody ELISA (OD 450nm)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    set(gca, 'YColor', 'k');
    
    % --- Panel A 美化 ---
    title('A. Viral Clearance & Antibody Response (ELISA)', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    grid off; box off;
    xticklabels([]); % 移除X轴标签，保持整洁
    xlim([0, max(time_vec)]);
    
    % 图例
    legend([h_virus, h_igm, h_igg, h_total], ...
           {'Viral Load', 'IgM ELISA', 'IgG ELISA', 'Total Ab ELISA'}, ...
           'Location', 'northeast', 'Box', 'off', 'FontSize', 10);

    % ==========================================================
    % Panel B: 亲和力成熟热图 (Log Kd: -31 -> -13)
    % ==========================================================
    ax2 = nexttile;
    
    % --- 绘制热图 ---
    % 使用 Log10 压缩频率值以显示细节，加微小量防止 log(0)
    h_map = imagesc(time_vec, 1:19, log10(affinity_matrix + 1e-8));
    
    % --- 关键：自定义 Y 轴 (Log Kd) ---
    % 矩阵有19行。Row 1 对应 -31，Row 19 对应 -13。
    % 我们需要设置 Y 轴刻度，并标记正确的数值。
    
    set(gca, 'YDir', 'normal'); % 确保 Row 1 (Kd -31) 在底部，Row 19 (Kd -13) 在顶部
    
    % 设置 3 个主要刻度：底部、中间、顶部
    ticks_y = [1, 10, 19];
    labels_y = {'-31', '-22', '-13'};
    
    yticks(ticks_y);
    yticklabels(labels_y);
    
    ylabel('Log (K_d)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time post-infection (Days)', 'FontSize', 12, 'FontWeight', 'bold');
    
    % --- 热图美化 ---
    colormap(ax2, 'parula'); % 或 'jet', 'turbo'
    c = colorbar;
    c.Label.String = 'Log_{10}(Frequency)';
    c.Label.FontSize = 11;
    
    title('B. Affinity Maturation Landscape (Log K_d Evolution)', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    
    % --- 叠加“高亲和力”趋势线 (可选) ---
    hold on;
    % 计算亲和力重心 (Weighted Center of Mass)
    % 这条线展示平均 Kd 随时间如何从 -31 移向 -13
    mean_kd_idx = zeros(size(time_vec));
    for t = 1:length(time_vec)
        col = affinity_matrix(:, t);
        if sum(col) > 0
            mean_kd_idx(t) = sum((1:19)' .* col) / sum(col);
        else
            mean_kd_idx(t) = 1;
        end
    end
    % 转化为实际 Kd 值用于 Check (不需要画出来，只需要画索引位置)
    plot(time_vec, mean_kd_idx, 'w--', 'LineWidth', 1.5);
    text(time_vec(end)*0.8, mean_kd_idx(end)+1.5, 'Mean Affinity', 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

    % 链接坐标轴
    linkaxes([ax1, ax2], 'x');
    xlim([0, 200]);
end