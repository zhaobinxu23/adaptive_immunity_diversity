figure('Position', [100, 100, 1600, 1200], 'Color', 'white');

for i = 1:100
    % 创建子图 - 使用更合理的布局（10×10可能太小）
    subplot(10, 10, i);
    
    % 绘制曲线，增加线宽并使用不同的颜色
    plot(t, y(:,1720+i), 'linewidth', 2.5, 'Color', [0, 0.4470, 0.7410]); % 蓝色
    
    % 设置坐标轴属性
    set(gca, 'FontSize', 8, 'FontWeight', 'bold'); % 增大字体并加粗
    grid on;
    grid minor; % 添加次要网格线
    
    % 设置标题和标签
    title(['List ', num2str(i)], 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time', 'FontSize', 9, 'FontWeight', 'bold');
    ylabel('Concentration', 'FontSize', 9, 'FontWeight', 'bold');
    
    % 优化坐标轴范围（根据数据调整）
    xlim([min(t) max(t)]); % 确保x轴覆盖所有数据
    
    % 添加边框
    box on;
end

% 调整子图间距
sgtitle('Concentration Profiles for All CD4+ T memory cell Lists', 'FontSize', 14, 'FontWeight', 'bold');

% 自动调整子图间距以避免重叠
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.8, 0.8]);

% clc
% clear
% 
% % 创建一个大图
% figure('Position', [100, 100, 1200, 800]);
% 
% for i = 1:100
%  
%     % 创建子图 - 2行5列排列
% subplot(10, 10, i);
% 
% 
% 
% plot(t,y(:,1720+i),'linewidth',2);
% hold on
% 
% 
% legend('Dynamics of memory CD4+ T cells');
% title(['CD4+ T cell list ', num2str(i)]);
% xlabel('Time');
% ylabel('Concentration');
% grid on;
%     
% hold off;
% end