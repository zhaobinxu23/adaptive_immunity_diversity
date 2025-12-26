for i = 1:100
data_new(i,:) = interp1(t_new,y_new(:,i+503),(0:10:2000));
data_new_IgG(i,:) = interp1(t_new,y_new(:,i+703),(0:10:2000));
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
time_points = [0:10:2000];
frequencies = data_kd;
[TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序

frequencies_2 = data_kd_IgG;



% 绘制三维频率分布图
figure;
surf(Values, TimePoints, frequencies_2, 'EdgeColor', 'none');
% hold on
% surf(Values, TimePoints, frequencies_2, 'EdgeColor', 'none');


% 设置 colormap
colormap(parula);

% 添加轴标签和标题
xlabel('log(Kd) Values');
ylabel('Time');
zlabel('Frequencies');
title('3D Frequency Distribution');

% 调整视角和视点
view(30, 45);

% 添加网格线
grid on;

