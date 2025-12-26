for i = 1:100
data_new(i,:) = interp1(t,y(:,i+1600),(0:10:1000));
data_new_IgG(i,:) = interp1(t,y(:,i+1700),(0:10:1000));
end

data_k_on = zeros(10, 101);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on(row_index, :) = data_k_on(row_index, :) + data_new(i, :);
end

data_k_on_IgG = zeros(10, 101);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on_IgG(row_index, :) = data_k_on_IgG(row_index, :) + data_new_IgG(i, :);
end



data_k_off = zeros(10, 101);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off(row_index, :) = data_k_off(row_index, :) + data_new(i, :);
end

data_k_off_IgG = zeros(10, 101);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off_IgG(row_index, :) = data_k_off_IgG(row_index, :) + data_new_IgG(i, :);
end


data_kd = zeros(19, 101);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:101
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end

data_kd_IgG = zeros(19, 101);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:101
        data_kd_IgG(row_index, j) = data_kd_IgG(row_index, j) + data_new_IgG(i, j);
    end
end

values = [-26:1:-8];
time_points = [0:10:1000];
frequencies = data_kd;
[TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序

% 绘制三维频率分布图
figure;
surf(Values, TimePoints, frequencies, 'EdgeColor', 'none');


% 设置 colormap
colormap(parula);

% 添加轴标签和标题
xlabel('Time Points');
ylabel('Values');
zlabel('Frequencies');
title('3D Frequency Distribution');

% 调整视角和视点
view(30, 45);

% 添加网格线
grid on;



% for i = 1:10
% data_new_t(i,:) = interp1(t,y(:,i+8600),((0:10:1000));
% data_new_IgG_t(i,:) = interp1(t,y(:,i+8610),((0:10:1000));
% end

% for i = 1:200
% data_new_d1(i,:) = interp1(t,y(:,i),(0:10:1000));
% data_new_d2(i,:) = interp1(t,y(:,i+200),(0:10:1000));
% data_new_d3(i,:) = interp1(t,y(:,i+400),(0:10:1000));
% data_new_d4(i,:) = interp1(t,y(:,i+600),(0:10:1000));
% end
% 
% 
% for i = 1:2
% data_new_d7(i,:) = interp1(t,y(:,i+2200),(0:10:1000));
% data_new_d8(i,:) = interp1(t,y(:,i+2202),(0:10:1000));
% end