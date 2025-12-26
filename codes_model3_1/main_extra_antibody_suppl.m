clc
clear

% 创建一个大图
figure('Position', [100, 100, 1200, 800]);

for i = 1:9
    % 创建子图 - 2行5列排列
    subplot(3, 3, i);
    
    x0(1) = 10; %% antibody
    x0(2) = 0; %% injected antibody
    x0(3) = 1; %% virus 
    x0(4) = 0; %% antibody_virus complex
    x0(5) = 0;  %% injected antibody-virus complex

    para(1) = 1e-5;
    para(2) = 1e-14;
    para(3) = 1.5;% feedback constant
    para(4) = 0.5;% virus-antibody complex decay constant
    para(5) = 0.5; % virus proliferation constant
    para(6) = 1e-5;
    para(7) = 1e-14;
    para(8) = 0.05; % monoclonal antibody decay constant

    [t,y]=ode15s(@pathway_model_extra_antibody,[0 28],x0,[],para);

    % 绘制第一段曲线

    plot(t,(y(:,3)),'linewidth',2);
    hold on

    x0(1) = interp1(t,y(:,1),28);
    x0(2) = 0.5*i*1e6;
    x0(3) = interp1(t,y(:,3),28);
    x0(4) = interp1(t,y(:,4),28);
    x0(5) = interp1(t,y(:,5),28);

    [t_new,y_new]=ode15s(@pathway_model_extra_antibody,[28 150],x0,[],para);
    
    % 绘制第二段曲线
    plot(t_new,(y_new(:,3)),'linewidth',2);

    
    % 添加图例和标题
    legend('Virus dynamics before therapy',  'Virus dynamics after therapy');
    title(['Iteration ', num2str(i)]);
    xlabel('Time');
    ylabel('Concentration');
    grid on;
    
    hold off;
end

% 添加总标题
sgtitle(' Model Simulation Results of Different Injection Dose');