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


% plot(t,y(:,1402),'linewidth',2);
% hold on


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

for iii = 1:10
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
 x0(1402) = 1e14;





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
para(19) = -0.02; % virus replication constant
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





[t_new, z] = ode15s(@pathway_model_many_antibody_immune_include_plasma,[1000 2000],x0,[],para,para_new,para_new_1,AA);

% plot(t, y(:,1402),'linewidth',2);
% hold on
% plot(t_new,y_new(:,1402),'linewidth',2);
% hold on

for i = 1:1000

strong_antibody_new(i) = interp1(t_new,z(:,251),1000+i)+interp1(t_new,z(:,261),1000+i)+interp1(t_new,z(:,262),1000+i)+interp1(t_new,z(:,271),1000+i)+interp1(t_new,z(:,272),1000+i)...
    +interp1(t_new,z(:,273),1000+i)+interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,282),1000+i)+interp1(t_new,z(:,283),1000+i)+interp1(t_new,z(:,284),1000+i)+interp1(t_new,z(:,291),1000+i)...
    +interp1(t_new,z(:,292),1000+i)+interp1(t_new,z(:,293),1000+i)+interp1(t_new,z(:,294),1000+i)+interp1(t_new,z(:,295),1000+i);

very_strong_antibody_new(i) = interp1(t_new,z(:,261),1000+i)+interp1(t_new,z(:,271),1000+i)+interp1(t_new,z(:,272),1000+i)...
    +interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,282),1000+i)+interp1(t_new,z(:,283),1000+i)+interp1(t_new,z(:,291),1000+i)...
    +interp1(t_new,z(:,292),1000+i)+interp1(t_new,z(:,293),1000+i)+interp1(t_new,z(:,294),1000+i);

super_strong_antibody_new(i) = interp1(t_new,z(:,271),1000+i)...
    +interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,282),1000+i)+interp1(t_new,z(:,291),1000+i)...
    +interp1(t_new,z(:,292),1000+i)+interp1(t_new,z(:,293),1000+i);


super_super_strong_antibody_new(i) = interp1(t_new,z(:,281),1000+i)+interp1(t_new,z(:,291),1000+i)...
    +interp1(t_new,z(:,292),1000+i);

super_super_super_strong_antibody_new(i) = interp1(t_new,z(:,291),1000+i);
end

for i = 1:1000
    cross_over_1(iii,i) = (interp1(t_new,z(:,251),1000+i)*(zz0(251))/x0(251) + interp1(t_new,z(:,261),1000+i)*(zz0(261))/x0(261)+interp1(t_new,z(:,262),1000+i)*(zz0(262))/x0(262)...
        +interp1(t_new,z(:,271),1000+i)*(zz0(271))/x0(271)+interp1(t_new,z(:,272),1000+i)*(zz0(272))/x0(272)+interp1(t_new,z(:,273),1000+i)*(zz0(273))/x0(273)...
        +interp1(t_new,z(:,281),1000+i)*(zz0(281))/x0(281)+interp1(t_new,z(:,282),1000+i)*(zz0(282))/x0(282)+interp1(t_new,z(:,283),1000+i)*(zz0(283))/x0(283)...
        +interp1(t_new,z(:,284),1000+i)*(zz0(284))/x0(284)+interp1(t_new,z(:,291),1000+i)*(zz0(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz0(292))/x0(292)...
        +interp1(t_new,z(:,293),1000+i)*(zz0(293))/x0(293)+interp1(t_new,z(:,294),1000+i)*(zz0(294))/x0(294)+interp1(t_new,z(:,295),1000+i)*(zz0(295))/x0(295))/strong_antibody_new(i);
    
    cross_over_2(iii,i) = (interp1(t_new,z(:,261),1000+i)*(zz1(261))/x0(261)...
        +interp1(t_new,z(:,271),1000+i)*(zz1(271))/x0(271)+interp1(t_new,z(:,272),1000+i)*(zz1(272))/x0(272)...
        +interp1(t_new,z(:,281),1000+i)*(zz1(281))/x0(281)+interp1(t_new,z(:,282),1000+i)*(zz1(282))/x0(282)+interp1(t_new,z(:,283),1000+i)*(zz1(283))/x0(283)...
        +interp1(t_new,z(:,291),1000+i)*(zz1(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz1(292))/x0(292)...
        +interp1(t_new,z(:,293),1000+i)*(zz1(293))/x0(293)+interp1(t_new,z(:,294),1000+i)*(zz1(294))/x0(294))/very_strong_antibody_new(i);

    cross_over_3(iii,i) = (interp1(t_new,z(:,271),1000+i)*(zz2(271))/x0(271)...
        +interp1(t_new,z(:,281),1000+i)*(zz2(281))/x0(281)+interp1(t_new,z(:,282),1000+i)*(zz2(282))/x0(282)...
        +interp1(t_new,z(:,291),1000+i)*(zz2(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz2(292))/x0(292)...
        +interp1(t_new,z(:,293),1000+i)*(zz2(293))/x0(293))/super_strong_antibody_new(i);

    cross_over_4(iii,i) = (interp1(t_new,z(:,281),1000+i)*(zz3(281))/x0(281)...
        +interp1(t_new,z(:,291),1000+i)*(zz3(291))/x0(291)+interp1(t_new,z(:,292),1000+i)*(zz3(292))/x0(292))/super_super_strong_antibody_new(i);


    cross_over_5(iii,i) = (interp1(t_new,z(:,291),1000+i)*(zz4(291))/x0(291))/super_super_super_strong_antibody_new(i);

    
    

%%
    cross_over_new_1(iii,i) = cross_over_1(iii,i)*strong_antibody_new(i)/(strong_antibody*(1-10^(-iii))+cross_over_1(iii,i)*strong_antibody_new(i));

    cross_over_new_2(iii,i) = cross_over_2(iii,i)*very_strong_antibody_new(i)/(very_strong_antibody*(1-10^(-iii))+cross_over_2(iii,i)*very_strong_antibody_new(i));

    cross_over_new_3(iii,i) = cross_over_3(iii,i)*super_strong_antibody_new(i)/(super_strong_antibody*(1-10^(-iii))+cross_over_3(iii,i)*super_strong_antibody_new(i));

    cross_over_new_4(iii,i) = cross_over_4(iii,i)*super_super_strong_antibody_new(i)/(super_super_strong_antibody*(1-10^(-iii))+cross_over_4(iii,i)*super_super_strong_antibody_new(i));

    cross_over_new_5(iii,i) = cross_over_5(iii,i)*super_super_super_strong_antibody_new(i)/(super_super_super_strong_antibody*(1-10^(-iii))+cross_over_5(iii,i)*super_super_super_strong_antibody_new(i));
end

ddd = interp1(t_new,z(:,1402),(1000:1:2000));
max_virus(iii) = max(ddd);
end



plot(cross_over_new_1(1,:));
hold on
plot(cross_over_new_1(2,:));
hold on
plot(cross_over_new_1(3,:));
hold on
plot(cross_over_new_1(4,:));
hold on
plot(cross_over_new_1(5,:));
hold on
plot(cross_over_new_1(6,:));
hold on
plot(cross_over_new_1(7,:));
hold on
plot(cross_over_new_1(8,:));
hold on
plot(cross_over_new_1(9,:));
hold on
plot(cross_over_new_1(10,:));
hold on
% 
% 
% 
% %%
% 
% plot(cross_over_1(1,:));
% hold on
% plot(cross_over_1(2,:));
% hold on
% plot(cross_over_1(3,:));
% hold on
% plot(cross_over_1(4,:));
% hold on
% 
% plot(cross_over_1(5,:));
% hold on
% 
% plot(cross_over_1(6,:));
% hold on
% 
% plot(cross_over_1(7,:));
% hold on
% plot(cross_over_1(8,:));
% hold on
% plot(cross_over_1(9,:));
% hold on
% plot(cross_over_1(10,:));
% hold on