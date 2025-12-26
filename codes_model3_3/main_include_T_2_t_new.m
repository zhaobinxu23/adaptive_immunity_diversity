% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
clc
clear

x0 = zeros(2606,1);

%% two B cell distributions

prob_C(1) = 1e-5;
prob_C(2) = 1-1e-5;

%%  IgM distribution
mu = -16.5; % 均值
sigma = 0.8; % 标准差


prob_A(1) = normcdf(-20.5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(-20.5+i-1, mu, sigma) - normcdf(-20.5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end


mu = 1.5; % 均值
sigma = 1; % 标准差


prob_B(1) = normcdf(-2.5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_B(i) = normcdf(-2.5+i-1, mu, sigma) - normcdf(-2.5+i-2, mu, sigma);
end
for i = 6:10
prob_B(i) = prob_B(11-i);
end

for i = 1:10
    for j = 1:10
        final_AA(10*(i-1)+j) =  prob_A(i)*prob_B(j);
    end
end

% %%  IgM distribution
% mu = 9; % 均值
% sigma = 1; % 标准差
% 
% 
% prob_A(1) = normcdf(mu-4, mu, sigma);
% % 计算累积分布概率
% for i = 2:5
% prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
% end
% for i = 6:10
% prob_A(i) = prob_A(11-i);
% end

A_1 = 1e15;
A_2 = 1e15;
T_1 = 1e10;
T_2 = 1e10;
E = 1e15;
D_1 = 1e8;
D_2 = 1e8;
D_3 = 1e8;
D_4 = 1e8;
C_1 = 3.6e13;
C_2 = 3.6e13;
F_1 = 1e13;
F_2 = 1e13;


% A_1 = 1e15;
% A_2 = 1e15;
% T_1 = 1e10;
% T_2 = 1e10;
% E = 1e15;
% D_1 = 1e8;
% D_2 = 1e8;
% D_3 = 1e8;
% D_4 = 1e8;
% C_1 = 1e12;
% C_2 = 1e12;
% F_1 = 5e12-2e8;
% F_2 = 5e12-2e8;



for i = 1:100
    for j = 1:2
        x0(2*(i-1)+j) = D_1*prob_C(j)*final_AA(i);
    end
end
for i = 1:100
    for j = 1:2
        x0(2*(i-1)+j+200) = D_2*prob_C(j)*final_AA(i);
    end
end
for i = 1:100
    for j = 1:2
        x0(2*(i-1)+j+400) = D_3*prob_C(j)*final_AA(i);
    end
end
for i = 1:100
    for j = 1:2
        x0(2*(i-1)+j+600) = D_4*prob_C(j)*final_AA(i);
    end
end




for i = 1:100
    x0(i+1600) = A_1*final_AA(i);
end

for i = 1:100
     x0(i+1700) = A_2*final_AA(i);
end

for i = 1:100
    x0(i+1800) = C_1*final_AA(i);
end

for i = 1:100
    x0(i+1900) = C_2*final_AA(i);
end


for i = 1:100
    x0(i+2000) = F_1*final_AA(i);
end

for i = 1:100
    x0(i+2100) = F_2*final_AA(i);
end


for i = 1:400
    x0(i+2200) = 0;
end


for i = 1:2
    x0(i+2600) = T_1*prob_C(i);
end

for i = 1:2
    x0(i+2602) = T_2*prob_C(i);
end
 x0(2605) = E;
 x0(2606) = 10;




para(1) = 1e-4; % T and C binding constant
para(2) = 1e6-0.6;  % 5e5-20.1
para(3) = 1e-16; % 1e-8; % A and E binding constant 
para(4) = 10/3.6-1;% 99
para(5) = 0.4e13; % pi1 0.019*1e15
para(6) = 0.8e8;% pi3 1.8e8 
para(7) = 0.008; % decay constant for A1  0.038
para(8) = 0.004;% decay constant for A2 0.019
para(9) = 0.016; % decay constant for T1   0.036
para(10) = 0.008;%% decay constant for T2  0.018
para(11) = 0.5;% decay constant for C
para(12) = 0.5;% decay constant for D


para(13) = 0.6; % pi1'  0.1 
para(14) = 0;% pi3'

para(15) = 0; % decay constant for E





para(16) = 7.2e13;% pi2  2e12
para(17) = 0.8;% virus_replication constant
para(18) = 1e5;% N
para(19) = 2;% regenerate constant for B cell
para(20) = 1;%% regeneration constant for T cell
para(21) = 2;% amplication constant for virus antigen

for j = 1:10
       para_new_1(10*(1-1)+j) = 10^(-21);
       para_new_1(10*(2-1)+j) = 10^(-20);
       para_new_1(10*(3-1)+j) = 10^(-19);
       para_new_1(10*(4-1)+j) = 10^(-18);
       para_new_1(10*(5-1)+j) = 10^(-17);
       para_new_1(10*(6-1)+j) = 10^(-16);
       para_new_1(10*(7-1)+j) = 10^(-15);
       para_new_1(10*(8-1)+j) = 10^(-14);
       para_new_1(10*(9-1)+j) = 10^(-13);
       para_new_1(10*(10-1)+j) = 10^(-12);
end



       
for j = 1:10
       para_new_2(10*(j-1)+1) = 1e-3;
       para_new_2(10*(j-1)+2) = 1e-2;
       para_new_2(10*(j-1)+3) = 1e-1;
       para_new_2(10*(j-1)+4) = 1e0;
       para_new_2(10*(j-1)+5) = 1e1;
       para_new_2(10*(j-1)+6) = 1e2;
       para_new_2(10*(j-1)+7) = 1e3;
       para_new_2(10*(j-1)+8) = 1e4;
       para_new_2(10*(j-1)+9) = 1e5;
       para_new_2(10*(j-1)+10) = 1e6;
end

para_new_3(1) = 1e-3;
para_new_3(2) = 1e-8;

para_new_4(1) = 1e-3;
para_new_4(2) = 1e6;

[t, y]=ode15s(@pathway_model_many_antibody_immune_include_T_complicated_2t_new,[0 1000],x0,[],para,para_new_1,para_new_2,para_new_3,para_new_4,final_AA,prob_C);

plot(t,y(:,2606),'linewidth',2);
hold on




%%


% x0 = zeros(2606,1);
% 
% %% two B cell distributions
% 
% prob_C(1) = 1e-5;
% prob_C(2) = 1-1e-5;
% 
% %%  IgM distribution
% mu = -16.5; % 均值
% sigma = 0.8; % 标准差
% 
% 
% prob_A(1) = normcdf(-20.5, mu, sigma);
% % 计算累积分布概率
% for i = 2:5
% prob_A(i) = normcdf(-20.5+i-1, mu, sigma) - normcdf(-20.5+i-2, mu, sigma);
% end
% for i = 6:10
% prob_A(i) = prob_A(11-i);
% end
% 
% 
% mu = 1.5; % 均值
% sigma = 1; % 标准差
% 
% 
% prob_B(1) = normcdf(-2.5, mu, sigma);
% % 计算累积分布概率
% for i = 2:5
% prob_B(i) = normcdf(-2.5+i-1, mu, sigma) - normcdf(-2.5+i-2, mu, sigma);
% end
% for i = 6:10
% prob_B(i) = prob_B(11-i);
% end
% 
% for i = 1:10
%     for j = 1:10
%         final_AA(10*(i-1)+j) =  prob_A(i)*prob_B(j);
%     end
% end
% 
% % %%  IgM distribution
% % mu = 9; % 均值
% % sigma = 1; % 标准差
% % 
% % 
% % prob_A(1) = normcdf(mu-4, mu, sigma);
% % % 计算累积分布概率
% % for i = 2:5
% % prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
% % end
% % for i = 6:10
% % prob_A(i) = prob_A(11-i);
% % end
% 
% A_1 = 1e15;
% A_2 = 1e15;
% T_1 = 1e10;
% T_2 = 1e10;
% E = 1e15;
% D_1 = 1e8;
% D_2 = 1e8;
% D_3 = 1e8;
% D_4 = 1e8;
% C_1 = 3.6e13;
% C_2 = 3.6e13;
% F_1 = 1e13;
% F_2 = 1e13;
% 
% 
% % A_1 = 1e15;
% % A_2 = 1e15;
% % T_1 = 1e10;
% % T_2 = 1e10;
% % E = 1e15;
% % D_1 = 1e8;
% % D_2 = 1e8;
% % D_3 = 1e8;
% % D_4 = 1e8;
% % C_1 = 1e12;
% % C_2 = 1e12;
% % F_1 = 5e12-2e8;
% % F_2 = 5e12-2e8;
% 
% 
% 
% for i = 1:100
%     for j = 1:2
%         x0(2*(i-1)+j) = D_1*prob_C(j)*final_AA(i);
%     end
% end
% for i = 1:100
%     for j = 1:2
%         x0(2*(i-1)+j+200) = D_2*prob_C(j)*final_AA(i);
%     end
% end
% for i = 1:100
%     for j = 1:2
%         x0(2*(i-1)+j+400) = D_3*prob_C(j)*final_AA(i);
%     end
% end
% for i = 1:100
%     for j = 1:2
%         x0(2*(i-1)+j+600) = D_4*prob_C(j)*final_AA(i);
%     end
% end
% 
% 
% 
% 
% for i = 1:100
%     x0(i+1600) = interp1(t,y(:,i+1600),1000);
% end
% 
% for i = 1:100
%      x0(i+1700) = interp1(t,y(:,i+1700),1000);
% end
% 
% for i = 1:100
%     x0(i+1800) = interp1(t,y(:,i+1800),1000);
% end
% 
% for i = 1:100
%     x0(i+1900) = interp1(t,y(:,i+1900),1000);
% end
% 
% 
% for i = 1:100
%     x0(i+2000) = F_1*final_AA(i);
% end
% 
% for i = 1:100
%     x0(i+2100) = F_2*final_AA(i);
% end
% 
% 
% for i = 1:400
%     x0(i+2200) = 0;
% end
% 
% 
% for i = 1:2
%     x0(i+2600) = T_1*prob_C(i);
% end
% 
% for i = 1:2
%     x0(i+2602) = T_2*prob_C(i);
% end
%  x0(2605) = E;
%  x0(2606) = 10;
% 
% 
% 
% 
% para(1) = 1e-4; % T and C binding constant
% para(2) = 1e6-0.6;  % 5e5-20.1
% para(3) = 1e-16; % 1e-8; % A and E binding constant 
% para(4) = 10/3.6-1;% 99
% para(5) = 0.4e13; % pi1 0.019*1e15
% para(6) = 0.8e8;% pi3 1.8e8 
% para(7) = 0.008; % decay constant for A1  0.038
% para(8) = 0.004;% decay constant for A2 0.019
% para(9) = 0.016; % decay constant for T1   0.036
% para(10) = 0.008;%% decay constant for T2  0.018
% para(11) = 0.5;% decay constant for C
% para(12) = 0.5;% decay constant for D
% 
% 
% para(13) = 0.6; % pi1'  0.1 
% para(14) = 0;% pi3'
% 
% para(15) = 0; % decay constant for E
% 
% 
% 
% 
% 
% para(16) = 7.2e13;% pi2  2e12
% para(17) = 0.8;% virus_replication constant
% para(18) = 1e5;% N
% para(19) = 2;% regenerate constant for B cell
% para(20) = 1;%% regeneration constant for T cell
% para(21) = 2;% amplication constant for virus antigen
% 
% for j = 1:10
%        para_new_1(10*(1-1)+j) = 10^(-21);
%        para_new_1(10*(2-1)+j) = 10^(-20);
%        para_new_1(10*(3-1)+j) = 10^(-19);
%        para_new_1(10*(4-1)+j) = 10^(-18);
%        para_new_1(10*(5-1)+j) = 10^(-17);
%        para_new_1(10*(6-1)+j) = 10^(-16);
%        para_new_1(10*(7-1)+j) = 10^(-15);
%        para_new_1(10*(8-1)+j) = 10^(-14);
%        para_new_1(10*(9-1)+j) = 10^(-13);
%        para_new_1(10*(10-1)+j) = 10^(-12);
% end
% 
% 
% 
%        
% for j = 1:10
%        para_new_2(10*(j-1)+1) = 1e-3;
%        para_new_2(10*(j-1)+2) = 1e-2;
%        para_new_2(10*(j-1)+3) = 1e-1;
%        para_new_2(10*(j-1)+4) = 1e0;
%        para_new_2(10*(j-1)+5) = 1e1;
%        para_new_2(10*(j-1)+6) = 1e2;
%        para_new_2(10*(j-1)+7) = 1e3;
%        para_new_2(10*(j-1)+8) = 1e4;
%        para_new_2(10*(j-1)+9) = 1e5;
%        para_new_2(10*(j-1)+10) = 1e6;
% end
% 
% para_new_3(1) = 1e-3;
% para_new_3(2) = 1e-8;
% 
% para_new_4(1) = 1e-3;
% para_new_4(2) = 1e6;
% 
% [t_new, y_new]=ode15s(@pathway_model_many_antibody_immune_include_T_complicated_2t_new,[1000 2000],x0,[],para,para_new_1,para_new_2,para_new_3,para_new_4,final_AA,prob_C);
% 
% plot(t_new,y_new(:,2606),'linewidth',2);
% hold on