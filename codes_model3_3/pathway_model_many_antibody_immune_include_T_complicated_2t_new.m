function F=pathway_model_many_antibody_immune_include_T_complicated_2t_new(t,y,para,para_new_1,para_new_2,para_new_3,para_new_4,AA,BB)

%  y = max(0,y);
F = zeros(2606,1);
F_new_t_1 = zeros(100,2);
F_new_c_1_g = zeros(100,2);
F_new_t_1_g = zeros(100,2);

F_new_t_2 = zeros(100,2);
F_new_c_2_g = zeros(100,2);
F_new_t_2_g = zeros(100,2);

F_new_t_3 = zeros(100,2);
F_new_c_3_g = zeros(100,2);
F_new_t_3_g = zeros(100,2);

F_new_t_4 = zeros(100,2);
F_new_c_4_g = zeros(100,2);
F_new_t_4_g = zeros(100,2);

F_new_t_5 = zeros(100,2);
F_new_c_5_g = zeros(100,2);
F_new_t_5_g = zeros(100,2);
F_new_c_5_g_2 = zeros(100,2);
F_new_t_5_g_2 = zeros(100,2);


F_new_t_6 = zeros(100,2);
F_new_c_6_g = zeros(100,2);
F_new_t_6_g = zeros(100,2);
F_new_c_6_g_2 = zeros(100,2);

F_new_t_7 = zeros(100,2);
F_new_c_7_g = zeros(100,2);
F_new_t_7_g = zeros(100,2);
F_new_t_7_g_2 = zeros(100,2);

F_new_t_8 = zeros(100,2);
F_new_c_8_g = zeros(100,2);
F_new_t_8_g = zeros(100,2);

for i = 1:100
    for j= 1:2
        F(2*(i-1)+j,1) = para(1)*y(2000+i)/para(18)*y(2600+j)-para(2)*y(2*(i-1)+j)-para(13)*y(2*(i-1)+j);
        F_new_c_1(i,j) = (-para(1)*y(2000+i)/para(18)*y(2600+j)+para(2)*y(2*(i-1)+j))*para(18);
        F_new_t_1(i,j) = -para(1)*y(2000+i)/para(18)*y(2600+j)+para(2)*y(2*(i-1)+j);
        F_new_c_1_g(i,j) = para(18)*para(19)*y(2*(i-1)+j);
        % F_new_c_1_g_2(i,j) = para(18)*para(19)*y(2*(i-1)+j)*0.05;
        F_new_t_1_g(i,j) = para(20)*y(2*(i-1)+j);
    end
end

for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+200,1) = para(1)*y(2000+i)/para(18)*y(2602+j)-para(2)*y(2*(i-1)+j+200)-para(13)*y(2*(i-1)+j+200);
        F_new_c_2(i,j) = (-para(1)*y(2000+i)/para(18)*y(2602+j)+para(2)*y(2*(i-1)+j+200))*para(18);
        F_new_t_2(i,j) = -para(1)*y(2000+i)/para(18)*y(2602+j)+para(2)*y(2*(i-1)+j+200);
        F_new_c_2_g(i,j) = para(18)*para(19)*y(2*(i-1)+j+200);
        % F_new_c_2_g_2(i,j) = para(18)*para(19)*y(2*(i-1)+j+200)*0.05;
        F_new_t_2_g(i,j) = para(20)*y(2*(i-1)+j+200);
    end
end


for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+400,1) = para(1)*y(2100+i)/para(18)*y(2600+j)-para(2)*y(2*(i-1)+j+400)-para(13)*y(2*(i-1)+j+400);
        F_new_c_3(i,j) = (-para(1)*y(2100+i)/para(18)*y(2600+j)+para(2)*y(2*(i-1)+j+400))*para(18);
        F_new_t_3(i,j) = -para(1)*y(2100+i)/para(18)*y(2600+j)+para(2)*y(2*(i-1)+j+400);
        F_new_c_3_g(i,j) = para(18)*para(19)*y(2*(i-1)+j+400);
        F_new_t_3_g(i,j) = para(20)*y(2*(i-1)+j+400);
    end
end

for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+600,1) = para(1)*y(2100+i)/para(18)*y(2602+j)-para(2)*y(2*(i-1)+j+600)-para(13)*y(2*(i-1)+j+600);
        F_new_c_4(i,j) = (-para(1)*y(2100+i)/para(18)*y(2602+j)+para(2)*y(2*(i-1)+j++600))*para(18);
        F_new_t_4(i,j) = -para(1)*y(2100+i)/para(18)*y(2602+j)+para(2)*y(2*(i-1)+j+600);
        F_new_c_4_g(i,j) = para(18)*para(19)*y(2*(i-1)+j+600);
        F_new_t_4_g(i,j) = para(20)*y(2*(i-1)+j+600);
    end
end




for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+800,1) = para_new_3(j)*y(2400+i)/para(18)*y(2600+j)-para_new_4(j)*y(2*(i-1)+j+800)-para(13)*y(2*(i-1)+j+800);
        F_new_c_5(i,j) = (-para_new_3(j)*y(2400+i)/para(18)*y(2600+j)+para_new_4(j)*y(2*(i-1)+j+800))*para(18);
        F_new_t_5(i,j) = -para_new_3(j)*y(2400+i)/para(18)*y(2600+j)+para_new_4(j)*y(2*(i-1)+j+800);
        F_new_c_5_g(i,j) = para(21)*para(18)*para(19)*y(2*(i-1)+j+800)*0.95;
        F_new_t_5_g(i,j) = para(21)*para(20)*y(2*(i-1)+j+800)*0.95;
        F_new_c_5_g_2(i,j) = para(21)*para(18)*para(19)*y(2*(i-1)+j+800)*0.05;
        F_new_t_5_g_2(i,j) = para(21)*para(20)*y(2*(i-1)+j+800)*0.05;
    end
end

for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+1000,1) = para_new_3(j)*y(2400+i)/para(18)*y(2602+j)-para_new_4(j)*y(2*(i-1)+j+1000)-para(13)*y(2*(i-1)+j+1000);
        F_new_c_6(i,j) = (-para_new_3(j)*y(2400+i)/para(18)*y(2602+j)+para_new_4(j)*y(2*(i-1)+j+1000))*para(18);
        F_new_t_6(i,j) = -para_new_3(j)*y(2400+i)/para(18)*y(2602+j)+para_new_4(j)*y(2*(i-1)+j+1000);
        
        F_new_c_6_g(i,j) = para(21)*para(18)*para(19)*y(2*(i-1)+j+1000)*0.95;
        F_new_t_6_g(i,j) = para(21)*para(20)*y(2*(i-1)+j+1000);
        F_new_c_6_g_2(i,j) = para(21)*para(18)*para(19)*y(2*(i-1)+j+1000)*0.05;
        
    end
end


for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+1200,1) = para_new_3(j)*y(2500+i)/para(18)*y(2600+j)-para_new_4(j)*y(2*(i-1)+j+1200)-para(13)*y(2*(i-1)+j+1200);
        F_new_c_7(i,j) = (-para_new_3(j)*y(2500+i)/para(18)*y(2600+j)+para_new_4(j)*y(2*(i-1)+j+1200))*para(18);
        F_new_t_7(i,j) = -para_new_3(j)*y(2500+i)/para(18)*y(2600+j)+para_new_4(j)*y(2*(i-1)+j+1200);
        
        F_new_c_7_g(i,j) = para(21)*0.5*para(18)*para(19)*y(2*(i-1)+j+1200);
        F_new_t_7_g(i,j) = para(21)*0.5*para(20)*y(2*(i-1)+j+1200)*0.95;      
        F_new_t_7_g_2(i,j) = para(21)*0.5*para(20)*y(2*(i-1)+j+1200)*0.05;
    end
end

for i = 1:100
    for j= 1:2
        F(2*(i-1)+j+1400,1) = para_new_3(j)*y(2500+i)/para(18)*y(2602+j)-para_new_4(j)*y(2*(i-1)+j+1400)-para(13)*y(2*(i-1)+j+1400);
        F_new_c_8(i,j) = (-para_new_3(j)*y(2500+i)/para(18)*y(2602+j)+para_new_4(j)*y(2*(i-1)+j+1400))*para(18);
        F_new_t_8(i,j) = -para_new_3(j)*y(2500+i)/para(18)*y(2602+j)+para_new_4(j)*y(2*(i-1)+j+1400);
        
        F_new_c_8_g(i,j) = para(21)*0.5*para(18)*para(19)*y(2*(i-1)+j+1400);
        F_new_t_8_g(i,j) = para(21)*0.5*para(20)*y(2*(i-1)+j+1400);
        
    end
end




for i = 1:100
    F(i+1600,1) = para(5)*AA(i)-para(3)*y(1600+i)*y(2605)+para(4)*y(1800+i)+sum(F_new_c_1_g(i,:)) + sum(F_new_c_2_g(i,:))-10*para_new_1(i)*y(1600+i)*y(2606)+para_new_2(i)*y(2200+i) + sum(F_new_c_5_g(i,:)) + sum(F_new_c_6_g(i,:))-para(7)*y(1600+i);
    F_new_E(i) = -para(3)*y(1600+i)*y(2605)+para(4)*y(1800+i);
    F_new_V(i) = -10*para_new_1(i)*y(1600+i)*y(2606)+para_new_2(i)*y(2200+i);
end


for i = 1:100
    F(i+1700,1) = -para(3)*y(1700+i)*y(2605)+para(4)*y(1900+i)+sum(F_new_c_3_g(i,:)) + sum(F_new_c_4_g(i,:))-para_new_1(i)*y(1700+i)*y(2606)+para_new_2(i)*y(2300+i) + sum(F_new_c_7_g(i,:)) + sum(F_new_c_8_g(i,:)) +sum(F_new_c_5_g_2(i,:)) + sum(F_new_c_6_g_2(i,:))-para(8)*y(1700+i);
    F_new_E_2(i) = -para(3)*y(1700+i)*y(2605)+para(4)*y(1900+i);
    F_new_V_2(i) = -para_new_1(i)*y(1700+i)*y(2606)+para_new_2(i)*y(2300+i);
end


for i = 1:100
    
    F(i+1800,1) = para(3)*y(1600+i)*y(2605)-para(4)*y(1800+i)-para(11)*y(i+1800)-para(12)*y(i+1800);
    
end


for i = 1:100
    
    F(i+1900,1) = para(3)*y(1700+i)*y(2605)-para(4)*y(1900+i)-para(11)*y(1900+i)-para(12)*y(i+1900);
    
end



for i = 1:100
    
    F(i+2000,1) = para(12)*y(i+1800)-para(13)*y(i+2000)+sum(F_new_c_1(i,:))+sum(F_new_c_2(i,:));
    
end

for i = 1:100
    
    F(i+2100,1) = para(12)*y(i+1900)-para(13)*y(i+2100)+sum(F_new_c_3(i,:))+sum(F_new_c_4(i,:));
    
end




for i = 1:100
    
    F(i+2200,1) = 10*para_new_1(i)*y(1600+i)*y(2606)-para_new_2(i)*y(2200+i)-para(11)*y(2200+i)-para(12)*y(2200+i);
    
end


for i = 1:100
    
    F(i+2300,1) = para_new_1(i)*y(1700+i)*y(2606)-para_new_2(i)*y(2300+i)-para(11)*y(2300+i)-para(12)*y(2300+i);
    
end


for i = 1:100
    
    F(i+2400,1) = para(12)*y(i+2200)-para(13)*y(i+2400)+sum(F_new_c_5(i,:))+sum(F_new_c_6(i,:));
    
end

for i = 1:100
    
    F(i+2500,1) = para(12)*y(i+2300)-para(13)*y(i+2500)+sum(F_new_c_7(i,:))+sum(F_new_c_8(i,:));
    
end





for i = 1:2
    
    F(i+2600,1) = para(6)*BB(i)+sum(F_new_t_1(:,i))+sum(F_new_t_3(:,i))+sum(F_new_t_5(:,i))+sum(F_new_t_7(:,i))-para(9)*y(2600+i)+sum(F_new_t_1_g(:,i))+sum(F_new_t_3_g(:,i))+sum(F_new_t_5_g(:,i))+sum(F_new_t_7_g(:,i));
    
end


for i = 1:2
    
    F(i+2602,1) = sum(F_new_t_2(:,i))+sum(F_new_t_4(:,i))+sum(F_new_t_6(:,i))+sum(F_new_t_8(:,i))-para(10)*y(2602+i)+sum(F_new_t_2_g(:,i))+sum(F_new_t_4_g(:,i))+sum(F_new_t_6_g(:,i))+sum(F_new_t_8_g(:,i))+ sum(F_new_t_5_g_2(:,i))+sum(F_new_t_7_g_2(:,i));
    
end


F(2605,1) = para(16) + sum(F_new_E)+sum(F_new_E_2);
F(2606,1) = para(17)*y(2606)+sum(F_new_V)+sum(F_new_V_2);




end

