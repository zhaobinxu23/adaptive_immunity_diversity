function F=pathway_model_many_antibody_immune_include_T_complicated_2b_new(t,y,para,para_new_1,para_new_2,para_new_3,para_new_4,AA,BB)

%  y = max(0,y);
F = zeros(1822,1);
F_new_t_1 = zeros(2,100);
F_new_c_1_g = zeros(2,100);
F_new_t_1_g = zeros(2,100);

F_new_t_2 = zeros(2,100);
F_new_c_2_g = zeros(2,100);
F_new_t_2_g = zeros(2,100);

F_new_t_3 = zeros(2,100);
F_new_c_3_g = zeros(2,100);
F_new_t_3_g = zeros(2,100);

F_new_t_4 = zeros(2,100);
F_new_c_4_g = zeros(2,100);
F_new_t_4_g = zeros(2,100);

F_new_t_5 = zeros(2,100);
F_new_c_5_g = zeros(2,100);
F_new_t_5_g = zeros(2,100);
F_new_c_5_g_2 = zeros(2,100);
F_new_t_5_g_2 = zeros(2,100);


F_new_t_6 = zeros(2,100);
F_new_c_6_g = zeros(2,100);
F_new_t_6_g = zeros(2,100);
F_new_c_6_g_2 = zeros(2,100);

F_new_t_7 = zeros(2,100);
F_new_c_7_g = zeros(2,100);
F_new_t_7_g = zeros(2,100);
F_new_t_7_g_2 = zeros(2,100);

F_new_t_8 = zeros(2,100);
F_new_c_8_g = zeros(2,100);
F_new_t_8_g = zeros(2,100);


for i = 1:2
    for j= 1:100
        F(100*(i-1)+j,1) = para(1)*y(1608+i)/para(18)*y(1620+j)-para(2)*y(100*(i-1)+j)-para(13)*y(100*(i-1)+j);
        F_new_c_1(i,j) = (-para(1)*y(1608+i)/para(18)*y(1620+j)+para(2)*y(100*(i-1)+j))*para(18);
        F_new_t_1(i,j) = -para(1)*y(1608+i)/para(18)*y(1620+j)+para(2)*y(100*(i-1)+j);
        F_new_c_1_g(i,j) = para(18)*para(19)*y(100*(i-1)+j);
        % F_new_c_1_g_2(i,j) = para(18)*para(19)*y(100*(i-1)+j)*0.05;
        F_new_t_1_g(i,j) = para(20)*y(100*(i-1)+j);
    end
end

for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+200,1) = para(1)*y(1608+i)/para(18)*y(1720+j)-para(2)*y(100*(i-1)+j+200)-para(13)*y(100*(i-1)+j+200);
        F_new_c_2(i,j) = (-para(1)*y(1608+i)/para(18)*y(1720+j)+para(2)*y(100*(i-1)+j+200))*para(18);
        F_new_t_2(i,j) = -para(1)*y(1608+i)/para(18)*y(1720+j)+para(2)*y(100*(i-1)+j+200);
        F_new_c_2_g(i,j) = para(18)*para(19)*y(100*(i-1)+j+200);
        % F_new_c_2_g_2(i,j) = para(18)*para(19)*y(100*(i-1)+j+200)*0.05;
        F_new_t_2_g(i,j) = para(20)*y(100*(i-1)+j+200);
    end
end


for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+400,1) = para(1)*y(1610+i)/para(18)*y(1620+j)-para(2)*y(100*(i-1)+j+400)-para(13)*y(100*(i-1)+j+400);
        F_new_c_3(i,j) = (-para(1)*y(1610+i)/para(18)*y(1620+j)+para(2)*y(100*(i-1)+j+400))*para(18);
        F_new_t_3(i,j) = -para(1)*y(1610+i)/para(18)*y(1620+j)+para(2)*y(100*(i-1)+j+400);
        F_new_c_3_g(i,j) = para(18)*para(19)*y(100*(i-1)+j+400);
        F_new_t_3_g(i,j) = para(20)*y(100*(i-1)+j+400);
    end
end

for i = 1:2
    for j = 1:100
        F(100*(i-1)+j+600,1) = para(1)*y(1610+i)/para(18)*y(1720+j)-para(2)*y(100*(i-1)+j+600)-para(13)*y(100*(i-1)+j+600);
        F_new_c_4(i,j) = (-para(1)*y(1610+i)/para(18)*y(1720+j)+para(2)*y(100*(i-1)+j++600))*para(18);
        F_new_t_4(i,j) = -para(1)*y(1610+i)/para(18)*y(1720+j)+para(2)*y(100*(i-1)+j+600);
        F_new_c_4_g(i,j) = para(18)*para(19)*y(100*(i-1)+j+600);
        F_new_t_4_g(i,j) = para(20)*y(100*(i-1)+j+600);
    end
end




for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+800,1) = para_new_3(j)*y(1616+i)/para(18)*y(1620+j)-para_new_4(j)*y(100*(i-1)+j+800)-para(13)*y(100*(i-1)+j+800);
        F_new_c_5(i,j) = (-para_new_3(j)*y(1616+i)/para(18)*y(1620+j)+para_new_4(j)*y(100*(i-1)+j+800))*para(18);
        F_new_t_5(i,j) = -para_new_3(j)*y(1616+i)/para(18)*y(1620+j)+para_new_4(j)*y(100*(i-1)+j+800);
        F_new_c_5_g(i,j) = para(21)*para(18)*para(19)*y(100*(i-1)+j+800)*0.95;
        F_new_t_5_g(i,j) = para(21)*para(20)*y(100*(i-1)+j+800)*0.95;
        F_new_c_5_g_2(i,j) = para(21)*para(18)*para(19)*y(100*(i-1)+j+800)*0.05;
        F_new_t_5_g_2(i,j) = para(21)*para(20)*y(100*(i-1)+j+800)*0.05;
    end
end

for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+1000,1) = para_new_3(j)*y(1616+i)/para(18)*y(1720+j)-para_new_4(j)*y(100*(i-1)+j+1000)-para(13)*y(100*(i-1)+j+1000);
        F_new_c_6(i,j) = (-para_new_3(j)*y(1616+i)/para(18)*y(1720+j)+para_new_4(j)*y(100*(i-1)+j+1000))*para(18);
        F_new_t_6(i,j) = -para_new_3(j)*y(1616+i)/para(18)*y(1720+j)+para_new_4(j)*y(100*(i-1)+j+1000);
        
        F_new_c_6_g(i,j) = para(21)*para(18)*para(19)*y(100*(i-1)+j+1000)*0.95;
        F_new_t_6_g(i,j) = para(21)*para(20)*y(100*(i-1)+j+1000);
        F_new_c_6_g_2(i,j) = para(21)*para(18)*para(19)*y(100*(i-1)+j+1000)*0.05;
        
    end
end


for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+1200,1) = para_new_3(j)*y(1618+i)/para(18)*y(1620+j)-para_new_4(j)*y(100*(i-1)+j+1200)-para(13)*y(100*(i-1)+j+1200);
        F_new_c_7(i,j) = (-para_new_3(j)*y(1618+i)/para(18)*y(1620+j)+para_new_4(j)*y(100*(i-1)+j+1200))*para(18);
        F_new_t_7(i,j) = -para_new_3(j)*y(1618+i)/para(18)*y(1620+j)+para_new_4(j)*y(100*(i-1)+j+1200);
        
        F_new_c_7_g(i,j) = para(21)*0.5*para(18)*para(19)*y(100*(i-1)+j+1200);
        F_new_t_7_g(i,j) = para(21)*0.5*para(20)*y(100*(i-1)+j+1200)*0.95;      
        F_new_t_7_g_2(i,j) = para(21)*0.5*para(20)*y(100*(i-1)+j+1200)*0.05;
    end
end

for i = 1:2
    for j= 1:100
        F(100*(i-1)+j+1400,1) = para_new_3(j)*y(1618+i)/para(18)*y(1720+j)-para_new_4(j)*y(100*(i-1)+j+1400)-para(13)*y(100*(i-1)+j+1400);
        F_new_c_8(i,j) = (-para_new_3(j)*y(1618+i)/para(18)*y(1720+j)+para_new_4(j)*y(100*(i-1)+j+1400))*para(18);
        F_new_t_8(i,j) = -para_new_3(j)*y(1618+i)/para(18)*y(1720+j)+para_new_4(j)*y(100*(i-1)+j+1400);
        
        F_new_c_8_g(i,j) = para(21)*0.5*para(18)*para(19)*y(100*(i-1)+j+1400);
        F_new_t_8_g(i,j) = para(21)*0.5*para(20)*y(100*(i-1)+j+1400);
        
    end
end




for i = 1:2
    F(i+1600,1) = para(5)*AA(i)-para(3)*y(1600+i)*y(1821)+para(4)*y(1604+i)+sum(F_new_c_1_g(i,:)) + sum(F_new_c_2_g(i,:))-10*para_new_1(i)*y(1600+i)*y(1822)+para_new_2(i)*y(1612+i) + sum(F_new_c_5_g(i,:)) + sum(F_new_c_6_g(i,:))-para(7)*y(1600+i);
    F_new_E(i) = -para(3)*y(1600+i)*y(1821)+para(4)*y(1604+i);
    F_new_V(i) = -10*para_new_1(i)*y(1600+i)*y(1822)+para_new_2(i)*y(1612+i);
end


for i = 1:2
    F(i+1602,1) = -para(3)*y(1602+i)*y(1821)+para(4)*y(1606+i)+sum(F_new_c_3_g(i,:)) + sum(F_new_c_4_g(i,:))-para_new_1(i)*y(1602+i)*y(1822)+para_new_2(i)*y(1614+i) + sum(F_new_c_7_g(i,:)) + sum(F_new_c_8_g(i,:)) +sum(F_new_c_5_g_2(i,:)) + sum(F_new_c_6_g_2(i,:))-para(8)*y(1602+i);
    F_new_E_2(i) = -para(3)*y(1602+i)*y(1821)+para(4)*y(1606+i);
    F_new_V_2(i) = -para_new_1(i)*y(1602+i)*y(1822)+para_new_2(i)*y(1614+i);
end


for i = 1:2
    
    F(i+1604,1) = para(3)*y(1600+i)*y(1821)-para(4)*y(1604+i)-para(11)*y(i+1604)-para(12)*y(i+1604);
    
end


for i = 1:2
    
    F(i+1606,1) = para(3)*y(1602+i)*y(1821)-para(4)*y(1606+i)-para(11)*y(1606+i)-para(12)*y(i+1606);
    
end



for i = 1:2
    
    F(i+1608,1) = para(12)*y(i+1604)-para(13)*y(i+1608)+sum(F_new_c_1(i,:))+sum(F_new_c_2(i,:));
    
end

for i = 1:2
    
    F(i+1610,1) = para(12)*y(i+1606)-para(13)*y(i+1610)+sum(F_new_c_3(i,:))+sum(F_new_c_4(i,:));
    
end




for i = 1:2
    
    F(i+1612,1) = 10*para_new_1(i)*y(1600+i)*y(1822)-para_new_2(i)*y(1612+i)-para(11)*y(1612+i)-para(12)*y(1612+i);
    
end


for i = 1:2
    
    F(i+1614,1) = para_new_1(i)*y(1602+i)*y(1822)-para_new_2(i)*y(1614+i)-para(11)*y(1614+i)-para(12)*y(1614+i);
    
end


for i = 1:2
    
    F(i+1616,1) = para(12)*y(i+1612)-para(13)*y(i+1616)+sum(F_new_c_5(i,:))+sum(F_new_c_6(i,:));
    
end

for i = 1:2
    
    F(i+1618,1) = para(12)*y(i+1614)-para(13)*y(i+1618)+sum(F_new_c_7(i,:))+sum(F_new_c_8(i,:));
    
end





for i = 1:100
    
    F(i+1620,1) = para(6)*BB(i)+sum(F_new_t_1(:,i))+sum(F_new_t_3(:,i))+sum(F_new_t_5(:,i))+sum(F_new_t_7(:,i))-para(9)*y(1620+i)+sum(F_new_t_1_g(:,i))+sum(F_new_t_3_g(:,i))+sum(F_new_t_5_g(:,i))+sum(F_new_t_7_g(:,i));
    
end


for i = 1:100
    
    F(i+1720,1) = sum(F_new_t_2(:,i))+sum(F_new_t_4(:,i))+sum(F_new_t_6(:,i))+sum(F_new_t_8(:,i))-para(10)*y(1720+i)+sum(F_new_t_2_g(:,i))+sum(F_new_t_4_g(:,i))+sum(F_new_t_6_g(:,i))+sum(F_new_t_8_g(:,i))+ sum(F_new_t_5_g_2(:,i))+sum(F_new_t_7_g_2(:,i));
    
end


F(1821,1) = para(16) + sum(F_new_E)+sum(F_new_E_2);
F(1822,1) = para(17)*y(1822)+sum(F_new_V)+sum(F_new_V_2);




end

