function F=pathway_model_many_antibody_immune_imprinting(t,y,para,para_new,para_new_1,AA)

 y = max(0,y);
 F = zeros(1702,1);

 %% new addition of 300 terms
 % orignial 201-300 BCR IgG 
 % 601-700
 %

 % 1402-1502  100 BCR IgG from IgM
 % 1502-1602  100 BCR IgG-E Complex
 % 1602-1702  100 BCR-IgG-V Complex




for i = 1:10
    for j = 1:10
        F(10*(i-1)+j,1) = para(4)*AA(10*(i-1)+j)-para(1)*y(10*(i-1)+j)*y(1401)+para(2)*y(10*(i-1)+j+400)+para(7)*y(10*(i-1)+j+400)...
            -5*para_new(i)*y(10*(i-1)+j)*y(1402)+para_new_1(j)*y(10*(i-1)+j+800)...
            +1.05*para(15)*(1-para(16))*(1-para(17)*y(10*(i-1)+j+800)/(y(10*(i-1)+j+800)+y(10*(i-1)+j)))*y(10*(i-1)+j+800)-para(5)*y(10*(i-1)+j);
      
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+100,1) = -para(1)*y(10*(i-1)+j+100)*y(1401)+para(2)*y(10*(i-1)+j+500)+para(11)*y(10*(i-1)+j+1200)...
            -5*para_new(i)*y(10*(i-1)+j+100)*y(1402)+para_new_1(j)*y(10*(i-1)+j+900)...
            -para(13)*y(10*(i-1)+j+100);
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+200,1) = -para(1)*y(10*(i-1)+j+200)*y(1401)+para(2)*y(10*(i-1)+j+600)+para(7)*y(10*(i-1)+j+600)...
            -para_new(i)*y(10*(i-1)+j+200)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1000)...
            +1.0*para(15)*(1-para(17)*(y(10*(i-1)+j+1000)+y(10*(i-1)+j+1602))/(y(10*(i-1)+j+1000)+y(10*(i-1)+j+200)+y(10*(i-1)+j+1602)+y(10*(i-1)+j+1402)))*y(10*(i-1)+j+1000)-para(6)*y(10*(i-1)+j+200);
        
    end
end

%% NEW ADDED TERM

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1402,1) = -para(1)*y(10*(i-1)+j+1402)*y(1401)+para(2)*y(10*(i-1)+j+1502)+para(7)*y(10*(i-1)+j+1502)...
            -para_new(i)*y(10*(i-1)+j+1402)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1602)...
            +1.0*para(15)*para(16)*(1-para(17)*y(10*(i-1)+j+800)/(y(10*(i-1)+j+800)+y(10*(i-1)+j)))*y(10*(i-1)+j+800)...
            +1.0*para(15)*(1-para(17)*(y(10*(i-1)+j+1000)+y(10*(i-1)+j+1602))/(y(10*(i-1)+j+1000)+y(10*(i-1)+j+200)+y(10*(i-1)+j+1602)+y(10*(i-1)+j+1402)))*y(10*(i-1)+j+1602)-para(6)*y(10*(i-1)+j+1402);
        
        
    end
end
%% 
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+300,1) =  -para(1)*y(10*(i-1)+j+300)*y(1401)+para(2)*y(10*(i-1)+j+700)+para(12)*y(10*(i-1)+j+1300)...
            -para_new(i)*y(10*(i-1)+j+300)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1100)...
            -para(14)*y(10*(i-1)+j+300);
        

    end
end


%% c1-c4

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+400,1) =  para(1)*y(10*(i-1)+j)*y(1401)-para(2)*y(10*(i-1)+j+400)-para(18)*y(10*(i-1)+j+400);
        F_E_1(i,j) = -para(1)*y(10*(i-1)+j)*y(1401)+ para(2)*y(10*(i-1)+j+400);

    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+500,1) = para(1)*y(10*(i-1)+j+100)*y(1401)-para(2)*y(10*(i-1)+j+500)-para(18)*y(10*(i-1)+j+500);
        F_E_2(i,j) = -para(1)*y(10*(i-1)+j+100)*y(1401)+para(2)*y(10*(i-1)+j+500);
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+600,1) =  para(1)*y(10*(i-1)+j+200)*y(1401)-para(2)*y(10*(i-1)+j+600)-para(18)*y(10*(i-1)+j+600);
        F_E_3(i,j) = -para(1)*y(10*(i-1)+j+200)*y(1401)+para(2)*y(10*(i-1)+j+600);

    end
end
%%
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1502,1) =  para(1)*y(10*(i-1)+j+1402)*y(1401)-para(2)*y(10*(i-1)+j+1502)-para(18)*y(10*(i-1)+j+1502);
        F_E_5(i,j) = -para(1)*y(10*(i-1)+j+1402)*y(1401)+para(2)*y(10*(i-1)+j+1502);

    end
end

%%
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+700,1) = para(1)*y(10*(i-1)+j+300)*y(1401)-para(2)*y(10*(i-1)+j+700)-para(18)*y(10*(i-1)+j+700);
        F_E_4(i,j) = -para(1)*y(10*(i-1)+j+300)*y(1401)+para(2)*y(10*(i-1)+j+700);
    end
end
%% c1'-c4'

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+800,1) = 5*para_new(i)*y(10*(i-1)+j)*y(1402)-para_new_1(j)*y(10*(i-1)+j+800)-para(18)*y(10*(i-1)+j+800);
        F_V_1(i,j) = -5*para_new(i)*y(10*(i-1)+j)*y(1402)+para_new_1(j)*y(10*(i-1)+j+800);    
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+900,1) = 5*para_new(i)*y(10*(i-1)+j+100)*y(1402)-para_new_1(j)*y(10*(i-1)+j+900)-para(18)*y(10*(i-1)+j+900);
        F_V_2(i,j) = -5*para_new(i)*y(10*(i-1)+j+100)*y(1402)+para_new_1(j)*y(10*(i-1)+j+900);         
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1000,1) = para_new(i)*y(10*(i-1)+j+200)*y(1402)-para_new_1(j)*y(10*(i-1)+j+1000)-para(18)*y(10*(i-1)+j+1000);
        F_V_3(i,j) = -para_new(i)*y(10*(i-1)+j+200)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1000);      
    end
end
%%
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1602,1) = para_new(i)*y(10*(i-1)+j+1402)*y(1402)-para_new_1(j)*y(10*(i-1)+j+1602)-para(18)*y(10*(i-1)+j+1602);
        F_V_5(i,j) = -para_new(i)*y(10*(i-1)+j+1402)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1602);      
    end
end

%%

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1100,1) = para_new(i)*y(10*(i-1)+j+300)*y(1402)-para_new_1(j)*y(10*(i-1)+j+1100)-para(18)*y(10*(i-1)+j+1100);
        F_V_4(i,j) = -para_new(i)*y(10*(i-1)+j+300)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1100);      
    end
end

%% plasma_M
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1200,1) = para(8)*y(10*(i-1)+j+400)+1.0*para(15)*(1-para(16))*para(17)*y(10*(i-1)+j+800)/(y(10*(i-1)+j+800)+y(10*(i-1)+j))*y(10*(i-1)+j+800)/para(20)-para(9)*y(10*(i-1)+j+1200);
             
    end
end

%% plasma_G
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1300,1) = para(8)*y(10*(i-1)+j+600)+1.0*para(15)*para(17)*(y(10*(i-1)+j+1000)+y(10*(i-1)+j+1602))/(y(10*(i-1)+j+1000)+y(10*(i-1)+j+200)+y(10*(i-1)+j+1602)+y(10*(i-1)+j+1402))*(y(10*(i-1)+j+1000)+y(10*(i-1)+j+1602))/para(20)-para(10)*y(10*(i-1)+j+1300);
             
    end
end





F(1401,1) = para(3)+sum(sum(F_E_1))+sum(sum(F_E_2))+sum(sum(F_E_3))+sum(sum(F_E_4))+sum(sum(F_E_5));
F(1402,1) = para(19)*y(1402)+sum(sum(F_V_1))+sum(sum(F_V_2))+sum(sum(F_V_3))+sum(sum(F_V_4))+sum(sum(F_V_5));




end

