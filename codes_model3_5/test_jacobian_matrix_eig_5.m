clc
clear
syms y1 y2 y3 y4 y5 


para(1) = 1.0; % k1
para(2) = 0.05; % d1
para(3) = 10;% β1
para(4) = 1e7;% km
para(5) = 1e-4;% α1
para(6) = 1e-7;% k2
para(7) = 0;% k-2
para(8) = 1e4;% ρ1
para(9) = 1;% k3
para(10) = 1.8;% k4
para(11) = 1e-3;% d2
para(12) = 0;% Π2
para(13) = 0.1;% ρ2
para(14) = 2e-4;% k5
para(15) = 1e0;% Π3
para(16) = 1e-3;% d3
para(17) = 100;

para(19) = 1e2;%% 1e0

eq1 = para(17) + para(1)*y1-para(2)*y1-para(3)*y4/(y4+para(4))*y1-para(5)*y5*y1 == 0;

eq2 = para(8)*(para(2)*y1+para(3)*y4/(y4+para(4))*y1+para(5)*y5*y1)-para(6)*y2*y4+para(7)*y3 == 0;

eq3 = para(6)*y2*y4-para(7)*y3-para(9)*y3 == 0;

eq4 = -para(6)*y2*y4+para(7)*y3+para(10)*y3-para(11)*y4+para(19) == 0;

eq5 = -para(13)*para(5)*y1*y5+para(14)*y3+para(15)-para(16)*y5 == 0;




sol = solve([eq1, eq2, eq3, eq4, eq5], [y1, y2, y3, y4, y5]);
% disp([sol.y1, sol.y2, sol.y3,sol.y4,sol.y5,sol.y6]);
x_num(:,1) = vpa(sol.y1, 5);
x_num(:,2) = vpa(sol.y2, 5);
x_num(:,3) = vpa(sol.y3, 5);
x_num(:,4) = vpa(sol.y4, 5);
x_num(:,5) = vpa(sol.y5, 5);



x0(1) = x_num(3,1); %% Tu  466.70152158166456501930952072144
x0(2) = x_num(3,2); %% V  9999.9948277755174785852432250977
x0(3) = x_num(3,3); %% C  1933403.043163299560546875
x0(4) = x_num(3,4); %% A  1933404043.163330078125
x0(5) = x_num(3,5);%% Tc  264321.40618126433498247251918656




xxx =  zeros(5,5);

xxx(1,1) = para(1)-para(2)-para(3)*x0(4)/(x0(4)+para(4));
xxx(1,4) = -para(4)*para(3)*x0(1)/(x0(4)+para(4))^2;
xxx(1,5) = -para(5)*x0(1);

xxx(2,1) = para(8)*(para(3)*x0(4)/(x0(4)+para(4))+para(5)*x0(5)+para(2));
xxx(2,2) = -para(6)*x0(4);
xxx(2,3) = para(7);
xxx(2,4) = -para(6)*x0(2)+para(8)*para(4)*para(3)*x0(1)/(x0(4)+para(4))^2;
xxx(2,5) = para(8)*para(5)*x0(1);

xxx(3,2) = para(6)*x0(4);
xxx(3,3) = -para(7)-para(9);
xxx(3,4) = para(6)*x0(2);

xxx(4,2) = -para(6)*x0(2);
xxx(4,3) = para(7)+para(10);
xxx(4,4) = -para(6)*x0(2)-para(11);


xxx(5,1) = -para(13)*para(5)*x0(5);
xxx(5,3) = para(14);
xxx(5,5) = -para(13)*para(5)*x0(1)-para(16);



eig(xxx)


