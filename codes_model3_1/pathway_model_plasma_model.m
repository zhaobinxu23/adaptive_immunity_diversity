function F=pathway_model_plasma_model(t,y,para)

 y = max(0,y);
F = zeros(6,1);

F(1,1) = para(12)*y(6)-para(1)*y(1)*y(5)+para(2)*y(3)-para(6)*y(1);

F(2,1) = para(3)*(1-para(13)*y(4)/(y(4)+y(2)))*y(4)-para(10)*y(2)*y(5)+para(2)*y(4)-para(5)*y(2);

F(3,1) = para(1)*y(1)*y(5)-para(2)*y(3)-para(7)*y(3);

F(4,1) = para(10)*y(2)*y(5)-para(2)*y(4)-para(8)*y(4);

F(5,1) = para(9)*y(5)-para(1)*y(1)*y(5)-para(10)*y(2)*y(5)+para(2)*(y(3)+y(4));

F(6,1) = para(4)*para(3)*para(13)*y(4)/(y(4)+y(2))*y(4)-para(11)*y(6);



end

