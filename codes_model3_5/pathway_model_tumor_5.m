function F=pathway_model_tumor_5(t,y,para)

F(1,1) = para(17) + para(1)*y(1)-para(2)*y(1)-para(3)*y(4)/(y(4)+para(4))*y(1)-para(5)*y(5)*y(1);

F(2,1) = para(8)*(para(2)*y(1)+para(3)*y(4)/(y(4)+para(4))*y(1)+para(5)*y(5)*y(1))-para(6)*y(2)*y(4)+para(7)*y(3);

F(3,1) = para(6)*y(2)*y(4)-para(7)*y(3)-para(9)*y(3);

F(4,1) = -para(6)*y(2)*y(4)+para(7)*y(3)+para(10)*y(3)-para(11)*y(4)+para(12);

F(5,1) = -para(13)*para(5)*y(1)*y(5)+para(14)*y(3)+para(15)-para(16)*y(5);


end





