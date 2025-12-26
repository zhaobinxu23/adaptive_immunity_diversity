function F=pathway_model_with_drug_inhib(t,y,para)
for i = 1: 6
    y(i) = max(0,y(i));
end


% 
F(1,1) = -para(1)*y(1)*y(2)+para(2)*y(3)-para(3)*y(1);

F(2,1) = -para(1)*y(1)*y(2)+para(2)*y(3)+para(5)*para(11)*y(3)+para(5)*y(2)-para(6)*y(2)*y(4)+para(7)*y(5);

F(3,1) = para(1)*y(1)*y(2) - para(2)*y(3) - para(6)*y(3)*y(4) + para(7)*y(6);

F(4,1) = -para(6)*(y(2)+y(3))*y(4) + para(7)*(y(5)+y(6)) + para(8) - para(9)*y(4) + para(10)*(y(5)+y(6));

F(5,1) = para(6)*y(2)*y(4) - para(7)*y(5)-para(4)*y(5);

F(6,1) = para(6)*y(3)*y(4) - para(7)*y(6)-para(4)*y(6);

end





