function F=pathway_model_4_antibody_immune_cellular_humoral(t,y,para)

% for i = 1: 4
%     y(i) = max(0,y(i));
% end

F(1,1) = -para(1)*y(1)*y(2)+para(2)*y(3)+para(3)*y(3)+para(7)-para(8)*y(1);

F(2,1) = para(5)*y(2)-para(1)*y(1)*y(2)+para(2)*y(3);

F(3,1) = para(1)*y(1)*y(2)-para(2)*y(3)-para(4)*y(3);


F(4,1) = para(6)*(y(2)+y(3))+para(9)-para(10)*y(4);
% F(4,1) = para(6)*(y(3))+para(9)-para(10)*y(4);
end




%created by the program testexcel_IL
