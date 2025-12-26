function F=pathway_model_4_antibody_immune_res_new(t,y,para)

for i = 1: 14
    y(i) = max(0,y(i));
end

F(1,1) = -para(1)*y(1)*y(2)+para(5)*y(3)+para(10)*y(3)-para(14)*y(1)*y(10)+para(15)*y(11)+para(10)*y(3)+para(13)*y(11)-para(9)*y(1);

F(2,1) = para(11)*y(2)-para(1)*y(1)*y(2)-para(1)*y(1)*y(2)+para(5)*y(3)-para(2)*y(4)*y(2)+para(6)*y(5)-para(3)*y(6)*y(2)+para(7)*y(7)-para(4)*y(8)*y(2)+para(8)*y(9);

F(3,1) = para(1)*y(1)*y(2)-para(5)*y(3)-para(12)*y(3);

F(4,1) = -para(2)*y(4)*y(2)+para(6)*y(5)+para(10)*y(5)-para(14)*y(4)*y(10)+para(15)*y(12)+para(10)*y(5)+para(13)*y(12)-para(9)*y(4);

F(5,1) = para(2)*y(4)*y(2)-para(6)*y(5)-para(12)*y(5);

F(6,1) = -para(3)*y(6)*y(2)+para(7)*y(7)+para(10)*y(7)-para(14)*y(6)*y(10)+para(15)*y(13)+para(10)*y(7)+para(13)*y(13)-para(9)*y(6);


F(7,1) = para(3)*y(6)*y(2)-para(7)*y(7)-para(12)*y(7);


F(8,1) = -para(4)*y(8)*y(2)+para(8)*y(9)+para(10)*y(9)-para(14)*y(8)*y(10)+para(15)*y(14)+para(10)*y(9)+para(13)*y(14)-para(9)*y(8);

F(9,1) = para(4)*y(8)*y(2)-para(8)*y(9)-para(12)*y(9);

F(10,1) = para(16)-para(14)*y(1)*y(10)+para(15)*y(11)-para(14)*y(4)*y(10)+para(15)*y(12)-para(14)*y(6)*y(10)+para(15)*y(13)-para(14)*y(8)*y(10)+para(15)*y(14);

F(11,1) = para(14)*y(1)*y(10)-para(15)*y(11)-para(12)*y(11);

F(12,1) = para(14)*y(4)*y(10)-para(15)*y(12)-para(12)*y(12);

F(13,1) = para(14)*y(6)*y(10)-para(15)*y(13)-para(12)*y(13);

F(14,1) = para(14)*y(8)*y(10)-para(15)*y(14)-para(12)*y(14);

end




%created by the program testexcel_IL
