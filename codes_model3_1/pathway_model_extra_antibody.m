function F=pathway_model_extra_antibody(t,y,para)
% the rate constant parameters include 5 antibodies
%% para(1) = 2e-5; para(2) = 1e-12; para(3) = 10; para(4) = 0.98; para(5) = 1.2
for i = 1:5
    y(i) = max(0,y(i));
end



    
F(1,1) = -para(1)*y(1)*y(3)+para(2)*y(4)+para(3)*y(4);
F(2,1) = -para(6)*y(2)*y(3)+para(7)*y(5)-para(8)*y(2);
F(3,1) = -para(1)*y(1)*y(3)+para(2)*y(4)-para(6)*y(2)*y(3)+para(7)*y(5)+para(5)*y(3);
F(4,1) = para(1)*y(1)*y(3)-para(2)*y(4)-para(4)*y(4);
F(5,1) = para(6)*y(2)*y(3)-para(7)*y(5)-para(4)*y(5);


end