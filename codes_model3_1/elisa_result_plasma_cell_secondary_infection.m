mutant_strain_lead_to_secondary_infection;

for i = 1:100
data_new_IgM(i,:) = interp1(t_new,y_new(:,i+100),(1000:1:1200));
data_new_IgG(i,:) = interp1(t_new,y_new(:,i+300),(1000:1:1200));
end

clear x0




for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 101:200
    x0(i) = data_new_IgM(i-100,k)/1e4;
end

for i = 301:400
    x0(i) = data_new_IgG(i-300,k)/1e4;
end

x0(1402) = 1e16;%% virus



para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;


[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_elisa(k) = 1e16-interp1(t,zz(:,1402),2);

end


for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 101:200
    x0(i) = data_new_IgM(i-100,k)/1e4;
end



x0(1402) = 1e16;%% virus






para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;



[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_IgM(k) = 1e16-interp1(t,zz(:,1402),2);

end


%% 
for k = 1: 201
    k

    x0 = zeros(1402,1);

for i = 301:400
    x0(i) = data_new_IgG(i-300,k)/1e4;
end



x0(1402) = 1e16;%% virus






para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 0; % replenish constant pi 1 
para(4) = 0;% replenish constant pi 2
para(5) = 0; % decay constant of BCR IgM
para(6) = 0;% decay constant of BCR IgG
para(7) = 0; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 0;% k2' feedback constant on PC cell regeneration
para(9) = 0; % decay constant of plasma Cell IgM
para(10) = 0;%% decay constant of plasma Cell IgG
para(11) = 0;% production constant of IgM
para(12) = 0;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.05;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 0; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;



[t, zz]=ode15s(@pathway_model_many_antibody_immune_include_plasma_elisa,[0 2],x0,[],para,para_new,para_new_1);

out_overall_IgG(k) = 1e16-interp1(t,zz(:,1402),2);

end

plot(out_overall_elisa)
hold on
plot(out_overall_IgM)
hold on
plot(out_overall_IgG)