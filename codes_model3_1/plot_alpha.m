%% D = x0（1）, V = x0(2), V' = x0(3), A = x0(4), C = x0(5)

clc
clear




 for ii = 1:100
     ii
  for jj = 1:100
    for i = 1:100
   

x0(1) = 0; %% A
x0(2) = 500; %% B
x0(3) = 0; %% C1
x0(4) = 0; %% C2
x0(5) = 1;%% V
x0(6) = 0; %% ASC cell


para(1) = 1e-7; % k1
para(2) = 1e-14; % k-1
para(3) = 2+0.01*ii;% k2
para(4) = 1e-3;% alpha
para(5) = 0.005;% d1
para(6) = 0.02;% d2
para(7) = 0.5;% d3
para(8) = 0.5;% d4
para(9) = 0.5+0.001*jj;% k3
para(10) = 1e-7;
para(11) = 0.01; %% decay ratio of ASC cell
para(12) = 1e5;
para(13) = 0.005*i;


 [t,y]=ode15s(@pathway_model_plasma_model,[0 80],x0,[],para);
 dd(i) = max(y(:,5));
 xx(ii,jj) = 0.005*(find(dd(:) == min(dd)));

     end
  end
end










