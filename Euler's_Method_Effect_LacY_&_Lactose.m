%% Euler's Method: LacY and Lactose Initial Concentration Biostability

clear all
clc

% Defining Constants
b = 1;
gamma = 1;
lambda = 0.2;
l_0 = 4;
p = 4;
alpha = 1;
l_ext = 2.5;

%%% Defining time step, simulation time, initialize matrices 
dt    = 0.01 ; % s 
tlast = 20 ;  % s
iterations = fix(tlast/dt) ;
time = dt*(0:iterations-1) ;

figure
handle1 = gcf
hold on
 
for ii = 1:4
    
l = [8,3,3,2];
lacY = [3,1.3,1.2,1] ;
l = l(ii);
lacY = lacY(ii);


for i = 1:iterations 
  l_tot(i) = l;
  lacY_tot(i) = lacY;
  dl = b*l_ext*lacY - gamma*l;
  dlacY = lambda + p.*(l.^4./(l.^4 + l_0^4)) - alpha*lacY;    
  l = l + dt*dl ;
  lacY = lacY + dt*dlacY ;
end 

figure(handle1)
plot(time,l_tot,'b')
plot(time,lacY_tot,'r')

end

figure(handle1)
title('Simulation of Temporal Evolution of Lactose [L] & Lactose Permease Gene [lacY]')
xlabel('Time (seconds)')
ylabel('Concentration of [L]/ Concentration of [LacY]')
legend('[L]','[lacY]','location','northwest')


%% Bifurcation Diagram: LacY vs. L_ext

clear all
clc

% Defining Constants

b = 1;
gamma = 1;
lambda = 0.2;
l_0 = 4;
p = 4;
alpha = 1;

%%% Defining time step, simulation time, initialize matrices 
dt    = 0.01 ; % s 
tlast = 20.01 ;  % s
iterations = fix(tlast/dt) ;
time = dt*(0:iterations-1) ;

ii = 0;
for j = 1:0.1:7

l = [8];
lacY = [3];
l_ext = j;
ii = ii + 1

for i = 1:iterations
  l_tot(i) = l;
  lacY_tot(i) = lacY;
  dl = b.*l_ext*lacY - gamma*l;
  dlacY = lambda + p.*(l.^4./(l.^4 + l_0^4)) - alpha*lacY;    
  l = l + dt*dl ;
  lacY = lacY + dt*dlacY ;
end

l_fin(ii) = lacY_tot(iterations)

end

ii = 0
for j = 1:0.1:7

l = [2];
lacY = [1];
l_ext = j;
ii = ii + 1;

for i = 1:iterations
  l_tot(i) = l;
  lacY_tot(i) = lacY;
  dl = b.*l_ext*lacY - gamma*l;
  dlacY = lambda + p.*(l.^4./(l.^4 + l_0^4)) - alpha*lacY;    
  l = l + dt*dl ;
  lacY = lacY + dt*dlacY ;
end

l_fin2(ii) = lacY_tot(iterations)

end

l_ext = 1:0.1:7;

figure
hold on
plot(l_ext,l_fin,'r')
plot(l_ext,l_fin2,'b')
title('Lactose Extracellular Concentration Effect on Lactose Permease Gene Concentration')
legend(' Initial Conditions [L] = 8 [LacY] = 3','Initial Conditions [L] = 2 [LacY] = 1','location','southeast')
xlabel ('Lactose Extracellular Conc. [Lext]')
ylabel('Lactose Permease Gene Conc. [LacY]')



%% Testing, if Plot is correct

clc
clear all

b = 1;
gamma = 1;
lambda = 0.2;
l_0 = 4;
p = 4;
alpha = 1;
% Individually changing l
l = 2;
% Individually chaning lac_y
lacY = 1;
% Individually changing l_ext
l_ext = 7;


%%% Defining time step, simulation time, initialize matrices 
dt    = 0.01 ; % s 
tlast = 20.01 ;  % s
iterations = fix(tlast/dt) ;
time = dt*(0:iterations-1) ;

for i = 1:iterations 
  l_tot(i) = l;
  lacY_tot(i) = lacY;
  dl = b.*l_ext*lacY - gamma*l;
  dlacY = lambda + p.*(l.^4./(l.^4 + l_0^4)) - alpha*lacY;    
  l = l + dt*dl ;
  lacY = lacY + dt*dlacY ;
end

figure
hold on
plot(time,lacY_tot,'b')
plot(time,l_tot,'r')
xlabel('Time (second)')
ylabel('Concentration of [LacY]/Concentration of [L]')
legend('[LacY]','[L]','location','northwest')

