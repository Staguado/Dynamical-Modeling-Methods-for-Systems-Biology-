%%

% By Santiago Taguado Menza
% July 31st, 2021
% Model of yeast glycolytic oscillations developed by Bier et al (1996)

%% Glucose & ATP Concentration Unsteady System

V_in = 0.36;
k1 = 0.02 ;
kp = 6 ;
Km = 12;

tlast = 500;
dt = 0.05;
 
iterations = round(tlast/dt) ; 
ATP_tot = zeros(iterations,1) ;
G_tot = zeros(iterations,1);

% Initial Conditions
ATP = 4;
G = 3;

for i = 1:iterations 
    
    ATP_tot(i) = ATP ;
    G_tot(i) = G;
    datpdt = 2*k1*G*ATP - ((kp*ATP)/(ATP+Km));
    dgdt =  V_in - k1*G*ATP;
    ATP = ATP + datpdt*dt;
    G = G + dgdt*dt;
end 
        
time = dt*(0:iterations-1)' ;

figure
hold on
plot(time, ATP_tot, 'b')
plot(time, G_tot,'r')
title('Glucose [G] & Adenosine 5-triphosphate [ATP] Concentrations')
xlabel('Time (s)')
ylabel('[G]/[ATP]')
legend('[ATP]','[G]','location','northeast')

%% Glucose & ATP Concentration Sustained Oscillatory Behavior

clear all
clc

V_in = 0.1;
k1 = 0.01 ;
kp = 3 ;
Km = 12;

tlast = 500;
dt = 0.05;
 
iterations = round(tlast/dt) ; 
ATP_tot = zeros(iterations,1) ;
G_tot = zeros(iterations,1);

% Initial Conditions
ATP = 4;
G = 3;

for i = 1:iterations 
    
    ATP_tot(i) = ATP ;
    G_tot(i) = G;
    datpdt = 2*k1*G*ATP - ((kp*ATP)/(ATP+Km));
    dgdt =  V_in - k1*G*ATP;
    ATP = ATP + datpdt*dt;
    G = G + dgdt*dt;
end 
        
time = dt*(0:iterations-1)' ;

figure
hold on
plot(time, ATP_tot, 'b')
plot(time, G_tot,'r')
title('Glucose [G] & Adenosine 5-triphosphate [ATP] Concentrations')
xlabel('Time (s)')
ylabel('[G]/[ATP]')
legend('[ATP]','[G]','location','northeast')

%% Glucose & ATP Concentration Damped Oscillatory Behavior: Concentration vs. Time & [ATP] vs. [G]

clear all
clc

V_in = 0.3;
k1 = 0.02 ;
kp = 6 ;
Km = 18;


tlast = 2000;
dt = 0.05;
 
iterations = round(tlast/dt) ; 
ATP_tot = zeros(iterations,1) ;
G_tot = zeros(iterations,1);

ATP = 4;
G = 3;

    for i = 1:iterations 
    ATP_tot(i) = ATP ;
    G_tot(i) = G;
    datpdt = 2*k1*G*ATP - ((kp*ATP)/(ATP+Km));
    dgdt =  V_in - k1*G*ATP;
    ATP = ATP + datpdt*dt;
    G = G + dgdt*dt;
    end 


    
time = dt*(0:iterations-1)' ;



figure
hold on
plot(time, ATP_tot, 'b')
plot(time, G_tot,'r')
title('Glucose [G] & Adenosine 5-triphosphate [ATP] Concentrations vs. time')
xlabel('Time (s)')
ylabel('[G]/[ATP]')
legend('[ATP]','[G]','location','northeast')

figure
hold on
plot(ATP_tot,G_tot,'g')
title('Glucose [G] & Adenosine 5-triphosphate [ATP] Concentrations')
xlabel('[ATP]')
ylabel('[G]')

%% Bifurcation Diagram of [G] & [ATP] vs. Glucose Transport Rate [Vin]

clear all
clc

k1 = 0.02 ;
kp = 6 ;
Km = 13;

tlast = 2000;
dt = 0.05;
 
iterations = round(tlast/dt) ; 
ATP_tot = zeros(iterations,1) ;
G_tot = zeros(iterations,1);

ATP = 4;
G = 3;
Vin = 0.1:0.1:1.4;

for j = 1:14
    
    V_in = Vin(j);
   
    for i = 1:iterations 
    ATP_tot(i) = ATP ;
    G_tot(i) = G;
    datpdt = 2*k1*G*ATP - ((kp*ATP)/(ATP+Km));
    dgdt =  V_in - k1*G*ATP;
    ATP = ATP + datpdt*dt;
    G = G + dgdt*dt;
    end 
    
    Gmax(j) = max(G_tot(750:1000));
    Gmin(j) = min(G_tot(750:1000));
    Amax(j) = max(ATP_tot(750:1000));
    Amin(j) = min(ATP_tot(750:1000));
    
end

    
time = dt*(0:iterations-1)' ;

figure
hold on
plot(Vin,Gmax,'b')
plot(Vin,Gmin,'b-.')
plot(Vin,Amax,'r')
plot(Vin,Amin,'r-.')
title('Bifurcation Diagram: Glucose [G] & Adenosine 5-triphosphate [ATP] Concentrations vs. Glucose Transport Rate (Vin)')
xlabel('V_i_n')
ylabel('[G]/[ATP]')
legend('[G] maxima','[G] minima','[ATP] maxima','[ATP] minima','location','northeast')
