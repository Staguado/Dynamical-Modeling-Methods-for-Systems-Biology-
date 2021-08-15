%%  Hudgkin - Huxley Model
% Created by Icahn School of Medicine at Mount Sinai - Prof. Eric Sobie
% Modified by Santiago Taguado on August 14th, 2021

%% Effects of Second Stimulus on the Hudgkin - Huxley Model
%    t                   time                    ms
%    V                   membrane potantial      mV
%    INa,IK,Il,Iion      ionic current           uA/cm2
%    Cm                  capacitance             uF/cm2

%%%%%%%%%%%%%%%%%%%%%
% Defining Constants
%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

% Physical constants
global F R T RTF 
F = 96.5;                   % Faraday constant, coulombs/mmol
R = 8.314;                  % gas constant, J/K
T_celsius = 6.3;            % Temperature in celsius
T = 273 + T_celsius ;       % absolute temperature, K 

RTF = R*T/F ;

% default concentrations for squid axon in sea water - mmol/l
global Nao Ko Nai Ki 
Nao = 491 ;
Ko = 20 ;
Nai = 50 ;
Ki = 400 ;

% Cell constant
global Cm 
Cm = 1 ;                            % membrane capacitance, uF/cm^2;

% Maximum channel conductances -- mS/cm^2
global GNa GK Gl ENa EK El 
GNa = 120;
GK = 36;
Gl = 0.3;

% Nernst potentials -- mV
ENa = RTF*log(Nao/Nai);
EK = RTF*log(Ko/Ki);
El =  -49;

%%%%%%%%%%%%%%%%%%%%%%
% Defining Simulation lenght and Creating For Loops
%%%%%%%%%%%%%%%%%%%%%

tend = 50 ; % end of simulation, ms
handle1 = gcf;

for j = 1:4

% Stimulus Delay    
stimdelay = 1 ;
% Stimulus Duration
stimdur = 5 ;
% Stimulus Strength & Amplication
stim_amp = -10 ;

% Start of 1st Stimulus
stim_start = stimdelay ;
% Start of 1st Stimulus
stim_end = stimdelay + stimdur;

% Start Times of 2nd Stimulus
start_time = [8,13,6.01,7];
% Finish Times of 2nd Stimulus
fin_time = [19,20,21,11];

% 2nd Stimulus Start Time
sec_stim_start = start_time(j);
% 2nd Stimuls Finish Time
sec_stim_end = fin_time(j);

simints = 5;

% Stimulus Intervals
intervals(1,:) = [0,stim_start] ;
intervals(2,:) = [stim_start,stim_end] ;
intervals(3,:) = [stim_end,sec_stim_start] ;
intervals(4,:) = [sec_stim_start,sec_stim_end];
intervals(5,:) = [sec_stim_end,tend];

% Defining Stimulus Application
Istim(1) = 0 ;
Istim(2) = stim_amp ;
Istim(3) = 0 ;
Istim(4) = stim_amp;
Istim(5) = 0;

% Initial Conditions
V = -60 ;
m = 0 ;
h = 0.6 ;
n = 0.3 ;

% State Variables
statevar_i = [V,m,h,n] ;

% % Simulate 60 seconds at rest before stimulus applied
[post,posstatevars] = ode15s(@dydt_hh,[0,60000],statevar_i,[],0) ;
statevar_i = posstatevars(end,:) ;

% Application of Ordinary Differential Equation Solver
t = 0 ;
statevars = statevar_i ;
for i=1:simints
    
  [post,posstatevars] = ode15s(@dydt_hh,intervals(i,:),statevar_i,[],Istim(i)) ;
  t = [t;post(2:end)] ;
  statevars = [statevars;posstatevars(2:end,:)] ;
  statevar_i = posstatevars(end,:) ;
  
end

outputcell = num2cell(statevars,1) ;

[V,m,h,n] = deal(outputcell{:}) ;

gNa = GNa*m.^3.*h;
INa = gNa.*(V-ENa);

gK = GK*n.^4; 
IK = gK.*(V-EK);

Il = Gl*(V-El) ;
Iion = INa + IK + Il ;

% Ploting Results
figure(handle1)
hold on
plot(t,V)
set(gca,'TickDir','Out')
xlabel('time (ms)')
ylabel('V_m (mV)')
title('Location of Second Stimulus Effect on Action Potentials')
legend('[8ms-13ms]','[13ms-20ms]','[6ms-21ms]','[7ms-11ms]','location','northeast');

end