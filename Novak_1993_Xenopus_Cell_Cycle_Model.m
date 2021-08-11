%% Novak-Tyson (1993) Model of Xenopus cell cycle with/without DNA Unreplicated Material
% As described in Sible & Tyson (2007)
% Created by: Icahn School Of Medicine at Mount Sinai, Professor Eric Sobie

% August 10th, 2021
% Modified by Santiago Taguado 


%% With DNA Replication

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1:  Define constants 

global k1 k3
k1 = 1 ;
k3 = 0.005 ;

global ka Ka Kb kc Kc kd Kd
global ke Ke Kf kg Kg kh Kh
ka = 0.02 ;
Ka = 0.1 ;
Kb = 1 ;
kc = 0.13 ;
Kc = 0.01 ;
kd = 0.13 ;
Kd = 1 ;
ke = 0.02 ;
Ke = 1 ;
Kf = 1 ;
kg = 0.02 ;
Kg = 0.01 ;
kh = 0.15 ;
Kh = 0.01 ;


global v2_1 v2_2 v25_1 v25_2 vwee_1 vwee_2
v2_1 = 0.005 ;
v2_2 = 0.25 ;
% % These modified to increase bistability range
% v25_1 = 0.017 ;
% v25_2 = 0.17 ;
v25_1 = 0.5*0.017 ;
v25_2 = 0.5*0.17 ;
vwee_1 = 0.01 ;
vwee_2 = 1 ;

global CDK_total cdc25_total wee1_total IE_total APC_total PPase DNA_total

CDK_total = 100 ;
% The modification of DNA can affect the cyclin and MPF oscillations
% significantly.
DNA_total = 6;
% % This modified to increase bistability range
%cdc25_total = 1 ;
cdc25_total = 5;
wee1_total = 1 ;
IE_total = 1 ;
APC_total = 1 ;
PPase = 1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2:  Define simulation time 

tlast = 1500 ; % min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3:  Initial conditions 

cyclin = 0 ;
MPF = 0 ;
preMPF = 0 ;
cdc25P = 0 ;
wee1P = 0 ;
IEP = 1 ;
APC = 1 ;

statevar_i = [cyclin,MPF,preMPF,cdc25P,wee1P,IEP,APC] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 4:  Run it!

[time,statevars] = ode15s(@dydt_novak_modified,[0,tlast],statevar_i) ;

cyclin = statevars(:,1) ;
MPF = statevars(:,2) ;
preMPF = statevars(:,3) ;
cdc25P = statevars(:,4) ;
wee1P = statevars(:,5) ;
IEP = statevars(:,6) ;
APC = statevars(:,7) ;

cyclin_tot = cyclin + MPF + preMPF ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5:  Plot/save results

% % Note:  this implementation only plots stable limit cycles
% % This ignores initial condition-dependent effects
dices = find(time > 1000) ;
time = time(dices) - time(dices(1)) ;
cyclin_tot = cyclin_tot(dices) ;
MPF = MPF(dices) ;


% kb and kf effects
DNA = DNA_total - cyclin;
kf = (DNA) + 0.1;
kb = (DNA) + 0.1;

% Plot of kb versus kf
figure
hold on
plot(kb,kf,'k','LineWidth',2.25)
xlabel('kb')
ylabel('kf')
title('Rate Constants Comparison')

% Plot of Time vs. Total Cyclin & Maturaing Promoting Concentration
figure
hold on
plot(time,cyclin_tot,'k') 
plot(time,MPF,'r')
xlabel('Time (min)')
ylabel('[Cyclin_T_o_t_a_l]/[MPF]')
title('With Replicated DNA: Total Cyclin [Cyclin_T_o_t_a_l] and Maturating Promoting Factor Concentration [MPF]')
legend('[Cyclin]_T_o_t_a_l','[MPF]','location','northeast')

% Plot of MPF vs. Cyclin total
figure
handle1 = gcf;
hold on
plot(MPF,cyclin_tot,'g','LineWidth',2.25)
xlabel('[MPF]')
ylabel('[Cyclin_t_o_t_a_l]')
set(gca,'TickDir','Out')
title('Bifurcatin Diagram: Maturating Promoting Factor Concentration [MPF] versus Cyclin')

% Plot of MPF vs. Wee1P/cdc25P Unreplicated DNA
figure
handle2 = gcf;
hold on
plot(MPF,(wee1P(dices)),'r');
plot(MPF,(cdc25P(dices)),'b');
xlabel('MPF')
ylabel('Enzyme')
title('Wee1P/Cdc25P versus Maturating Promoting Factor Concentration')


%% Without DNA Replication

%%% Novak-Tyson (1993) model of Xenopus cell cycle
%%% As described in Sible & Tyson (2007)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1:  Define constants 

global k1 k3
k1 = 1 ;
k3 = 0.005 ;

global ka Ka kb Kb kc Kc kd Kd
global ke Ke kf Kf kg Kg kh Kh
ka = 0.02 ;
Ka = 0.1 ;
kb = 0.1 ;
Kb = 1 ;
kc = 0.13 ;
Kc = 0.01 ;
kd = 0.13 ;
Kd = 1 ;
ke = 0.02 ;
Ke = 1 ;
kf = 0.1 ;
Kf = 1 ;
kg = 0.02 ;
Kg = 0.01 ;
kh = 0.15 ;
Kh = 0.01 ;

global v2_1 v2_2 v25_1 v25_2 vwee_1 vwee_2
v2_1 = 0.005 ;
v2_2 = 0.25 ;
% % These modified to increase bistability range
% v25_1 = 0.017 ;
% v25_2 = 0.17 ;
v25_1 = 0.5*0.017 ;
v25_2 = 0.5*0.17 ;
vwee_1 = 0.01 ;
vwee_2 = 1 ;

global CDK_total cdc25_total wee1_total IE_total APC_total PPase
CDK_total = 100 ;
% % This modified to increase bistability range
% cdc25_total = 1 ;
cdc25_total = 5 ;
wee1_total = 1 ;
IE_total = 1 ;
APC_total = 1 ;
PPase = 1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2:  Define simulation time 

tlast = 1500 ; % min

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3:  Initial conditions 

cyclin = 0 ;
MPF = 0 ;
preMPF = 0 ;
cdc25P = 0 ;
wee1P = 0 ;
IEP = 1 ;
APC = 1 ;

statevar_i = [cyclin,MPF,preMPF,cdc25P,wee1P,IEP,APC] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 4:  Run it!

[time,statevars] = ode15s(@dydt_novak,[0,tlast],statevar_i) ;

cyclin = statevars(:,1) ;
MPF = statevars(:,2) ;
preMPF = statevars(:,3) ;
cdc25P = statevars(:,4) ;
wee1P = statevars(:,5) ;
IEP = statevars(:,6) ;
APC = statevars(:,7) ;

cyclin_tot = cyclin + MPF + preMPF ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5:  Plot/save results
% % Note:  this implementation only plots stable limit cycles
% % This ignores initial condition-dependent effects
dices = find(time > 1000) ;
time = time(dices) - time(dices(1)) ;
cyclin_tot = cyclin_tot(dices) ;
MPF = MPF(dices) ;

figure
hold on
plot(time,cyclin_tot,'k') 
plot(time,MPF,'r')
xlabel('Time (min)')
ylabel('[Cyclin_T_o_t_a_l]/[MPF]')
title('Without Replicated DNA: Total Cyclin [Cyclin_T_o_t_a_l] and Maturating Promoting Factor Concentration [MPF]')
legend('[Cyclin]_T_o_t_a_l','[MPF]','location','northeast')


% Run while previous section figures are open
figure(handle1)
hold on
plot(MPF,cyclin_tot,'k','LineWidth',2.25)
set(gca,'TickDir','Out')
legend('With Unreplicated DNA','Without Unreplicted DNA','location','southeast')


% Run while previous section figures are open
figure(handle2)
hold on
plot(MPF,(wee1P(dices)),'ro');
plot(MPF,(cdc25P(dices)),'bo');
legend('wee1P with DNA','CDC25P with DNA','wee1P control','CDC25P control','location','northeast')
