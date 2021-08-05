%%
% By Santiago Taguado
% August 4th, 2021
% Rate Balance Plots & Bifurfaction Diagrams

%% Rate & Bifurcation Diagram: Linear Feedback plus saturating back reaction effect of Kmb Changes

clear all
clc

Astar = 0:0.01:1 ;
S = 1 ;
kplus = 2 ;
kf = 30 ;
kminus = 5 ;
Kmb = 0:0.01:4;

FR = (kplus*S+kf.*Astar).*(1-Astar);
figure
handle1 = gcf ;
hold on 
plot(Astar,FR,'r','LineWidth',2)
set(gca,'TickDir','Out')

figure
handle2 = gcf ;
hold on

for i=1:length(Kmb)
  BR = kminus*((Astar)./(Astar + Kmb(i))) ; 
  figure(handle1)
  plot(Astar,BR,'b','LineWidth',2)
  
  crossings = [] ;
  difference = FR-BR ;
  for iii=2:length(FR)
    if (sign(difference(iii)) ~= sign(difference(iii-1)))
      crossings = [crossings,iii] ;
    end
  end
  figure(handle2)
  plot(Kmb(i),Astar(crossings),'bo')
  
end

figure(handle1)
set(gca,'TickDir','Out')
xlabel('[ A* ]')
ylabel('Rates')
title('Forward Reaction (FR) vs. Back Reaction (BR) Diagram')
legend('FR','BR','location','northeast')

figure(handle2)
set(gca,'TickDir','Out')
xlabel('Kmb')
ylabel('Steady-state [A*]/[A]')
title('Bifurcation Diagram: Kmb vs. Steady State [ A* ]')


%% Rate & Bifurcation Diagram: Linear Feedback plus saturating back reaction effect of Stimulus Changes

clear all
clc

Astar = 0:0.01:1 ;
S = 0:0.01:4 ;
kplus = 2 ;
kf = 30 ;
kminus = 5 ;
Kmb = 1;

BR = kminus*((Astar)./(Astar + Kmb)); 
figure
handle3 = gcf ;
hold on 
plot(Astar,BR,'r','LineWidth',2)
set(gca,'TickDir','Out')

figure
handle4 = gcf ;
hold on

for i=1:length(S)
  FR = (kplus*S(i)+kf.*Astar).*(1-Astar);
  figure(handle3)
  plot(Astar,FR,'b','LineWidth',2)
  
  crossings = [] ;
  difference = FR-BR ;
  for iii=2:length(FR)
    if (sign(difference(iii)) ~= sign(difference(iii-1)))
      crossings = [crossings,iii] ;
    end
  end
  figure(handle4)
  plot(S(i),Astar(crossings),'bo')
  
end

figure(handle3)
set(gca,'TickDir','Out')
xlabel('[A*]')
ylabel('Rates')
title('Linear Feedback plus Saturating Back Reaction Rate Diagram')
legend('BR','FR','location','northeast')


figure(handle4)
set(gca,'TickDir','Out')
xlabel('Stimulus [S]')
ylabel('Steady-state [A*]/[A]')
title('Linear Feedback plus saturating back reaction Bifurcation Diagram')


%% Rate & Bifurcation Diagram: Ultrasensitive Autocatalytic Feedback Loop & Saturating Back Reaction effect of Stimulus Change

clear all
clc

Astar = 0:0.01:1 ;
S = 0:0.02:4 ;
kplus = 0.5 ;
kf = 30 ;
Kmb = 0.5 ;
kminus = 5 ;
Kmf = 0.5;
h = 4 ; % exponent

BR = kminus*((Astar)./(Astar + Kmb));
figure
handle5 = gcf ;
hold on 
plot(Astar,BR,'r','LineWidth',2)
set(gca,'TickDir','Out')

figure
handle6 = gcf ;
hold on

for i=1:length(S)
  FR = (kplus*S(i)+kf*(Astar.^h./(Astar.^h+Kmf^h))).*(1-Astar) ;
  figure(handle5)
  plot(Astar,FR,'b','LineWidth',2)
  
  crossings = [] ;
  difference = FR-BR ;
  for iii=2:length(FR)
    if (sign(difference(iii)) ~= sign(difference(iii-1)))
      crossings = [crossings,iii] ;
    end
  end
  figure(handle6)
  plot(S(i),Astar(crossings),'bo')
  
end

figure(handle5)
set(gca,'TickDir','Out')
xlabel('[A*]/[A]')
ylabel('Rates')
title('Ultrasensitive Autocatalytic Feedback Loop & Saturating Back Reaction Rate Diagram')
legend('BR','FR','location','northeast')


figure(handle6)
set(gca,'TickDir','Out')
xlabel('Stimulus [S]')
ylabel('Steady-state [A*]/[A]')
title('Ultrasensitive Autocatalytic Feedback Loop & Saturating Back Reaction Bifurcation Diagram')


%% Rate & Bifuraction Diagram: Ultrasensitive Autocatalytic Feedback Loop effect of Stimulus Change

clear all
clc

Astar = 0:0.01:1 ;
S = 0:0.02:4 ;
kplus = 0.5 ;
kf = 30 ;
Kmb = 0.5 ;
kminus = 5 ;
Kmf = 0.5;
h = 4 ; % exponent

BR = kminus*(Astar);
figure
handle9 = gcf ;
hold on 
plot(Astar,BR,'r','LineWidth',2)
set(gca,'TickDir','Out')

figure
handle10 = gcf ;
hold on

for i=1:length(S)
  FR = (kplus*S(i)+kf*(Astar.^h./(Astar.^h+Kmf^h))).*(1-Astar) ; 
  figure(handle9)
  plot(Astar,FR,'b','LineWidth',2)

  crossings = [] ;
  difference = FR-BR ;
  for iii=2:length(FR)
    if (sign(difference(iii)) ~= sign(difference(iii-1)))
      crossings = [crossings,iii] ;
    end
  end
  figure(handle10)
  plot(S(i),Astar(crossings),'bo')
  
end

figure(handle9)
set(gca,'TickDir','Out')
xlabel('[A*]/[A]')
ylabel('Rates')
title('Ultrasensitive Autocatalytic Feedback Loop Rate Diagram')
legend('BR','FR','location','northeast')


figure(handle10)
set(gca,'TickDir','Out')
xlabel('Stimulus [S]')
ylabel('Steady-state [A*]/[A]')
title('Ultrasensitive Autocatalytic Feedback Loop Bifurcation Diagram')


%% Nullclines: Simple Two-Variable Model of E.Coli Operon

clear all
clc

% Differential Equation 1
% dl = b*l_ext*lacY - gamma*l;
    
% Differential Equation 2
% dlacY = lambda + p.*(l.^4./(l.^4 + l_0^4)) - alpha*lacY; 

b = 1;
gamma = 1;
lambda = 0.2;
l_0 = 4;
p = 4;
alpha = 1;
l_ext = 2.5;
l = 0:0.1:30;

% First Equation Rearranged for lacY
lac_y_1 = ((gamma.*l)/(b*l_ext));

% Second Equation Rearranged
lac_y_2 = ((lambda + p.*(l.^4./(l.^4 + l_0^4))))/alpha;

figure
hold on
plot(l,lac_y_1,'b')
plot(l,lac_y_2,'r')
ylabel('Lactose Permease Gene Concentration')
xlabel('Lactose Concentration')
title('Nullcline Diagram: Two-Variable Model of E. Coli Operon')
legend('Differential Eq. 1','Differential Eq. 2','location','northeast')


