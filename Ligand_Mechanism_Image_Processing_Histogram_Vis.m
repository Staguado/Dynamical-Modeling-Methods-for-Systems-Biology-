%% 

% By Santiago Taguado
% July 18, 2021
% Ligand vs. Ligand-Receptor Mechanism, Image Data Analysis, & Histogram Visualization

%% Ligand vs. Ligand-Receptor Plot

clear all
clc

colors = 'krgbmc';
R_tot = 20;
L = 0:0.01:200;
Kd = 10:20:90;

figure
hold on

for i = 1:length(Kd)
KD = Kd(i);
LR = (R_tot * L) ./ (KD+L);
LR_dat(:,i) = LR;
plot(L,LR,colors(i))
end

xlabel('[Ligand] (uM)')
ylabel('[Ligand-receptor] (nM')

%% From Imaging to Data Analysis

clear all
clc

% The experiment of Calcium Effect in the Ventricular Myocyte of
% a Rat to see the cell is affected by calcium concentration

% Transpose of Image
imagesc(flash4')
colorbar

% Flash Region
flash_transpose = flash4;
flash = flash_transpose(200:300,:);
F_flash = mean(flash);

% No Flash Region
no_flash_1 = flash_transpose(301:512,:);
F_no_flash = mean(no_flash_1);

% Control Region
con_flash = flash_transpose(70:100,:);
F_knot = mean(con_flash);

% Calculating Calcium Concenrations for Flash
R_flash = F_flash./F_knot;

% Calculating Calcium Concentrations for No Flash
R_noflash = F_no_flash./F_knot;

% Plotting
Kd = 1000;
Ca = 150;
Ca_con_flash = (R_flash .* Kd)./(Kd/Ca - R_flash + 1);
Ca_con_noflash = (R_noflash .* Kd)./(Kd/Ca - R_noflash + 1);
time = 1:1:634;

figure
hold on
plot(time,Ca_con_flash, 'r')
plot(time, Ca_con_noflash,'b')
xlabel('Time(s)')
ylabel('[Ca] nM')
legend('Flash', 'No Flash')

%% Analysis of Data: Comparing alcohol consumption to cancer rates

clear all
clc
load sampledata2.mat

patient_age = data(:,1);
patient_drink = data(:,2);
patients_cancer = data(:,3);

% Patients with cancer
sum(patients_cancer > 0);
% Patients who drink more than 6 and less than 5 drinks per week
sum(patient_drink > 6);
sum(patient_drink < 5)
% Patients who drink more than 5 drinks per week
sum(patient_drink > 5)
% Patients over the age of 50 and under the age of 20
sum(patient_age > 50)
sum(patient_age < 20)

% Histogram of Data

can_ind = find(patients_cancer == 1);
no_can_ind = find(patients_cancer == 0);

patient_age_cancer = patient_age(can_ind);
patient_age_nocancer = patient_age(no_can_ind);

edges = [15 25 35 45 55 65 75];
figure
histogram(patient_age_cancer,edges)

% No cancer Histogram
edges = [15 25 35 45 55 65];
histogram(patient_age_nocancer,edges)

% Relating Patient Age and cancer
sum(patient_age_cancer < 30)
mean(patient_age_cancer > 45)

% Relating Patient Age and cancer (2)
mean(patient_age_cancer > 35)

% Relating Patient Age and cancer (3)
mean(patient_age < 25)*100

% Cancer Rates for patients who drink more than 3 times a week
ind = find(patient_drink > 3)
mean(patients_cancer(ind) == 1)*100
