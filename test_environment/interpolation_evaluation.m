%% One Dimensional Interpolation Decision
% David Dobbie

clc
clear

set(0,'DefaultAxesXGrid','on')

loadedExperimental = load('testing_data\T2Distributions.mat');
T2_dist_cells = loadedExperimental.T2_distribution;
extractedT2Data = cell2mat(T2_dist_cells);
experimentalfT2 = extractedT2Data(:,2:2:end);
experimentalT2Axis = extractedT2Data(:,1)/1000;


set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)





Ny = 100; % number of bins in relaxation time grids
tE = 200e-6;  % sample interval
T2 = logspace(log10(2*tE),3,Ny); %form T2 domain, use log since will be small
T2_coarse = logspace(log10(2*tE),1,10);

true_raw = csvread('validator_models/model3_true.csv');

%true_raw = [experimentalT2Axis ,experimentalfT2(:,2)]

%{
interpol_linear = interp1(true_raw(:,1), true_raw(:,2), T2, 'linear',0);
interpol_nearest = interp1(true_raw(:,1), true_raw(:,2), T2, 'nearest',0);
interpol_next = interp1(true_raw(:,1), true_raw(:,2), T2, 'next',0);
%interpol_previous = interp1(experimentalT2Axis, experimentalfT2(:,1), T2, 'previous',0);
interpol_pchip = interp1(true_raw(:,1), true_raw(:,2), T2, 'pchip',0);
interpol_makima = interp1(true_raw(:,1), true_raw(:,2), T2, 'makima',0);
interpol_spline = interp1(true_raw(:,1), true_raw(:,2), T2, 'spline',0);
%}

coarse_f = interp1(true_raw(:,1), true_raw(:,2), T2_coarse, 'linear',0);
true_raw = coarse_f;

interpol_linear = interp1(T2_coarse, coarse_f, T2, 'linear',0);
interpol_nearest = interp1(T2_coarse, coarse_f, T2, 'nearest',0);
interpol_next = interp1(T2_coarse, coarse_f, T2, 'next',0);
%interpol_previous = interp1(experimentalT2Axis, experimentalfT2(:,1), T2, 'previous',0);
interpol_pchip = interp1(T2_coarse, coarse_f, T2, 'pchip',0);
interpol_makima = interp1(T2_coarse, coarse_f, T2, 'makima',0);
interpol_spline = interp1(T2_coarse, coarse_f, T2, 'spline',0);
    

figure(1)
clf
subplot(3,2,1)
hold on
p1 = plot(T2_coarse, coarse_f);
p1.LineWidth = 2;
p1.LineStyle = '-';

subplot(3,2,2)
hold on
p2 = plot(T2_coarse, coarse_f);
p2.LineWidth = 2;
p2.LineStyle = '-';

subplot(3,2,3)
hold on
p3 = plot(T2_coarse, coarse_f);
p3.LineWidth = 2;
p3.LineStyle = '-';
subplot(3,2,4)
hold on
p4 = plot(T2_coarse, coarse_f);
p4.LineWidth = 2;
p4.LineStyle = '-';
subplot(3,2,5)
hold on
p5 = plot(T2_coarse, coarse_f);
p5.LineWidth = 2;
p5.LineStyle = '-';
subplot(3,2,6)
hold on
p6 = plot(T2_coarse, coarse_f);
p6.LineWidth = 2;
p6.LineStyle = '-';
set(gca, 'XScale', 'log') 



subplot(3,2,1)
p2 = plot(T2, interpol_linear);
p2.LineWidth = 2;
p2.LineStyle = '--';
set(gca, 'XScale', 'log') 
xlim([4e-4 max(T2)])
legend('True','Linear Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

subplot(3,2,2)
p3 = plot(T2, interpol_nearest);
p3.LineWidth = 2;
p3.LineStyle = ':';
set(gca, 'XScale', 'log') 
xlim([4e-4 max(T2)])
legend('True','Nearest Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

subplot(3,2,3)
p6 = plot(T2, interpol_pchip);
p6.LineWidth = 2;
p6.LineStyle = '-.';
set(gca, 'XScale', 'log')
xlim([4e-4 max(T2)])
legend('True','PCHIP Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

subplot(3,2,4)
p7 = plot(T2, interpol_makima);
p7.LineWidth = 2;
p7.LineStyle = '--';
set(gca, 'XScale', 'log') 
xlim([4e-4 max(T2)])
legend('True','Makima Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

subplot(3,2,5)
p8 = plot(T2, interpol_spline);
p8.LineWidth = 2;
p8.LineStyle = ':';
set(gca, 'XScale', 'log') 
xlim([4e-4 max(T2)])
legend('True','Spline Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

subplot(3,2,6)
p8 = plot(T2, interpol_next);
p8.LineWidth = 2;
p8.LineStyle = ':';

set(gca, 'XScale', 'log')
xlim([4e-4 max(T2)])
legend('True','Next Point Interp.','Location','NorthEast')
xlabel('$T_2$')
ylabel('$f(T_2)$')

hold off
set(gca, 'XScale', 'log') 
