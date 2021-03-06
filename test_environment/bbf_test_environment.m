%% Test Environment for BBF Estimation from T2 Relaxation 
% David Dobbie
% Victoria University of Wellington
% 

%% Aim: 
% This script consolidates all of the testing of the BBF estimator
% techniques in to one system. Each of the estimators are in selfcontained
% function files in this file path.
% This uses the experimental data and cross validation for the entire test
% environment.


clc
clear

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)
set(groot,'defaultLineLineWidth',2)


%% Step 0: intialise variables

%-------- load high quality experimental data

loadedExperimental = load('testing_data\T2Distributions.mat');
T2_dist_cells = loadedExperimental.T2_distribution;
extractedT2Data = cell2mat(T2_dist_cells);
experimentalfT2 = extractedT2Data(:,2:2:end);
experimentalT2Axis = extractedT2Data(:,1)/1000; %assuming all have the same T2 axes spacing


N2 = 2500; % number of data points in each dimension
Ny = 30; % number of bins in relaxation time grids
sing_val=5; %sets how many singular values we compress to
tE = 200e-6;  % sample interval
T2 = logspace(log10(2*tE),2,Ny); %form T2 domain, use log since will be small
tau2 = (1:N2)'*tE;  %forms measurement arrays, time tau2 domains
K2 = exp(-tau2 * (1./T2) ); % simple T2 relaxation kernel


[tot_Prior  del] = size(experimentalfT2');

% interpolate the experimental data to be compatible with domain var
interpol_exp_fT2 = zeros(tot_Prior, Ny);
%interpol_exp_fT2 = interp1(experimentalT2Axis, experimentalfT2, T2, 'pchip',0)';


for idx = 1:tot_Prior
    interpol_exp_fT2(:,idx) = interp1(experimentalT2Axis, experimentalfT2(:,idx), T2, 'pchip',0);
    
    
end



figure(1)
clf
subplot(2,1,1)
plot(experimentalT2Axis, experimentalfT2);
set(gca, 'XScale', 'log') 
subplot(2,1,2)
plot(T2, interpol_exp_fT2);
set(gca, 'XScale', 'log') 




SNR = 10;

%% Step 1: Generate the Bound Fluid Integrals

%Tc = 10;

RMSE_all = zeros(length(T2),5);

for Tc_index = 1:length(T2)
    
%for Tc_index = 10    
    %Tc = 33e-3    
    %Tc = T2(Tc_index);
    Tc = 33e-3;
    
    Tc_index;

    

    bfv_sharp = zeros(Ny ,1);
    for idx = 1:Ny
        if T2(idx)<Tc
            bfv_sharp(idx) = 1;
        end
    end

    % this function comes from Gruber 2013, the tapered area used.
    bfv_tapered = 1 - ((0.7213 / Tc)*tanh(1.572*Tc*(1./T2 + 0.4087/Tc))./ (1./T2 + 0.4087/Tc))';

    figure(44)
    clf
    hold on
    p1 = plot(T2,bfv_sharp);
    p1.LineWidth = 2;
    p2 = plot(T2, bfv_tapered);
    p2.LineWidth = 2;
    hold off
    set(gca, 'XScale', 'log') 
    grid on
    legend('BFV Sharp', 'BFV Tapered')
    xlabel('$T_2$')
    ylabel('$I(T_2)$')
    xlim([min(T2) max(T2)])
    

    %% Step 2: Begin Leave One Out Cross Validation

    test_count = 10;

    bff_est_bayes_sharp_error = zeros(test_count*tot_Prior,1);
    bff_est_bayes_tapered_error = zeros(test_count*tot_Prior,1);
    bff_est_ilt_error = zeros(test_count*tot_Prior,1);
    bff_est_iltx_error = zeros(test_count*tot_Prior,1);

    bff_est_eht_error = zeros(test_count*tot_Prior,1);


    bff_est_bayes_sharp_compute_time = zeros(test_count*tot_Prior,1);
    bff_est_bayes_tapered_compute_time = zeros(test_count*tot_Prior,1);
    bff_est_ilt_compute_time = zeros(test_count*tot_Prior,1);
    bff_est_iltx_compute_time = zeros(test_count*tot_Prior,1);
    bff_est_eht_compute_time = zeros(test_count*tot_Prior,1);    
    
    
    for test_indx = 1:test_count
        test_indx
        for indx = 1:tot_Prior %for each prior density function
            indx
            % the one out answer that is being estimated
            answer_oneOut = interpol_exp_fT2(:,indx);
            prior = interpol_exp_fT2;
            prior(:,indx) = [];

            % this is valid a noise is known
            %n_sigma =  trapz(answer_oneOut)./SNR_lin;
            
            signal_power = sum((K2*answer_oneOut).^2) / N2;
            n_sigma = sqrt(signal_power/SNR);

            noise = n_sigma*normrnd(0, 1, [N2 ,1]); %assumes AWGN
            m = K2*answer_oneOut + noise; % generates simulated data



            % generate estimates of the bff
            [bff_est_bayes_sharp, compute_time_bayes_sharp]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior);
            [bff_est_bayes_tapered, compute_time_bayes_tapered]= bayes_estimator(bfv_tapered, m, K2, n_sigma, T2, tE, prior);        
            [bff_est_ilt, compute_time_ilt]= ilt_estimator(bfv_tapered, m, K2, n_sigma, T2, tE);              
            [bff_est_iltx, compute_time_iltx]= iltx_estimator(bfv_tapered, m, K2, n_sigma, T2, tE, tau2); 

            [bff_est_eht, compute_time_eht] = tapered_area(Tc, m, n_sigma, T2, tE, tau2);           




            bff_actual = bff_answer(answer_oneOut, bfv_sharp);


            bff_est_bayes_sharp_error((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp - bff_actual);
            bff_est_bayes_tapered_error((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_tapered - bff_actual); 
            bff_est_ilt_error((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_ilt - bff_actual); 
            bff_est_iltx_error((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_iltx - bff_actual); 
            bff_est_eht_error((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_eht - bff_actual); 
            
            bff_est_bayes_sharp_compute_time((test_indx-1)*(tot_Prior) + indx) = compute_time_bayes_sharp;
            bff_est_bayes_tapered_compute_time((test_indx-1)*(tot_Prior) + indx) = compute_time_bayes_tapered;
            bff_est_ilt_compute_time((test_indx-1)*(tot_Prior) + indx) = compute_time_ilt;
            bff_est_iltx_compute_time((test_indx-1)*(tot_Prior) + indx) = compute_time_iltx; 
            bff_est_eht_compute_time((test_indx-1)*(tot_Prior) + indx) = compute_time_eht;             
            
            
        end

    end

    %% Step 3 - Display results of the error made

    figure(2)
    clf
    hold on
    
    cdfplot(bff_est_bayes_sharp_error);
    cdfplot(bff_est_bayes_tapered_error);
    cdfplot(bff_est_ilt_error);
    cdfplot(bff_est_iltx_error);
    cdfplot(bff_est_eht_error);
    title('')
    legend('Bayes Sharp BFV', 'Bayes Tapered BFV', 'ILT Reproduced', ...
        'ILT+ Reproduced', 'EHT Reproduced','location','SouthEast')
    xlabel('BFF Absolute Error')
    ylabel('Prob(Error \textless  X)')
    xlim([0 0.4])
    hold off

    figure(3)
    clf
    subplot(2,2,1)
    hist(bff_est_bayes_sharp_error);
    legend('Bayes Sharp BFV')
    xlim([0 1])    
    subplot(2,2,2)    
    hist(bff_est_ilt_error);
    legend('ILT est.')
    xlim([0 1])
    subplot(2,2,3)     
    hist(bff_est_iltx_error);
    legend('ILT+ est')
    xlim([0 1])
    subplot(2,2,4)        
    hist(bff_est_eht_error);
    legend('EHT est.')
    xlim([0 1])


    

    figure(4)
    clf
    hold on
    
    cdfplot(bff_est_bayes_sharp_compute_time);
    
    
    lin = plot(bff_est_bayes_sharp_compute_time) ;

    
    cdfplot(bff_est_ilt_compute_time);
        delete(lin);
    cdfplot(bff_est_iltx_compute_time);
    cdfplot(bff_est_eht_compute_time);
    title('')
    legend('Bayes Method', 'ILT Reproduced', ...
        'ILT+ Reproduced', 'EHT Reproduced','location','SouthEast')
    xlabel('Compute Time (s)')
    ylabel('Prob(Time \textless  X)')
    hold off
    

    
    
    
    
    
    %root mean squared error
    bayes_sharp_RMSE = (mean((bff_est_bayes_sharp_error).^2))^.5
    bayes_tapered_RMSE = (mean((bff_est_bayes_tapered_error).^2))^.5
    ilt_RMSE = (mean((bff_est_ilt_error).^2))^.5
    iltx_RMSE = (mean((bff_est_iltx_error).^2))^.5
    eht_RMSE = (mean((bff_est_eht_error).^2))^.5

    RMSE_all(Tc_index,1) = bayes_sharp_RMSE;
    RMSE_all(Tc_index,2) = bayes_tapered_RMSE;
    RMSE_all(Tc_index,3) = ilt_RMSE;
    RMSE_all(Tc_index,4) = iltx_RMSE;
    RMSE_all(Tc_index,5) = eht_RMSE;
end

figure(4)
clf
hold on
plot(T2, RMSE_all(:,1:4))
% legend('Bayes Sharp BFV', 'Bayes Tapered BFV', 'ILT Reproduced', ...
%         'ILT+ Reproduced', 'EHT Reproduced','location','SouthEast')
xlabel('$T_c$')   
ylabel('BFF RMSE')
set(gca, 'XScale', 'log') 
xlim([2*tE 10])
ylim([0 0.2])

colormap flag

p1 = plot([33e-3 33e-3], [0 0.2], '--')
p1.Color = [0 0 0]
p2 = plot([90e-3 90e-3], [0 0.2], '--')
p2.Color = [0 0.5 0]


grid on
legend('Bayes Sharp BFV', 'Bayes Tapered BFV', 'ILT Reproduced', ...
        'ILT+ Reproduced','location','SouthEast')


%th=annotation('textarrow',[8e-3,33e-3],[0.12,0.12],'String','T_c = 33ms');

%te=annotation('textarrow',[8e-3,90e-3],[0.12,0.12],'String','T_c = 90ms');





