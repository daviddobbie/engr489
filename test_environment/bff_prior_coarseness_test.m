%% Test Environment for BBF Estimation from T2 Relaxation 
% Exmaining the coarseness of the prior distribution
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

    %% Step 2: Begin Leave One Out Cross Validation

    test_count = 50;

    bff_est_iltx_results = zeros(test_count*tot_Prior,1);    
    bff_est_bayes_sharp_error_29prior = zeros(test_count*tot_Prior,1);
    bff_est_bayes_sharp_error_20prior = zeros(test_count*tot_Prior,1);
    bff_est_bayes_sharp_error_10prior = zeros(test_count*tot_Prior,1);
    bff_est_bayes_sharp_error_5prior = zeros(test_count*tot_Prior,1);
    bff_est_bayes_sharp_error_1prior = zeros(test_count*tot_Prior,1);
    
    
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

            random_indices= randperm(tot_Prior -1);

            indx_rand_prior_20 = random_indices(1:20);
            indx_rand_prior_10 = random_indices(1:10);
            indx_rand_prior_5 = random_indices(1:5);
            indx_rand_prior_1 = random_indices(1);
            
            
            
            
            prior30 = prior;
            prior20 = prior(:, indx_rand_prior_20);
            prior10 = prior(:, indx_rand_prior_10);
            prior5 = prior(:, indx_rand_prior_5);            
            prior1 = prior(:, indx_rand_prior_1);
            
            

            % generate estimates of the bff
            [bff_est_bayes_sharp_29]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior30);       
            [bff_est_bayes_sharp_20]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior20);       
            [bff_est_bayes_sharp_10]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior10);       
            [bff_est_bayes_sharp_5]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior5);       
            [bff_est_bayes_sharp_1]= bayes_estimator(bfv_sharp, m, K2, n_sigma, T2, tE, prior1);       
            bff_est_iltx = iltx_estimator(bfv_tapered, m, K2, n_sigma, T2, tE, tau2); 




            
            bff_actual = bff_answer(answer_oneOut, bfv_sharp);

            bff_est_iltx_results((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_iltx - bff_actual); 
            bff_est_bayes_sharp_error_29prior((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp_29 - bff_actual); 
            bff_est_bayes_sharp_error_20prior((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp_20 - bff_actual); 
            bff_est_bayes_sharp_error_10prior((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp_10 - bff_actual); 
            bff_est_bayes_sharp_error_5prior((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp_5 - bff_actual); 
            bff_est_bayes_sharp_error_1prior((test_indx-1)*(tot_Prior) + indx) = abs(bff_est_bayes_sharp_1 - bff_actual); 

            
            
        end

    end

    %% Step 3 - Display results of the error made

    figure(2)
    clf
    hold on
    
    cdfplot(bff_est_bayes_sharp_error_29prior);
    cdfplot(bff_est_bayes_sharp_error_20prior);    
    cdfplot(bff_est_bayes_sharp_error_10prior);      
    cdfplot(bff_est_bayes_sharp_error_5prior);      
    cdfplot(bff_est_bayes_sharp_error_1prior);      
    
    cdfplot(bff_est_iltx_results);  
    
    title('')
    xlabel('BFF Absolute Error')
    ylabel('Prob(Error \textless  X)')
    xlim([0 0.4])
    hold off

    legend('Prior of 29','Prior of 20','Prior of 10','Prior of 5', ...
        'Prior of 1', 'ILT+ Reproduced')
    
    
    
    
    
    %root mean squared error
    bayes_sharp_RMSE = (mean((bff_est_bayes_sharp_error).^2))^.5
    bayes_tapered_RMSE = (mean((bff_est_bayes_tapered_error).^2))^.5

    RMSE_all(Tc_index,1) = bayes_sharp_RMSE;
    RMSE_all(Tc_index,2) = bayes_tapered_RMSE;
end

figure(4)
clf
hold on
plot(T2, RMSE_all)
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





