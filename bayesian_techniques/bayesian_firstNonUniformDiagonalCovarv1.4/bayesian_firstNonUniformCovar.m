%% David Dobbie
% Victoria University of Wellington
% based on:
% Recreating "Bayesian NMR Relaxomtery" 29 Aug - 1 Sept 2017 
% Bayesian NMR Relaxometry, paper 6
% 
% Paul Teal

%Aim: Estimating quantities like fraction of bound fluids from NMR
%relaxomtery measurements, checking different priors

% this file implements 1-out cross validation

% algorithm goes as follows:
%   load up f mean prior from high quality data
%   transpose K matrix with weighting fo different variances of f
%   iterate through these different variances to find where rmse is
%       minmised
%   plot variance, bias, standard deviation


clc
clear

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)


% load high quality experimental data
loadedExperimental = load('testing_data\T2Distributions.mat');
T2_dist_cells = loadedExperimental.T2_distribution;
extracedT2Data = cell2mat(T2_dist_cells);
experimentalfT2 = extracedT2Data(:,2:2:end);
experimentalT2Axis = extracedT2Data(:,1); %assuming all have the same T2 axes spacing


figure(23)
clf
plot(experimentalT2Axis, experimentalfT2)
xlabel('time [s]')
ylabel('$f(T_2)$')
set(gca, 'XScale', 'log')



%loading model 2 from Gruber 2013 paper 2
density_funcload = load('model2.csv');
%density_funcload(:,2) = density_funcload(:,2)
[C,ia,ic]  = unique(density_funcload(:,1)),'stable';
density_funcload = density_funcload(ia,:);
%density_funcload = density_funcload(:~non_unique,:);
%{
figure(3)
clf
plot(density_funcload(:,1), density_funcload(:,2))
set(gca, 'XScale', 'log')
%}
%% Step 0: intialise variables

% number of data points in each dimension
N2 = 1000;
% number of bins in relaxation time grids
Ny = 30;      
%sets how many singular values we compress to
sing_val=5; %no singular values
tE = 200e-6;
%tE = 200e-6; % sample interval
T2 = logspace(-3,4,Ny); %form T2 domain, use log since will be small
%T2 = logspace(-5,1,Ny);
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel


%f_answer = density_funcload;
f_answer = interp1(density_funcload(:,1),density_funcload(:,2),T2,'pchip',0)';
f_answer = f_answer ./ trapz(f_answer);
porosity = trapz(f_answer);

%{
f_answer = zeros(Ny,1)
for idx = 1:Ny
    if T2(idx)>0.09
        f_answer(idx) = 1;
    end
end
f_answer = ones(Ny,1);
f_answer = 1*f_answer./sum(f_answer);
%}

figure(3)
clf
plot(T2, f_answer)
set(gca, 'XScale', 'log')



noise_mean = 0;
n_std_dev = 0.2;


%% make integral transforms

Tc = 90e-3;

% make integral transform for density
porosity_density_g_vector = ones(Ny, 1);

% make integral transfrom for bfv
bfv_density_g_vector = zeros(Ny ,1);
for idx = 1:Ny
    if T2(idx)<Tc
        bfv_density_g_vector(idx) = 1;
    end
end

transform1 = porosity_density_g_vector;
transform2 = bfv_density_g_vector;
transform3 = 1 - bfv_density_g_vector;



figure(1)
clf
hold on
plot(T2, transform1);
plot(T2, transform2);
plot(T2, transform3);

set(gca, 'XScale', 'log')
legend('g0', 'g1', 'g2')
hold off
grid on

f_mean_prior =mean(experimentalfT2')'; %mean of all previous experimental data
f_mean_prior_interpolated = interp1(experimentalT2Axis, f_mean_prior, T2, 'spline')';
f_mean_prior_interpolated = f_mean_prior_interpolated ./trapz(f_mean_prior_interpolated);

%f_answer = f_mean_prior_interpolated;

%f_mean_prior_interpolated = f_answer;
%f_mean_prior_interpolated = zeros(Ny,1);
%{
figure(8)
clf
hold on
%plot(experimentalT2Axis, f_mean_prior);
plot(T2, f_answer);
plot(T2, f_mean_prior_interpolated);
set(gca, 'XScale', 'log') 
hold off
%}

[tot_Prior  del] = size(experimentalfT2')
alph_length = 200;

transformResults1 = zeros(alph_length,2);
transformPredict1 = zeros(alph_length,1);
rmseTransform1 = zeros(1,alph_length);
transformResults2 = zeros(alph_length,2);
transformPredict2 = zeros(alph_length,1);
rmseTransform2 = zeros(1,alph_length);
transformResults3 = zeros(alph_length,2);
transformPredict3 = zeros(alph_length,1);
rmseTransform3 = zeros(1,alph_length);


alpha_estimations = zeros(tot_Prior,1);

% take one out for each of the set of answer density functions
for indx = 1:tot_Prior
    indx
    tic
    % the one out answer that is being estimated
    answer_oneOut = experimentalfT2(:,indx);
    prior = experimentalfT2;
    prior(:,indx) = [];
    
    n_std_dev = 0.2 .* trapz(answer_oneOut) %normalise for unknown distribution, SNR is known

    
    
    % the prior of other results that is being used
    f_uncertain_prior = std(prior')';
    f_uncertain_prior_interpol = interp1(experimentalT2Axis, f_uncertain_prior, T2, 'pchip',0)';
    %f_uncertain_prior_interpol = ones(Ny,1);
    
    f_mean_prior =mean(prior')'; %mean of all previous experimental data
    f_mean_prior_interpolated = interp1(experimentalT2Axis, f_mean_prior, T2, 'pchip',0)';
    f_answer = interp1(experimentalT2Axis, answer_oneOut, T2, 'pchip',0)';
    
    % plot each one-out versus prior
    figure(8)
    clf
    hold on
    %plot(experimentalT2Axis, f_mean_prior);
    plot(T2, f_answer);
    errorbar(T2, f_mean_prior_interpolated, f_uncertain_prior_interpol);
    
    set(gca, 'XScale', 'log') 
    hold off
       
    est_alpha = estimateAlphaFromPrior(n_std_dev, f_uncertain_prior_interpol);
    alpha_estimations(indx) = est_alpha;
    
    
    [alph transformResults1Individual, transformPredict1Individual, rmseTransform1Individual] ...
    = bayesianEstimateIntegralTransform(transform1, f_answer, ...
    n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated, f_uncertain_prior_interpol, alph_length);

    [alph transformResults2Individual, transformPredict2Individual, rmseTransform2Individual] ...
        = bayesianEstimateIntegralTransform(transform2, f_answer, ...
        n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated, f_uncertain_prior_interpol, alph_length);

    [alph transformResults3Individual, transformPredict3Individual, rmseTransform3Individual] ...
        = bayesianEstimateIntegralTransform(transform3, f_answer,...
        n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated,f_uncertain_prior_interpol, alph_length);

    
    transformResults1 = transformResults1 + transformResults1Individual;
    transformPredict1 = transformPredict1 + transformPredict1Individual;
    rmseTransform1 = rmseTransform1 + rmseTransform1Individual;
    
        
    transformResults2 = transformResults2 + transformResults2Individual;
    transformPredict2 = transformPredict2 + transformPredict2Individual;
    rmseTransform2 = rmseTransform2 + rmseTransform2Individual;
    
        
    transformResults3 = transformResults3 + transformResults3Individual;
    transformPredict3 = transformPredict3 + transformPredict3Individual;
    rmseTransform3 = rmseTransform3 + rmseTransform3Individual;
    toc
    
end


% normalise the average of all of the experiments

transformResults1 = transformResults1 ./ tot_Prior;
transformPredict1 = transformPredict1 ./ tot_Prior;
rmseTransform1 = rmseTransform1 ./ tot_Prior;

transformResults2 = transformResults2 ./ tot_Prior;
transformPredict2 = transformPredict2 ./ tot_Prior;
rmseTransform2 = rmseTransform2 ./ tot_Prior;


transformResults3 = transformResults3 ./ tot_Prior;
transformPredict3 = transformPredict3 ./ tot_Prior;
rmseTransform3 = rmseTransform3 ./ tot_Prior;

mean_est_alpha = mean(alpha_estimations);
std_est_alpha = std(alpha_estimations);

%{
[alph transformResults1, transformPredict1, rmseTransform1] ...
    = bayesianEstimateIntegralTransform(transform1, f_answer, ...
    n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated);

[alph transformResults2, transformPredict2, rmseTransform2] ...
    = bayesianEstimateIntegralTransform(transform2, f_answer, ...
    n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated);

[alph transformResults3, transformPredict3, rmseTransform3] ...
    = bayesianEstimateIntegralTransform(transform3, f_answer,...
    n_std_dev, noise_mean, T2, K2, tau2, f_mean_prior_interpolated);
%}



figure(2)
clf
set(gca, 'XScale', 'log')
hold on
plot(alph,transformResults1(:,1));
plot(alph,transformResults2(:,1));
plot(alph,transformResults3(:,1));




y1 = get(gca,'ylim');
l1 = line([mean_est_alpha mean_est_alpha],y1);
l1.Color = [0 .5 .5]
l1.LineStyle = '-';
l2 = line([mean_est_alpha-std_est_alpha mean_est_alpha-std_est_alpha],y1);
l2.Color = [0 .5 .5]
l2.LineStyle = '--';
l3 = line([mean_est_alpha+std_est_alpha mean_est_alpha+std_est_alpha],y1);
l3.Color = [0 .5 .5]
l3.LineStyle = '--';

hold off
xlabel('$\alpha$')
ylabel('$\langle I \rangle$')
%ylim([0 10])
grid on
legend('g0', 'g1', 'g2')



% plot variance of the estimator
figure(4)
clf
subplot(1,2,1)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on
plot(alph, transformPredict1)
plot(alph, transformPredict2)
plot(alph, transformPredict3)

y1 = get(gca,'ylim');
l1 = line([mean_est_alpha mean_est_alpha],y1);
l1.Color = [0 .5 .5]
l1.LineStyle = '-';
l2 = line([mean_est_alpha-std_est_alpha mean_est_alpha-std_est_alpha],y1);
l2.Color = [0 .5 .5]
l2.LineStyle = '--';
l3 = line([mean_est_alpha+std_est_alpha mean_est_alpha+std_est_alpha],y1);
l3.Color = [0 .5 .5]
l3.LineStyle = '--';

hold off
xlabel('$\alpha$')
ylabel('$\hat{\sigma_I}$ Computed')
grid on
ylim([10e-2 10e1])

subplot(1,2,2)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on
plot(alph, transformResults1(:,2))
plot(alph, transformResults2(:,2))
plot(alph, transformResults3(:,2))

y1 = get(gca,'ylim');
l1 = line([mean_est_alpha mean_est_alpha],y1);
l1.Color = [0 .5 .5]
l1.LineStyle = '-';
l2 = line([mean_est_alpha-std_est_alpha mean_est_alpha-std_est_alpha],y1);
l2.Color = [0 .5 .5]
l2.LineStyle = '--';
l3 = line([mean_est_alpha+std_est_alpha mean_est_alpha+std_est_alpha],y1);
l3.Color = [0 .5 .5]
l3.LineStyle = '--';


hold off

xlabel('$\alpha$')
ylabel('$\sigma_I$ Empirical' )
grid on
%ylim([10e-6 10e1])


figure(5)
clf
hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
plot(alph, rmseTransform1)
plot(alph, rmseTransform2)
plot(alph, rmseTransform3)


y1 = get(gca,'ylim');
l1 = line([mean_est_alpha mean_est_alpha],y1);
l1.Color = [0 .5 .5]
l1.LineStyle = '-';
l2 = line([mean_est_alpha-std_est_alpha mean_est_alpha-std_est_alpha],y1);
l2.Color = [0 .5 .5]
l2.LineStyle = '--';
l3 = line([mean_est_alpha+std_est_alpha mean_est_alpha+std_est_alpha],y1);
l3.Color = [0 .5 .5]
l3.LineStyle = '--';


hold off
xlabel('$\alpha$')
ylabel('rmse I Bayesian' )
grid on
%ylim([10e-3 10e2])





%% functions

% Aims to calculate an estimate of alpha only fro mthe prior data. The goal
% is to find an estimation of the variance of the density function that
% best minimises the rmse
% INPUTS:
%    n_std_dev: standard deviation of the measurement noise
%    covar_f_est: the estimated covariance of the density function
% OUTPUTS:
%    alpha_est: the alpha estimate
%
function [alpha_est] = estimateAlphaFromPrior(n_std_dev, covar_f_est)
    N =length(covar_f_est);
    alpha_est = n_std_dev/(norm(covar_f_est)/N)
end



% Calculates the two parameters of a normally distributed general integral
% transform given measurement data. Implements equation 19.
% INPUTS: 
%    g: the integral transform in the T2 domain
%    m: measured data vector
%    K: exponential kernel mapping from the T2 to time domain
%    Cf: covariance of density function f
%    Cn: covariance of measurement noise
%    mu_f: mean of density function
% OUTPUTS:
%    mean: mean of the integral transform estimation
%    std_dev: standard deviation of integral transform estimation
function [mean std_dev] = calcNormIntegralTransformGivenMeasured(g, m, K, Cf, Cn, mu_f, T2)

    R = Cf * K' * inv(K* Cf * K' + Cn);
    %{
    [diff, idx] = min(abs(T2-mu_f));
    mu_f_distrib = zeros(length(T2),1); mu_f_distrib(idx,1) = 1;
    %}
    mu_f_distrib = zeros(length(T2),1); %a zero prior
    mu_f_distrib = mu_f;
    %{
    figure(44)
    clf
    plot(T2, R * (m - K*mu_f_distrib) - mu_f_distrib)
    set(gca, 'XScale', 'log')
    ylim([0 0.1])
    grid on
    %}
    
    %mean = g' * (  R * (m - K*mu_f_distrib) + mu_f_distrib);
    mean = g' * (  R * (m - K*mu_f_distrib) + mu_f_distrib);
    std_dev = sqrt(g'  *   (Cf - R*K*Cf')  *  g);
    %std_dev = g'  *   inv(   inv(Cf) + K' *inv(Cn) * K  )  *  g;
end

% Performs the prediction of the resulting integral transform with multiple
% initialiasations of the noise. This is used the assess the performance of
% the algorithm
% INPUTS: 
%    g: 
% OUTPUTS:
function [alpha_axis, intTransform_givenalpha, intTransform_computed_uncertainty, rmse_bayesian]  = ...
    bayesianEstimateIntegralTransform(g_intTransform,  f_answer, n_std_dev, noise_mean, T2, K2, tau2, ...
    prior_f, covar_prior_f, alph_length)
    % test intergation over whole timeline window
    % test different regularisation coefficients, alpha.
    N2 = length(tau2);
    Ny = length(T2);
    
    alpha_length = alph_length;
    alpha_axis = logspace(-5,5,alpha_length);
    num_attempts = 20;

    intTransform_results = zeros(alpha_length,num_attempts);
    intTransform_computed_uncertainty = zeros(alpha_length,1);

    intTransform_givenalpha = zeros(alpha_length,2); %[mean, uncertainty]

    for alph_idx = 1:alpha_length
        alpha = alpha_axis(alph_idx);
        alph_idx;
        for idx_attempt = 1:num_attempts
            % init measured data
            noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
            m = K2*f_answer + noise;  

            Cf = (covar_prior_f + n_std_dev).*eye(Ny)./(alpha);
            Cn = (n_std_dev)^2*eye(N2);

            %mu_f = mean((T2'.*f_answer)');
            mu_f =exp((log(T2))*f_answer);

            % create reference delta density at certain frequency point
            [diff, idx] = min(abs(T2-mu_f));
            mu_f_distrib = zeros(Ny,1); mu_f_distrib(idx,1) = 1;
            
            %{
            figure(99)
            clf
            hold on
            plot(tau2,m)
            plot(tau2,K2*mu_f_distrib)        
            ylim([-0.5 1.5])
            %}
            
            [intTrans_mean, intTrans_uncertainty] = calcNormIntegralTransformGivenMeasured(g_intTransform, m,K2, Cf, Cn, prior_f, T2);
            actual_intTransform = (g_intTransform'*f_answer);
            
            %{
            figure(44)
            hold on
            plot(T2, f_answer, '-r')
            set(gca, 'XScale', 'log')        
            %}

            intTransform_results(alph_idx,idx_attempt) = intTrans_mean;
            intTransform_computed_uncertainty(alph_idx) = intTrans_uncertainty;
        end

    end

    intTransform_givenalpha(:,1) = mean(intTransform_results')';
    intTransform_givenalpha(:,2) = std(intTransform_results')';

    % plot rmse of the porosity estimator
    rmse_bayesian = sqrt(mean(((intTransform_results - actual_intTransform).^2)'));    
    
    %actual_intTransform
    
    %rmse_bayesian = rmse_bayesian ./ abs(mean(intTransform_results')); %normalise to nrmse
    
    
    

end





