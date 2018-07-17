%% David Dobbie
% Victoria University of Wellington
% Recreating "Bayesian NMR Relaxomtery" 29 Aug - 1 Sept 2017 
% Bayesian NMR Relaxometry, paper 6
% 
% Paul Teal

%Aim: Estimating quantities like fraction of bound fluids from NMR
%relaxomtery measurements

% algorithm goes as follows:
% 


clc
clf
clear

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)


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
tE = 500e-6;
%tE = 200e-6; % sample interval
T2 = logspace(-3,0.3,Ny); %form T2 domain, use log since will be small
%T2 = logspace(-5,1,Ny);
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel

%{
f_answer = density_funcload;
f_answer = interp1(density_funcload(:,1),density_funcload(:,2),T2,'pchip')';
porosity = trapz(f_answer);
%}

f_answer = zeros(Ny,1)
for idx = 1:Ny
    if T2(idx)>0.09
        f_answer(idx) = 1;
    end
end
f_answer = 1*f_answer./sum(f_answer);

figure(3)
clf
plot(T2, f_answer)
set(gca, 'XScale', 'log')



noise_mean = 0;
n_std_dev = 0.02;


%% make integral transforms

Tc = 33e-3;

% make integral transform for density
porosity_density_g_vector = ones(Ny, 1);

% make integral transfrom for bfv
bfv_density_g_vector = zeros(Ny ,1);
for idx = 1:Ny
    if T2(idx)<Tc
        bfv_density_g_vector(idx) = 1;
    end
end

% test intergation over whole timeline window
% test different regularisation coefficients, alpha.
alpha_length = 50;
alpha_axis = logspace(-5,5,alpha_length);
num_attempts = 50;

porosity_results = zeros(alpha_length,num_attempts);
porosity_computed_uncertainty = zeros(alpha_length,1);

porosity_givenalpha = zeros(alpha_length,2) %[mean, uncertainty]

for alph_idx = 1:alpha_length
    alpha = alpha_axis(alph_idx);
    alph_idx
    for idx_attempt = 1:num_attempts
        % init measured data
        noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
        m = K2*f_answer + noise;  

        Cf = (n_std_dev)^2*eye(Ny)./(alpha);
        Cn = (n_std_dev)^2*eye(N2);
        
        %mu_f = mean((T2'.*f_answer)');
        mu_f =exp((log(T2))*f_answer);
        
        % create reference delta density at certain frequuency point
        [diff, idx] = min(abs(T2-mu_f));
        mu_f_distrib = zeros(Ny,1); mu_f_distrib(idx,1) = 1;

        figure(99)
        clf
        hold on
        plot(tau2,m)
        plot(tau2,K2*mu_f_distrib)        
        ylim([-0.5 1.5])

        [poro_mean poro_uncertainty] = calcNormIntegralTransformGivenMeasured(porosity_density_g_vector, m,K2, Cf, Cn, mu_f, T2);
        actual_porosity = (porosity_density_g_vector'*f_answer);
        
        figure(44)
        hold on
        plot(T2, f_answer, '-r')
        set(gca, 'XScale', 'log')        
        
        
        porosity_results(alph_idx,idx_attempt) = poro_mean;
        porosity_computed_uncertainty(alph_idx) = poro_uncertainty;
    end
    
end

porosity_givenalpha(:,1) = mean(porosity_results')';
porosity_givenalpha(:,2) = std(porosity_results')';

figure(2)
clf
plot(alpha_axis,porosity_givenalpha(:,1));
set(gca, 'XScale', 'log')
xlabel('$\alpha$')
ylabel('$\langle \phi \rangle$')
ylim([0 1.2])
grid on

% plot variance of the estimator
figure(4)
clf
subplot(1,2,1)
plot(alpha_axis, porosity_computed_uncertainty)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('$\alpha$')
ylabel('$\hat{\sigma_\phi}$ Computed')
grid on
ylim([10e-6 10e0])

subplot(1,2,2)
plot(alpha_axis, porosity_givenalpha(:,2))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('$\alpha$')
ylabel('$\sigma_\phi$ Empirical' )
grid on
ylim([10e-6 10e0])

% plot rmse of the porosity estimator

rmse_bayesian = sqrt(mean(((porosity_results - actual_porosity).^2)'));


figure(5)
clf

plot(alpha_axis, rmse_bayesian)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('$\alpha$')
ylabel('RMSE Porosity Bayesian' )
grid on
ylim([10e-4 10e0])


%{
[bfv_mean bfv_uncertainty] = calcNormIntegralTransformGivenMeasured(bfv_density_g_vector, m,K2, Cf, Cn, mu_f)
actual_bfv = (bfv_density_g_vector'*f_answer)

est_bff = bfv_mean ./ poro_mean
actual_bff = (bfv_density_g_vector'*f_answer) / (porosity_density_g_vector'*f_answer)
%}

%% functions

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
    
    figure(44)
    clf
    %plot(T2, R * (m - K*mu_f_distrib)+ mu_f_distrib)
    plot(T2, R * (m - K*mu_f_distrib) - mu_f_distrib)
    set(gca, 'XScale', 'log')
    ylim([0 0.1])
    grid on
    %ylim([-2 2])
    
    
    %mean = g' * (  R * (m - K*mu_f_distrib) + mu_f_distrib);
    mean = g' * (  R * (m - K*mu_f_distrib) - mu_f_distrib);
    std_dev = sqrt(g'  *   (Cf - R*K*Cf')  *  g);
    %std_dev = g'  *   inv(   inv(Cf) + K' *inv(Cn) * K  )  *  g;
end





