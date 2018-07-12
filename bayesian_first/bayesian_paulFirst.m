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
tE = 10e-6;
%tE = 200e-6; % sample interval
T2 = logspace(log10(300e-6),log10(3),Ny); %form T2 domain, use log since will be small
%T2 = logspace(-5,1,Ny);
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel

f_answer = density_funcload;
f_answer = interp1(density_funcload(:,1),density_funcload(:,2),T2,'pchip')';
porosity = trapz(f_answer);

figure(3)
clf
plot(T2, f_answer)
set(gca, 'XScale', 'log')



noise_mean = 0;
n_std_dev = 0.2;

Tc = 0.033;

% make integral transform for density
porosity_density_g_vector = ones(Ny, 1);

% make integral transfrom for bfv
bfv_density_g_vector = zeros(Ny ,1);
for idx = 1:Ny
    if T2(idx)<Tc
        bfv_density_g_vector(idx) = 1;
    end
end


% init measured data
noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
m = K2*f_answer + noise;  

Cf = n_std_dev*eye(Ny);
Cn = n_std_dev* eye(N2);
mu_f = exp(log(T2)*f_answer) * ones(Ny,1)

[poro_mean poro_uncertainty] = calcNormIntegralTransformGivenMeasured(porosity_density_g_vector, m,K2, Cf, Cn, mu_f)
[bfv_mean bfv_uncertainty] = calcNormIntegralTransformGivenMeasured(bfv_density_g_vector, m,K2, Cf, Cn, mu_f)

est_bff = bfv_mean ./ poro_mean
actual_bff = (bfv_density_g_vector'*f_answer) / (porosity_density_g_vector'*f_answer)



%% functions

% Calculates the two parameters of a normally distributed general integral
% transform given measurement data. Implements equation 19.
% INPUTS: 
%    g: the integral transform in the T2 domain
%    m: measured data vector
%    K: exponential kernel mapping from the T2 to time domain
%    Cf: covariance of desnity function f
%    Cn: covariance of measurement noise
%    mu_f: mean of density function
% OUTPUTS:
%    mean: mean of the integral transform estimation
%    std_dev: standard deviation of integral transform estimation
function [mean std_dev] = calcNormIntegralTransformGivenMeasured(g, m, K, Cf, Cn, mu_f)

    R = Cf * K' * inv(K* Cf * K' + Cn);
    
    mean = g' * (R * (m - K*mu_f) - mu_f);
    std_dev = g'  *   inv(   inv(Cf) + K' *inv(Cn) * K  )  *  g;
end





