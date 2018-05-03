%% David Dobbie
% Victoria University of Wellington
% Recreating paper 3
% Estimation of petrophysical and fluid properties using integral transforms in
% nuclear magnetic resonance)
%
% 
% F. K. Gruber et al / Journal of Magnetic Resonance 228 (2013) 104-115

%Aim: Using specific integral transforms to extract linear functionals that
%        describe the distrubtion function of T2 relaxation times without 
%        actually finding the distribution. This also includes aquiring the uncertainty 
%        of them for analysis.
%           Parts of the paper are
%           1) Estimate mu-th moment of distribution S2.1
%           2) Estimate tapered areas of the distribution S2.1
%           3) Extend to non-poolarised data S2.2
%           4) Other linear functionals to be estimated S2.3
%               a) moments of specified region  compute Mellin operator and Expo
%               Haar Transform
%               (using mixes of other transforms)


close all
clc
clf
clear all
set(0,'defaultTextInterpreter','latex');
 


%% init variables
% number of data points in each dimension
N2 = 5000;
% number of bins in relaxation time grids
Ny = 100;      
%sets how many singular values we compress to
sing_val=20; %no singular values
tE = 500*1e-6; % sample interval
T2 = logspace(log10(2*tE),1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel


%generate the density fn
T2_mean1 = 0.1
T2_var1 = 0.04

T2_mean2 = 0.006
T2_var2 = 0.09

% formation of distribution
f_answer = 1*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + 0.66*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';


porosity = trapz(f_answer);
f_answer = f_answer./porosity; % normalise to unity porosity

porosity = trapz(f_answer)
%delta distribut
%f_answer = zeros(Ny,1);
%f_answer(500) = 1;

figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f(T_2)$')
title('Correct Density Function of $T_2$');

p = ones(Ny,1); %% Full polarization
p = 1 - 2*exp(-tau2 .* (1./T2) ); % inversion recovery
p = 1 - exp(-tau2 * (1./T2) ); % saturation recovery


% generate the noise
noise_mean = 0.0;
n_std_dev = 0.02;
noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);

m = K2*f_answer + noise;






%% function definitions:

