%% David Dobbie
% Victoria University of Wellington
% Recreating paper 3
% Estimation of petrophysical and fluid properties using integral transforms in
% nuclear meagnetic resonance)
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

%generate the density fn
T2_mean = 0.3
T2_var = 0.04

% normal dist answer
f_answer = normpdf(log10(T2), log10(T2_mean), sqrt(T2_var))';

%{
%recration MODEL C pg 26
T2_mean = 0.018
T2_var = 0.43

% normal dist answer
f_answer = normpdf(T2, log10(T2_mean), sqrt(T2_var))';

%}

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

n_stddev = 0.2


%% function definitions:

