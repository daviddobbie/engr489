%% David Dobbie
% Victoria University of Wellington
% Recreating paper 2 (Mellin Transform of CPMG data)
%
% 
% L. Venkataramanan et al/ journal of Magnetic Resonance 206 (2010) 20-31

%Aim: Computing moments of transverse T2 relaxation time from measured CPMG
%       data. Will implement in one dimension.


close all
clc
clf
clear all

set(1,'defaultTextInterpreter','latex');

%% init variables
% number of data points in each dimension
N2 = 5000;
% number of bins in relaxation time grids
Ny = 1000;      
%sets how many singular values we compress to
sing_val=20; %no singular values
tE = 500e-6; % sample interval
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data kernel


%recreation of MODEL D pg 26

%generate the density fn
T2_mean = 0.1;
T2_var = 0.1;
f_answer = normpdf(log10(T2), log10(T2_mean), T2_var)'/400;


figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f_{log_{10}}(T_2)({log_{10}}(T_2))$')
title('Correct Density Function of $T_2$');


% generate measured data from density function

n_stddev = 0.02;
n = normrnd(0,n_stddev,[N2,1]);

M = K2*f_answer + n;

figure(2)
plot (tau2, M);
xlabel('$t(s)$')
ylabel('$M(t)$')
title('Simulated Noisy Data M(t), $\sigma=0.2$');


% compression of the measurment data for computation
[U2, S2, V2] = svd(K2);

if sing_val < Ny
    S2c = S2(1:sing_val,1:sing_val);
    U2c = U2(:,1:sing_val);
    V2c = V2(:,1:sing_val);
else
    S2c = S2(1:sing_val,:);
    U2c = U2(:,1:sing_val);
    V2c = V2(:,:);
end

%set new compressed kernels
K2 = S2c*V2c';
M_comp = (U2c'*M);

figure(3)
plot(M_comp) %compressed M
xlabel('$Data points$')
ylabel('$M_{compressed}(t)$')
title('Compressed Simulated Noisy Data M(t), $\sigma=0.2$');



