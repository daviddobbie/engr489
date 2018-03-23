%% David Dobbie
% Recreating paper 1 (Solving Fredholm Integrals of the first Kind With
% Tensor Product Strcuture in 2 and 2.5 Dimensions)

%Aim: Extract the exponential time constants from a summation of several
%       type of them.

clc;
clf;
clear;


time_constants = [1, 4, 8];
exp_weightings = [1, 3, 5];

time = 0:0.01:2;
init_individual_exp = [];



for t=time
    init_individual_exp = [init_individual_exp; exp_weightings.*exp(-time_constants*t)];
end

data= sum(init_individual_exp'); %adds expoential functions together


% generate the noise
noise_mean = 0;
noise_std_dev = 0.05;
noise = normrnd(noise_mean, noise_std_dev, 1,length(time));

data = noise + data

figure(1)
plot(time, data);
title('Measured Signal (Sum of different exponentials)')
xlabel('Time [s]')
ylabel('Amplitude')