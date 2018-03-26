%% David Dobbie
% Recreating paper 1 (Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions)

%Aim: Extract the exponential time constants from a summation of several
%       type of them.

clc
clf
clear


time_constants = [1, 4];
exp_weightings = [1, 3];

time = 0:0.05:2;
init_individual_exp = [];



for t=time
    init_individual_exp = [init_individual_exp; exp_weightings.*exp(-time_constants*t)];
end

data= sum(init_individual_exp'); %adds expoential functions together


% generate the noise
noise_mean = 0;
noise_std_dev = 0.05;
noise = normrnd(noise_mean, noise_std_dev, 1,length(time));

data = noise + data;

figure(1)
plot(time, data);
title('Measured Signal (Sum of different exponentials)')
xlabel('Time [s]')
ylabel('Amplitude')


%% Step 1 Compression
% skipped at this point

K_1 = exp(-1./time)';

K_2 = ones(size(data'));

K_0 = kron(K_1, K_2);
%% Step 2 Optimisation

alpha = 0.01;
m = K_0;



c = ones(size(m));

alpha_hist = [];

figure(2)
hold on

f = c;

%this is the method that prevents it being divergent
for i=1:10
    k_square = K_0*K_0';    
    c = inv(k_square + alpha*eye(size(k_square)))*m;
    plot(c)
    
    alpha =  sqrt(size(c,2)) / norm(c); %implement eq 41
    
    alpha_hist = [alpha_hist; alpha];
    
    
end
hold off
%{
%this is the divergent method
for i=1:50   
    f = K_0'*c;
    
    c = (K_0 * f - m)/(-alpha);
    
    plot(time, c)
    
    alpha = length(c) / norm(c);
    
    alpha_hist = [alpha_hist; alpha];

    
end
%}
f = K_0.*c

var = K_0'*((K_0*K_0' + alpha*eye(size(k_square)))^-2)*K_0


figure(3)
plot(f)


figure(4)
plot(alpha_hist)
title('\alpha change over each computation')
xlabel('Iteration Number')
ylabel('\alpha')