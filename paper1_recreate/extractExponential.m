%% David Dobbie
% Recreating paper 1 (Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions)
%
% 1053-587X(02)03282-8

%Aim: Extract the exponential time constants from a summation of several
%       type of them.




clc
clf
clear

%% Paul Teal FLINT init variables
N1 = 50;       % number of data points in each dimension
N2 = 1000;

Nx = 100;      % number of bins in relaxation time grids
Ny = 101;      

tau1min = 1e-4;
tau1max = 10;
deltatau2 = 3.5e-4;

T1 = logspace(-2,1,Nx); %form T1 domain, use log since will be small
T2 = logspace(-2,1,Ny); %form T2 domain, use log since will be small

[T2a,T1a] = meshgrid(log10(T2),log10(T1));

%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*deltatau2;  
% measurement time array (note, no t=0 data available)
tau1 = logspace(log10(tau1min),log10(tau1max),N1)'; 

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data
K1 = 1-2*exp(-tau1 *(1./T1) );  % T1 relaxation data
%%


time_constants = [-0.2, -0.05];
exp_weightings = [1, 4];

time = tau2;
init_individual_exp = [];




for tc= time_constants
    init_individual_exp = [init_individual_exp; exp(time/tc)'];
    %%init_individual_exp = [init_individual_exp; normpdf(time, tc, 0.15) ];

end




data= sum(init_individual_exp); %adds exponetial functions together




%m = (K1.*data)'; %for pdf of the it

%K0 = kron(K1, K2);

K0 = K2;

% generate the noise
noise_mean = 0;
noise_std_dev = 0.05;
noise = normrnd(noise_mean, noise_std_dev, 1,length(time));
%m = data'
m = (noise + data)';

figure(1)
hold on
plot(time, m,'r');
%%plot (time,K2,'b');
plot(time,data, 'g');
hold off
legend("measured","noiseless data");
title('Measured Signal Function (Sum of different exponentials)')
xlabel('Time [s]')
ylabel('Amplitude')



%% Step 1 Compression
% skipped at this point


%% Step 2 Optimisation

alpha = 100;




c = ones(size(m));

alpha_hist = [];

figure(2)
clf
hold on
f = c;
title('c_r vector changing at each iteration');
%this is the method that prevents it being divergent
for i=1:20
    %k_square = K0*K0'; 
    stepFnMatrix = (heaviside(K0'*c)).* eye(Ny,Ny);
    k_square = K0 *  stepFnMatrix * K0';       %recreate eq 30
    %made symmetric and semi-positive definite
    c = inv(k_square + alpha*eye(size(k_square)))*m; %eq 29
    plot(time, c)

    alpha =  sqrt(size(c,2)) / norm(c); %implement eq 41
    
    alpha_hist = [alpha_hist; alpha];
    
    
end
hold off
%{
%this is the divergent method
for i=1:50   
    f = K0'*c;
    
    c = (K0 * f - m)/(-alpha);
    
    plot(time, c)
    
    alpha = length(c) / norm(c);
    
    alpha_hist = [alpha_hist; alpha];

    
end
%}
f = K0'*c;

var = K0'*((K0*K0' + alpha*eye(size(k_square)))^-2)*K0;


figure(3)
clf
hold on
set(gca, 'XScale', 'log')
plot(T2, 10.^(f)./2,'r')
hold off
legend("Density Function")
title('Probability Density Function for T_2')
ylabel('f(T_2)')
xlabel('T_2')

trapz(log10(T2), 10.^(f)./2)

figure(4)
plot(alpha_hist)
title('\alpha change over each computation')
xlabel('Iteration Number')
ylabel('\alpha')