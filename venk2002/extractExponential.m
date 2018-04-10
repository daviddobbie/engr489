%% David Dobbie
% Recreating paper 1 (Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions)
%
% 1053-587X(02)03282-8
% Venkataramanan et al Solving Fredholm Integrals of the First Kind . 2002

%Aim: Extract the exponential time constants from a summation of several
%       type of them.



close all
clc
clf
clear all

%% Paul Teal FLINT init variables
N1 = 50;       % number of data points in each dimension
N2 = 1000;

Nx = 100;      % number of bins in relaxation time grids
Ny = 101;      

%sets how we compress the data points.
% NOTE that:
%    N2 > trunc2 >= Ny
%    N1 > trunc1 >= Nx
trunc1=15;
trunc2=25;

tau1min = 1e-4;
tau1max = 10;
deltatau2 = 3.5e-4;

T1 = logspace(-3,1,Nx); %form T1 domain, use log since will be small
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small

[T2a,T1a] = meshgrid(log10(T2),log10(T1));

%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*deltatau2;  
% measurement time array (note, no t=0 data available)
tau1 = logspace(log10(tau1min),log10(tau1max),N1)'; 

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data
K1 = 1-2*exp(-tau1 *(1./T1) );  % T1 relaxation data
%%


time_constants_y = [-0.3, -0.05];
time_constants_x = [-0.3, -0.3];
exp_weightings_y = [1, 4];

time = tau2;
initexpY = [];

timeX = tau1;
initexpX = [];

for tc= time_constants_y
    initexpY = [initexpY; exp(time/tc)'];
    initexpX = [initexpX; exp(timeX/tc)'];
    %%init_individual_exp = [init_individual_exp; normpdf(time, tc, 0.15) ];

end



initexpY= sum(initexpY); %adds exponetial functions together
initexpX = sum(initexpX);

data = initexpX'*initexpY;

%m = (K1.*data)'; %for pdf of the it

%K0 = kron(K1, K2);



% generate the noise
noise_mean = 0;
noise_std_dev = 0.1;
noise = normrnd(noise_mean, noise_std_dev, length(timeX),length(time));
%m = data'
m = (noise + data)';

figure(1)
hold on
mesh(m);
%%plot (time,K2,'b');
plot(time,data, 'g');
hold off
title('Measured Signal Function (Sum of different exponentials)')
xlabel('Time 1  \tau_1 [s]')
ylabel('Time 2  \tau_2 [s]')
zlabel('Amplitude')



%% Step 1 Compression
% implement TSVD (truncated singular value decomposition), get the desired
% rank for only the largest singular values in the system.


%svd resultsin sorted by magnitude singular values. i.e we only have to
%truncate to s1 rowsxcols.
[U1, S1, V1] = svd(K1);
[U2, S2, V2] = svd(K2);




%only leave trunc number of largest values. This removes small weigthed
%components that have little bearing on actual data.

%Use if statements for truncation.
if trunc1 < Nx
    S1c = S1(1:trunc1,1:trunc1);
    U1c = U1(:,1:trunc1);
    V1c = V1(:,1:trunc1);
else
    S1c = S1(1:trunc1,:);    
    U1c = U1(:,1:trunc1);
    V1c = V1(:,:);
end

if trunc2 < Ny
    S2c = S2(1:trunc2,1:trunc2);
    U2c = U2(:,1:trunc2);
    V2c = V2(:,1:trunc2);
else
    S2c = S2(1:trunc2,:);
    U2c = U2(:,1:trunc2);
    V2c = V2(:,:);
    
end







%set new compressed kernels
K1 = S1c*V1c';
K2 = S2c*V2c';


K0 = kron(K1,K2);

m = (U1c'*m'*U2c);


%% Step 2 Optimisation

alpha = 100;




c = ones([trunc1*trunc2  1]);

% lexiographical sorting of the m matrix
mlex = sortrows(reshape(m,1,[]).');

alpha_hist = [];

figure(2)
clf
hold on
f = c;
title('c_r vector changing at each iteration');


% note that at more compressed values, we get more chance of divergent
% alpha, this needs to be considered. Takes more effort!

%this is the method that prevents it being divergent
for i=1:20
    %k_square = K0*K0'; 
    stepFnMatrix = (heaviside(K0'*c))'.* eye(Nx*Ny);
    k_square = K0 *stepFnMatrix * K0';       %recreate eq 30
    %made symmetric and semi-positive definite
    
    % =========> ERROR HERE
    c = inv(k_square + alpha*eye(trunc1*trunc2))*mlex; %eq 29
    % ===========
    
    plot(c)

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

f = reshape(f,Nx,Ny); %back to 2D representation
%var = K0'*((K0*K0' + alpha*eye(size(k_square)))^-2)*K0;


figure(3)
clf
hold on
set(gca, 'XScale', 'log')
mesh(f)
hold off
legend("Density Function")
title('Probability Density Function for T_2')
ylabel('f(T_2)')
xlabel('T_2')

%trapz(log10(T2), 10.^(f)./2)

figure(4)
plot(alpha_hist)
title('\alpha change over each computation')
xlabel('Iteration Number')
ylabel('\alpha')