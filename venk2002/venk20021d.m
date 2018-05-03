%% David Dobbie
% Victoria University of Wellington
% Recreating paper 1 (Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions)
%
% 1053-587X(02)03282-8
% Venkataramanan et al Solving Fredholm Integrals of the First Kind . 2002

%Aim: Extract the exponential time constants from a summation of several
%       type of them.




clc
clf
clear all
set(0,'defaultTextInterpreter','latex');

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

T1 = logspace(log10(deltatau2),log10(tau1max),Nx); %form T1 domain, use log since will be small
T2 = logspace(log10(deltatau2),log10(tau1max),Ny); %form T2 domain, use log since will be small

%[T2a,T1a] = meshgrid(log10(T2),log10(T1));

%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*deltatau2;  
% measurement time array (note, no t=0 data available)
tau1 = logspace(log10(tau1min),log10(tau1max),N1)'; 

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data
K1 = 1-2*exp(-tau1 *(1./T1) );  % T1 relaxation data
%%


%generate the density fn
T2_mean1 = 0.1
T2_var1 = 0.004

T2_mean2 = 0.05
T2_var2 = 0.005


f_answer = 2*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + 0.5*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';


f_answer = f_answer./trapz(f_answer); % normalise to unity porosity

figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2 [s]$')
ylabel('$f(T_2)$')
title('Correct Density Function of $T_2$');



%m = (K1.*data)'; %for pdf of the it

%K0 = kron(K1, K2);



% generate the noise
noise_mean = 0.0;
n_std_dev = 0.02;
noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
%m = data'
m = K2*f_answer + noise;

m_orig = m;

figure(2)
hold on
%%plot (time,K2,'b');
plot(tau2,m);
hold off
title('Measured Signal Function (Sum of different exponentials)')
xlabel('Time 2  $ \tau_2 $ [s]')
ylabel('M(t)')



%% Step 1 Compression
% implement TSVD (truncated singular value decomposition), get the desired
% rank for only the largest singular values in the system.

%{
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

m = (m'*U2c);

figure(3)
hold on
%%plot (time,K2,'b');
plot(1:trunc2,m);
hold off
title('Compressed Measured Signal')
xlabel('Time 2  $ \tau_2 $ [s]')
ylabel('$M_c(t)$')

%}

[U2, S2, V2] = svd(m);
S2(trunc2:N2) = 0;
m = U2*S2*V2';





%% Step 2 Optimisation

alpha = 1000;




c = ones([length(m)  1]);

% lexiographical sorting of the m matrix
mlex = sortrows(reshape(m,1,[]).');
m = mlex;
%K2 = sortrows(reshape(K2,Ny,N2).');

alpha_hist = [];

figure(4)
clf
hold on
f = c;
title('$c_r$ vector changing at each iteration');


% note that at more compressed values, we get more chance of divergent
% alpha, this needs to be considered. Takes more effort!

%this is the method that prevents it being divergent
for i=1:20
    %k_square = K2*K2'; 
    stepFnMatrix = (heaviside(K2'*c))'.* eye(Ny);
    k_square = K2 *stepFnMatrix * K2';       %recreate eq 30
    %made symmetric and semi-positive definite
    
    c = inv(k_square + alpha*eye(length(m)))*m; %eq 29
    
    plot(c)

    alpha = sqrt(size(c,2))/ norm(c); %implement eq 41
    
    %alpha =1000;
    
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
f = K2'*c;
%f = f;
f = f ./ trapz(f); %normalise to unity
%f = reshape(f,Nx,Ny); %back to 2D representation
%var = K0'*((K0*K0' + alpha*eye(size(k_square)))^-2)*K0;







figure(5)
plot(alpha_hist)
title('$\alpha$ change over each computation')
xlabel('Iteration Number')
ylabel('$\alpha$')

figure(6)
clf
hold on
set(gca, 'XScale', 'log')
plot(T2,f);
plot(T2,f_answer);
hold off
lgd = legend("$\hat{f}(t)$", "f(t)");
set(lgd,'Interpreter','latex');
title('Probability Density Function for $T_2$')
ylabel('$f(T_2)$')
xlabel('$T_2 [s]$')


figure(7)
clf
hold on
plot(tau2,K2*f);
plot(tau2,m_orig);
hold off
lgd = legend("$ \hat{M}(t) $", "M(t)");
set(lgd,'Interpreter','latex');
title('$ \hat{M}(t) $ given acquired density $\hat{f}(T_2)$')
ylabel('$ \hat{M}(t) $')
xlabel('$t [s]$')
