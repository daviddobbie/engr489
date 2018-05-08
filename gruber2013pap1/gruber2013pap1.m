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


%close all
clc
clf
clear all
set(0,'defaultTextInterpreter','latex');
 


%% init variables
% number of data points in each dimension
N2 = 8000;
% number of bins in relaxation time grids
Ny = 1000;      
%sets how many singular values we compress to
sing_val=20; %no singular values
tE = 200e-6; % sample interval
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel

% RECREATE MODEL 4 from Gruber et al

%generate the density fn
T2_mean1 = 0.066
T2_var1 = 0.012

T2_mean2 = 0.2
T2_var2 = 0.022

% formation of distribution
f_answer = .25*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + .05*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';


porosity = trapz(f_answer);
f_answer = f_answer./porosity; % normalise to unity porosity

%delta distribut
%f_answer = zeros(Ny,1);
%f_answer(500) = 1;


Tc = 0.033; % estimated Tc of 33 ms (section 3.2)

figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f(T_2)$')
title('Correct Density Function of $T_2$');
hold on
%plot([Tc Tc],get(gca,'ylim'))

% p = ones(Ny,1); %% Full polarization
% p = 1 - 2*exp(-tau2 .* (1./T2) ); % inversion recovery
% p = 1 - exp(-tau2 * (1./T2) ); % saturation recovery



% generate the noise
noise_mean = 0;
n_std_dev = 0.2;

results_leng = 100;
results = zeros(1,results_leng);
results_var = zeros(1,results_leng);

for i = 1:results_leng

    noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);

    m = K2*f_answer + noise;

    %{
    figure(2)
    plot (tau2, m);
    xlabel('$t(s)$')
    ylabel('$M(t)$')
    title('Simulated Noisy Data M(t)');
    %}
    
  

    [A, A_sigma] = exponentialHaarTransform(tE, m, n_std_dev, Tc, T2, tau2);
    results(i) = A;  
    results_var(i) = A_sigma;
    
end    
   
figure(4)
h1 = histogram(results,'BinWidth', 0.05/7);
h1.BinLimits =([0.65 0.9]);
title('Estimated Tapered Area of Distribution')
xlabel('Tapered Area')
ylabel('Frequency')

taperedAreaMean_mean = mean(results)
taperedAreaMean_stddev = std(results)


taperedAreaVariance_mean = mean(results_var)

%% function definitions:

% Discretised Exponential Haar Transform. Used for computing the tapered
% areas
% INPUTS: 
%    m = N x 1 measurement vector
%    Tc = the bound and fluid fraction point (in time)
%    tE = sample interval of m
%    sigma_n = standard deviation of the noise
%    T2 = T2 relaxation axis
% OUTPUTS:
%    area = the T2 tapered area (where K(T2) > f(T2))
%    area_uncert = the T2 uncertainty of the tapered area

function [area area_uncert] = exponentialHaarTransform(tE, M, sigma_n, Tc, T2,tau2)
    area = 0;
    area_uncert = 0;

    % set up in table 1 of paper
    C = 0.7213 / Tc;
    alpha = 1.572*Tc;
    beta = 0.4087 / Tc;
    gamma = 1./T2 + beta;
    N = length(M);

    % implementing (-1)^n as vector operator
    t_tot = N * tE; % time axis
    nLen = ceil(t_tot / (2*alpha)); % the number of n values in length of measured data
    
    n_term_small = (-1).^(0:nLen-1);
    
    n_term_long = repelem(n_term_small,ceil(N/nLen));
    n_term = n_term_long(1:N);
    
    k = C.*n_term'.*exp(-beta * tau2);
    
    K = (C./gamma).*tanh(alpha*gamma);
    
    
    figure(7)
    hold on
    plot (T2, K);
    set(gca, 'XScale', 'log')
    xlabel('$T_2(s)$')
    ylabel('$K(T_2,T_C)$')
    title('The Exponential Haar Transform Kernel $K(T_2,T_c = 0.033s)$');
    hold off
    
    
    figure(8)
    hold on
    plot (tau2, k);
    xlabel('$t(s)$')
    ylabel('$K(t,T_c)$')
    title('The Exponential Haar Transform Kernel $k(t,T_c = 0.033s)$');
    hold off
    
    area = tE * (k' * M); % eq5
    area_uncert = (sigma_n)^2 *tE * ((tE * k')*k); % eq6
    
    %{
        figure(3)
        plot (tau2,k);
        xlabel('$t(s)$')
        ylabel('$M(t)$')
        title('n vect');
    %}

end

% Discretised Mellin transform. Assumes that m is discretised along tE
% sample interval (adapted from paper Mellin Transform Venk. 2010)
% INPUTS: 
%    m = N x 1 measurement vector
%    omega = the omega-th moment being calculated
%    tE = sample interval of m
%    poro = porosity (linear scaling function for mean)
% OUTPUTS:
%    the T2 moment
%    variance of T2 moment

function [moment var] = mellinTransform(m, omega, tE, poro, sigma_p, sigma_n);
       N = length(m);
       moment = 0;
       var = 0;

    if omega==0
        moment = 1;
        var = 0;
    elseif omega > 0
        tau_min = tE^omega; %eq 19a
        k = tau_min / gamma(omega+1); %eq 19a
        
        n = 2:1:N-1;
        
        % eq 19c-e
        delta_1 = (0.5 * tau_min * (2^omega - 1^omega) );
        delta_mid = 0.5 * tau_min * ((n+1).^omega - (n-1).^omega);
        delta_N = (0.5 * tau_min * (N^omega - (N-1)^omega) );
        delta = [delta_1 delta_mid delta_N];
        
        % note that erroneous values are apparent for a negative
        % measurement (leads to complex result from log)
        moment = k + 1/(gamma(omega + 1)*poro) * (delta*m); % eq18
        
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           omega
        end
        
        %eq 23
        var = ((delta.^2)*(delta.^2)'/gamma(omega+1)^2)*(sigma_n/poro)^2;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        return;
    elseif -1 < omega && omega < 0  %implement eq 22
        
        tau_min = tE^omega; %eq 19a
        k = tau_min / gamma(omega+1); %eq 19a
        
        n = 2:1:N-1;
        
        % eq 19c-e
        delta_1 = (0.5 * tau_min * (2^omega - 1^omega) );
        delta_mid = 0.5 * tau_min * ((n+1).^omega - (n-1).^omega);
        delta_N = (0.5 * tau_min * (N^omega - (N-1)^omega) );
        delta = [delta_1 delta_mid delta_N];
          
        %estimate 1st derivate of M at 0 with polyfit
        coeffc = polyfit(1:100, m(1:100)', 1);
        a1 = (coeffc(1))/(tE);
        m_est = coeffc(1).*(1:100) + coeffc(2);
        a1_stddev = std(abs(m(1:100) - m_est')); %standard deviation of slope
        
        const = ((a1*omega) / (omega+1)) * tau_min^((omega+1)/omega);
        
        moment = k + (1/(gamma(omega + 1)*poro)) * (const + delta*m);
        
        
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
        const + delta*m
        const
        delta*m
        omega
        end
        
        if moment < 5e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            %disp('delta * m')
            (delta*m);
            %disp('1/gamma(omg + 1)')
            1/(gamma(omega + 1)*poro);
            moment = NaN;
        end
        
        
        var = (((delta.^2)*(delta.^2)')/(gamma(omega+1))^2)*(sigma_n/poro)^2 + var;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        var = var + ((omega*tau_min^((omega+1)/omega))/gamma(omega + 2)) * (a1_stddev / poro)^2;
        return;
    else %allows outside case, shouldn't be called
        moment = 0;
        var = 0;
    end
end

