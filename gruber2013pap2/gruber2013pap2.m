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
sing_val=10; %no singular values
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

% generate the noise
noise_mean = 0;
n_std_dev = 0.2;

omega = linspace(-0.5,1,12);
T_cutoff = [0.01 0.1 1];




noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);

m = K2*f_answer + noise;  




% estimate tapered areas for different cut off times
tpdAreasVect = ones(size(T_cutoff,2), 2);
tpdAreaKern = ones(size(T_cutoff,2), Ny);

indx = 1;
for Tc = T_cutoff
    [tpd, tpd_var] = exponentialHaarTransform(tE, m, n_std_dev, Tc, T2,tau2);
    tpdAreasVect(indx,:) = [tpd, tpd_var];
    
    
    tpdAreaKern(indx,:) = exponetialHaarTransformKernel(Tc,T2);
    indx = indx + 1;
end

tpdAreasVect



% estimate different moments of the T2 distribution
momentVect = ones(size(omega,2), 2);
momentKern = ones(size(omega,2), Ny);


indx = 1;
for w = omega
    
    kern = T2.^w;
    momentKern(indx,:) = kern;
    
    
    [mom, mom_var] = mellinTransform(m, w, tE, porosity, 0, n_std_dev);
    momentVect(indx,:) = [mom sqrt(mom_var)]; %%since we use uncertainty
    indx = indx + 1;
end

momentVect

% create compressed m vector of values for optimisation
[m_comp, k_comp] = compressData(m,K2,10);

m_comp = m_comp'; % compressed m vector


M_opt = [m_comp; momentVect(:,1) ; tpdAreasVect(:,1)];


L_opt = [k_comp ; momentKern ; tpdAreaKern];


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
    n = floor(tau2/(2*alpha));
    n_term = ((-1).^n)';    
    k = C.*n_term'.*exp(-beta * tau2);
    area = tE * (k' * M); % eq5
    
    area_uncert = (sigma_n)^2 *tE * ((tE * k')*k); % eq6
    kernel = (C./gamma).*tanh(alpha*gamma);
end


function [kern] = exponetialHaarTransformKernel(Tc, T2)
    % set up in table 1 of paper
    C = 0.7213 / Tc;
    alpha = 1.572*Tc;
    beta = 0.4087 / Tc;
    gamma = 1./T2 + beta;
    kern = (C./gamma).*tanh(alpha*gamma);
end



% Calculate the actual tapered area of the distribution (fully polarised).
% uses the exponetial Haar transform.
% INPUTS: 
%    f = Ny x 1 desnity function vector
%    Tc = the bound and fluid fraction point (in time)
%    T2 = T2 relaxation axis
% OUTPUTS:
%    area = the T2 tapered area (where K(T2) > f(T2))
%  

function [area] = actualTaperedArea(f, Tc, T2)
    % set up in table 1 of paper
    C = 0.7213 / Tc;
    alpha = 1.572*Tc;
    beta = 0.4087 / Tc;
    gamma = 1./T2 + beta;   
    K = (C./gamma).*tanh(alpha*gamma);
    area = K * f;
end


% Compress the measured data to a set number of singular values. This
% involves extracting the largest singular values that make up the inherent
% behaviour of the system
% INPUTS: 
%    M = N x 1 measurement vector
%    K1 = input kernel
%    sing_val = number of singular values we compress to
% OUTPUTS:
%    M_compressed = the compressed measurement vector
%    k_compressed = the compressed kernel vector
% 

function [M_compressed K_compressed] = compressData(M,K,sing_val)
    N = length(M);

    %svd resultsin sorted by magnitude singular values. i.e we only have to
    %truncate to s1 rowsxcols.
    [U2, S2, V2] = svd(K);

    %only leave trunc number of largest values. This removes small weigthed
    %components that have little bearing on actual data.

    %Use if statements for truncation.
    if sing_val < N
        S2c = S2(1:sing_val,1:sing_val);
        U2c = U2(:,1:sing_val);
        V2c = V2(:,1:sing_val);
    else
        S2c = S2(1:sing_val,:);
        U2c = U2(:,1:sing_val);
        V2c = V2(:,:);

    end

    %set new compressed kernels
    K_compressed = S2c*V2c';
    M_compressed = (M'*U2c);
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
        var = (((delta.^2)*(delta.^2)')/(gamma(omega+1))^2)*(sigma_n/poro)^2 + var;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        var = var + ((omega*tau_min^((omega+1)/omega))/gamma(omega + 2)) * (a1_stddev / poro)^2;
        return;
    else %allows outside case, shouldn't be called
        moment = 0;
        var = 0;
    end
end


