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
set(0,'defaultTextInterpreter','latex');
 


%% init variables
% number of data points in each dimension
N2 = 5000;
% number of bins in relaxation time grids
Ny = 1000;      
%sets how many singular values we compress to
sing_val=20; %no singular values
tE = 500*1e-6; % sample interval
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

%recreation of MODEL D pg 26

%generate the density fn
T2_mean = 0.3
T2_var = 0.2

% normal dist answer
f_answer = normpdf(log10(T2), log10(T2_mean), T2_var)'/1700;

%delta distribut
%f_answer = zeros(Ny,1);
%f_answer(500) = 1;

figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f_{log_{10}}(T_2)({log_{10}}(T_2))$')
title('Correct Density Function of $T_2$');


% generate measured data from density function

n_stddev = 0.2;
porosity = trapz(f_answer);
n_stddev = n_stddev * porosity

N_histleng = 50;

hist_data = zeros(N_histleng,2);
for hist_indx = 1:N_histleng
    [m, v] = meanAndVarEstimation(N2, Ny, sing_val, tE, f_answer, n_stddev);
    hist_data(hist_indx,1) = m;
    hist_data(hist_indx,2) = v;
end
figure(5)
h1 = histogram(hist_data(:,1), 25)

estmean_avg = mean(hist_data(:,1))
estmean_med = median(hist_data(:,1))

xlabel('$T_{2,LM}$ Mellin Transfom (sec)')
ylabel('Frequency')
title('Mellin Transform Estimated  $\langle T_2 \rangle$ Relaxation Times')

figure(6)
h2 = histogram(hist_data(:,2), 25)

estvar_avg = mean(hist_data(:,2))
estvar_med = median(hist_data(:,2))

xlabel('$\sigma^{2}_{log_{10} T_2}$ Mellin Transfom (sec)')
ylabel('Frequency')
title('Mellin Transform Estimated $\sigma^{2}_{log_{10} \langle T_2 \rangle}$ variance Relaxation Times')
%% function definitions:

% Discretised Mellin transform. Assumes that m is discretised along tE
% sample interval
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
        
%         hold on
%         plot (delta)
%         hold off
        
        %delta = [(0.5*tau_min(2.^omega-1.^omega)) delta (0.5*tau_min*(N.^omega-(N-1).^omega))];
        
        %omega
        
        % note that erroneous values are apparent for a negative
        % measurement (leads to complex result from log)
        moment = k + 1/(gamma(omega + 1)*poro) * (delta*m); % eq18
        
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           omega
        end
        
        %{
        if moment < 1e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            %disp('delta * m');
            (delta*m);
            %disp('1/gamma(omg + 1)');
            1/(gamma(omega + 1)*poro);
            moment = NaN;
        end
           %}
        
        %eq 23
        var = ((delta.^2)*(delta.^2)'/gamma(omega+1)^2)*(sigma_n/poro)^2;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        return;
    elseif -1 < omega && omega < 0  %implement eq 22
        
        tau_min = tE^omega; %eq 19a
        k = tau_min/gamma(omega+1); %eq 19a 
        
        n = 2:1:N-1;
        delta_1 = (0.5 * tau_min * (2^omega - 1^omega) );
        delta_mid = 0.5 * tau_min * ((n+1).^omega - (n-1).^omega);
        delta_N = (0.5 * tau_min * (N^omega - (N-1)^omega) );
        delta = [delta_1 delta_mid delta_N];
        %{
            hold on
            plot (delta)
            hold off
        %} 
        
        %estimate 1st derivate of M at 0 with polyfit
        coeffc = polyfit(1:200, m(1:200)', 1);
        a1 = (coeffc(1))/(tE);
        %a1 = (m(2) - m(1))/(tE); % gradient of m at t=0 approximation
        
        a1_stddev = 0; %standard deviation of slope
        
        const = ((a1*omega) / (omega+1)) * tau_min^((omega+1)/omega);
        
        moment = k + (1/(gamma(omega + 1)*poro)) * (const + delta*m);
        
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           omega
        end
        
        
        if moment < 1e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            %disp('delta * m')
            (delta*m);
            %disp('1/gamma(omg + 1)')
            1/(gamma(omega + 1)*poro);
            moment = NaN;
        end
        
        
        var = (((delta.^2)*(delta.^2)')/(gamma(omega+1))^2)*(sigma_n/poro)^2;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        var = var + ((omega*tau_min^((omega+1)/omega))/gamma(omega + 2)) * (a1_stddev / poro)^2;
        return;
    else %allows outside case, shouldn't be called
        moment = 0;
        var = 0;
    end
end




function [momentT2, varT2] = meanAndVarEstimation(N2, Ny, sing_val, tE, f_answer, n_stddev)
DEBUG = 0;
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data kernel

n = n_stddev*normrnd(0,1,[N2,1]);

M = K2*f_answer + n;
%M = M .* heaviside(M);
if DEBUG
    figure(2)
    plot (tau2, M);
    xlabel('$t(s)$')
    ylabel('$M(t)$')
    title('Simulated Noisy Data M(t)');
end

% compression of the measurment data for computation
[U2, S2, V2] = svd(M);

S2c = zeros(size(S2));
if sing_val < Ny
    S2c(1:sing_val,1) = S2(1:sing_val,1);
else
    S2c(1:sing_val,1) = S2(1:sing_val,1);

end

%set new compressed kernels
M_comp = (U2*S2c*V2');

if DEBUG
    figure(3)
    plot(M_comp) %compressed M
    xlabel('Data points')
    ylabel('$M_{compressed}(t)$')
    title('Compressed Simulated Noisy Data M(t), $\sigma_{\epsilon}=0.2$');
end

porosity_answer = trapz(f_answer);

%mom = mellinTransform(M_comp, 1, tE, porosity_answer, 0.001 ,n_stddev);

% use G(omega) = ln<(T_2)^omega>, plot it

omega_axis = linspace(-0.5,1,30)';

result_axis = zeros(length(omega_axis), 2);
J = zeros(1, length(omega_axis));
Sigma_G = eye(length(omega_axis));
indx_covar_mat = 1;

% printing the index for delta vector
%{
figure(66)
xlabel("index - i")
ylabel("$\Delta_i$ values")
%}

for inx = 1:length(omega_axis)
    omg = omega_axis(inx);
    % certain about porosity
    [mom var] = mellinTransform(M_comp, omg, tE , porosity_answer, 0.001 ,n_stddev);
    J(indx_covar_mat) = 1/mom;
    
    Sigma_G(indx_covar_mat,indx_covar_mat) = var;
    
    result_axis(indx_covar_mat,:) = [log(mom)  var];
    
    indx_covar_mat = indx_covar_mat +1; 
    
end 

G_covar = J * Sigma_G * J';

[omega_axis result_axis(:,1)];

r_a = result_axis(:,1);

% clean up invalid values
n1_r_a = r_a(~isnan(r_a));
n2_r_a = n1_r_a(~isinf(n1_r_a));

% clean up invalid values
n1_omega_axis = omega_axis(~isnan(r_a));
n2_omega_axis = n1_omega_axis(~isinf(n1_r_a));

if DEBUG
    figure(4)
    plot(n2_omega_axis, n2_r_a(:,1))
    xlabel("$\omega^{th}$ moment")
    ylabel("$G(\omega)$ Mellin transform of data")
    title("$ G(\omega) \equiv ln \langle T_2^\omega \rangle $")
end

%quadratic fit on Mellin transform (2nd order polynomial)
coeffc = polyfit(n2_omega_axis, n2_r_a, 2);
var_lnT2 = abs(coeffc(1)/2);
moment_lnT2 = coeffc(2);

momentT2 = exp(moment_lnT2);
varT2 = momentT2 * var_lnT2;
return;
end




% Estimates the porosity and its uncertainty - correlates to area under T2
% dist.
function [poro, poro_uncert] = calcPorosity();
end