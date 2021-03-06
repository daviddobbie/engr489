%% David Dobbie
% Victoria University of Wellington
% Recreating paper 2 (Mellin Transform of CPMG data)
%
% 
% L. Venkataramanan et al/ journal of Magnetic Resonance 206 (2010) 20-31

%Aim: Computing moments of transverse T2 relaxation time from measured CPMG
%       data. Will implement in one dimension.

addpath('../polyfit3/');



%close all
clc
clf
clear all
set(0,'defaultTextInterpreter','latex');
 
density_funcload = load('modelD.csv');

%% init variables
% number of data points in each dimension
N2 = 5000;
% number of bins in relaxation time grids
Ny = 100;      
%sets how many singular values we compress to
sing_val=20; %no singular values
tE = 500e-6; % sample interval
T2 = logspace(log10(2*tE),1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

%recreation of MODEL D pg 26

f_answer = density_funcload
f_answer = interp1(density_funcload(:,1),density_funcload(:,2),T2,'pchip')'
f_answer = f_answer./trapz(f_answer); % normalise to unity porosity
porosity =1;

%f_answer = zeros(Ny, 1);
%f_answer(75) = 1;

%{
%generate the density fn
T2_mean = 0.3
T2_var = 0.04

% normal dist answer
f_answer = normpdf(log10(T2), log10(T2_mean), sqrt(T2_var))';
%}
%{
%recration MODEL C pg 26
T2_mean = 0.018
T2_var = 0.43

% normal dist answer
f_answer = normpdf(T2, log10(T2_mean), sqrt(T2_var))';

%}

%porosity = trapz(f_answer);
%f_answer = f_answer./porosity; % normalise to unity porosity

%porosity = trapz(f_answer)
%delta distribut
%f_answer = zeros(Ny,1);
%f_answer(45) = 1;

figure(1)
clf
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f(T_2)$')
title('Correct Density Function of $T_2$');

n_stddev = 0.02


% generate measured data from density function
N_histleng = 300;

hist_data = zeros(N_histleng,2);
for hist_indx = 1:N_histleng
    [m, v] = meanAndVarEstimation(N2, Ny, sing_val, tE, f_answer, n_stddev);
    hist_data(hist_indx,1) = m;
    hist_data(hist_indx,2) = v;
end
figure(5)
clf
h1 = histogram(hist_data(:,1), 'BinWidth', 0.005);
%xlim([0.150 0.350]);
mean_actual = 10^((log10(T2))*f_answer)
estmean_avg = mean(hist_data(:,1))
estmean_stddev = 100*std(hist_data(:,2)) / estmean_avg
%estmean_med = median(hist_data(:,1))

xlabel('$T_{2,LM}$ Mellin Transfom (sec)')
ylabel('Frequency')
title('Mellin Transform Estimated  $\langle T_2 \rangle$ Relaxation Times')

figure(6)
clf
h2 = histogram(hist_data(:,2), 'BinWidth', 0.01);
%xlim([0.02 0.4]);
estvar_avg = mean(hist_data(:,2))
estvar_stddev = 100*std(hist_data(:,2))/estvar_avg
%estvar_med = median(hist_data(:,2))

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
        var = ((delta.^2)*(delta.^2)'/(gamma(omega+1)^2))*(sigma_n/poro)^2;
        %var = ((delta)*(delta)'/(gamma(omega+1)^2))*(sigma_n/poro)^2;
        var= var + ((moment - k)^2)*(sigma_p/poro)^2;
        return;
    elseif -1 < omega && omega < 0  %implement eq 22
        %{
        for indx = 2:3;
            if m(indx) > 1;
               m(indx) = 1;
            end
        end
        %}
        tau_min = tE^omega; %eq 19a
        k = tau_min / gamma(omega+1); %eq 19a
        
        n = 2:1:N-1;
        
        % eq 19c-e
        delta_1 = (0.5 * tau_min * (2^omega - 1^omega) );
        delta_mid = 0.5 * tau_min * ((n+1).^omega - (n-1).^omega);
        delta_N = (0.5 * tau_min * (N.^omega - (N-1).^omega) );
        delta = [delta_1 delta_mid delta_N];
          %{
                 hold on
         plot (delta)
         hold off
%}
        %estimate 1st derivate of M at 0 with polyfit
        coeffc = polyfit(1:100, m(1:100)', 1);
        a1 = (coeffc(1))/(tE)
        m_est = coeffc(1).*(1:100) + coeffc(2);
        a1_stddev = std(abs(m(1:100) - m_est')); %standard deviation of slope
        
        %a1 = -1;
        
        const = ((a1*omega) / (omega+1)) * tau_min^((omega+1)/omega);
        
        moment = k + (1/(gamma(omega + 1)*poro)) * (const + delta*m);
        
        
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           moment
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




function [momentT2, varT2] = meanAndVarEstimation(N2, Ny, sing_val, tE, f_answer, n_stddev)
DEBUG = 0;
T2 = logspace(log10(2*tE),1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data kernel

n = n_stddev*randn(N2,1); %generates the measurement AWG noise

M = K2*f_answer + n;
%M = M .* heaviside(M);


if DEBUG
    figure(2)
    clf
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
    clf
    plot(M_comp) %compressed M
    xlabel('Data points')
    ylabel('$M_{compressed}(t)$')
    title('Compressed Simulated Noisy Data M(t), $\sigma_{\epsilon}=0.2$');
end

porosity_answer = trapz(f_answer);

%mom = mellinTransform(M_comp, 1, tE, porosity_answer, 0.001 ,n_stddev);

% use G(omega) = ln<(T_2)^omega>, plot it

% initialie the omega_axis
omega_axis = linspace(-0.5,1,30)';

result_axis = zeros(length(omega_axis), 2);
J = zeros(1, length(omega_axis));
Sigma_G = eye(length(omega_axis));
indx_covar_mat = 1;

% printing the index for delta vector
%{
figure(66)
clf
xlabel("index - i")
ylabel("$\Delta_i$ values")
%}

for inx = 1:length(omega_axis)
    omg = omega_axis(inx);
    % certain about porosity
    [mom var] = mellinTransform(M_comp, omg, tE , porosity_answer, 0.2 ,n_stddev);
    J(indx_covar_mat) = 1/mom;
    
    Sigma_G(indx_covar_mat,indx_covar_mat) = var;
    
    result_axis(indx_covar_mat,:) = [log10(mom)  var];
    
    indx_covar_mat = indx_covar_mat +1; 
    
end 

G_covar = J * Sigma_G * J';

[omega_axis result_axis(:,1)];

r_a = result_axis(:,1);
var_r_a = result_axis(:,2);


% clean up invalid values
n1_r_a = r_a(~isnan(r_a));
n2_r_a = n1_r_a(~isinf(n1_r_a));

% clean up invalid values
n1_var_r_a = var_r_a(~isnan(r_a));
n2_var_r_a = n1_var_r_a(~isinf(n1_var_r_a));

% clean up invalid values
n1_omega_axis = omega_axis(~isnan(r_a));
n2_omega_axis = n1_omega_axis(~isinf(n1_r_a));



%{
n2_omega_axis = omega_axis;
n2_r_a = r_a;
n2_var_r_a = var_r_a;
%}
%if DEBUG
    figure(4)
    clf
    errorbar(n2_omega_axis, n2_r_a(:,1), n2_var_r_a(:,1))
    %plot(n2_omega_axis, n2_r_a(:,1))
    xlabel("$\omega^{th}$ moment")
    ylabel("$G(\omega)$ Mellin transform of data")
    title("$ G(\omega) \equiv ln \langle T_2^\omega \rangle $")
    ylim([-2 2])
    
%end

%quadratic fit on Mellin transform (2nd order polynomial)

coeffc = polyfit(n2_omega_axis, n2_r_a, 2);

% added weighting least squares for values that have a significantly
% smaller variance involved
coeffc = polyfit3(n2_omega_axis, n2_r_a , 2 , [], 1./n2_var_r_a);

var_lnT2 = abs(coeffc(1));
moment_lnT2 = coeffc(2);

momentT2 = 10^(moment_lnT2);
varT2 =  (var_lnT2);
return;
end




% Estimates the porosity and its uncertainty - correlates to area under T2
% dist.
function [poro, poro_uncert] = calcPorosity();
end