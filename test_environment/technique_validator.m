%% technique validator

% aim: to validate the recreated techniques to strengthen the comparative
% result of them. Gruber 2013 More Accurate Estimate is recreated


clc
clear

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)


N2 = 3000; % number of data points in each dimension
Ny = 30; % number of bins in relaxation time grids
sing_val=5; %sets how many singular values we compress to
tE = 200e-6;  % sample interval
T2 = logspace(log10(2*tE),log10(tE * (N2-4)),Ny); %form T2 domain, use log since will be small
tau2 = (1:N2)'*tE;  %forms measurement arrays, time tau2 domains
K2 = exp(-tau2 * (1./T2) ); % simple T2 relaxation kernel






% loaded data from Gruber 2013 

model1_true_raw = csvread('validator_models/model1_true.csv');
model1_estilt_raw = csvread('validator_models/model1_estILT.csv');
model1_estiltx_raw = csvread('validator_models/model1_estILT+.csv');
model1_true = interp1(model1_true_raw(:,1), model1_true_raw(:,2), T2, 'pchip',0)';
model1_estilt_ans = interp1(model1_estilt_raw(:,1), model1_estilt_raw(:,2), T2, 'pchip',0)';
model1_estiltx_ans = interp1(model1_estiltx_raw(:,1), model1_estiltx_raw(:,2), T2, 'pchip',0)';


%n_sigma = 0.8*trapz(T2, model1_true);
%n_sigma = sqrt(0.01*(trapz(T2, (model1_true).^2)));
%n_sigma = 0.5*trapz(T2, model1_true);
n_sigma= 0.1*trapz(model1_true);

figure(1)
clf
subplot(2,1,1)
hold on
plot(T2,model1_true)
plot(T2,model1_estilt_ans)
plot(T2,model1_estiltx_ans)
hold off
set(gca, 'XScale', 'log') 
xlim([1e-4 10])
legend('True','Estimated ILT','Estimated ILT+')


bfv_tapered = zeros(Ny,1);

num_results =100;

model1_estiltx_results = zeros(Ny, num_results);
model1_estilt_results = zeros(Ny, num_results);

%n_sigma = n_sigma * trapz(T2, model1_true)
figure(2)
clf
p1 = plot(T2,model1_estilt_ans, 'b--')
p1.LineWidth = 2;
set(gca, 'XScale', 'log') 
xlim([1e-4 10])
ylim([0 0.025])
      

omega = linspace(-0.5,1,12);
T_cutoff = [0.01 0.1 1];
    
plot_actual_omega(model1_true, T2, omega);

actual_tapered_area = plot_actual_tapered_area(model1_true, T2, T_cutoff);



for indx = 1:num_results
    noise = n_sigma*randn(N2-30,1); %assumes AWGN
    m_long = K2(31:N2,:)*model1_true + noise; % generates simulated data
    
    short_pulse_count = 10;
    m_short_meausrements = zeros(30,11);
    for spc_indx = 1:short_pulse_count+1 %include the long pulse here
       noise_short =n_sigma*randn(30,1);
       m_short = K2(1:30,:)*model1_true + noise_short;
       m_short_meausrements(:,spc_indx) = m_short;
    end
    m_short_mean = mean(m_short_meausrements')';
    m = [  m_short_mean  ;  m_long];
    
    weight_short_pulse = ones(N2,1);
    weight_short_pulse(1:30) = sqrt(11)* ones(30,1);
    m_weighted_for_short_pulse = weight_short_pulse .* m;
    K_weighted_for_short_pulse = weight_short_pulse .* K2;    
    
   %{
    figure(4)
    clf
    plot(tau2,m)
    %}

    
    [bff, compute, model1_estiltx] = iltx_estimator(bfv_tapered, ...
    m, K2, n_sigma, T2, tE, tau2); 
    
    
    
    
    %[bff, compute, model1_estilt] = ilt_estimator(bfv_tapered, m, K2, n_sigma, T2, tE);  
    [bff, compute, model1_estilt] = ilt_estimator(bfv_tapered, m_weighted_for_short_pulse ...
        , K_weighted_for_short_pulse, n_sigma, T2, tE);  
    %{
        figure(2)
        hold on

        plot(T2, model1_estilt)
        set(gca, 'XScale', 'log') 
        xlim([1e-4 10])
        ylim([0 0.025])
%}

        
    
    
    model1_estiltx_results(:,indx) = model1_estiltx;
    model1_estilt_results(:,indx) = model1_estilt;
    
end


plot_overall_est_omega(model1_estiltx, T2, omega);
est_tapered_area = plot_est_tapered_area(model1_estiltx, T2, T_cutoff);




model1_estiltx = mean(model1_estiltx_results')';

model1_estilt = mean(model1_estilt_results')';

%model1_estiltx =  mean(model1_estiltx_results')';

RMSE_ilt_techn = ((mean((model1_estilt - model1_estilt_ans).^2))^.5)

% returns percentage of error, deviation form the actual simulation
NRMSE_ilt_techn = 100*(RMSE_ilt_techn)/(mean((model1_estilt_ans).^2))^.5

MSE_ilt_techn = (mean((model1_estilt - model1_estilt_ans).^2));
MSD = mean((model1_estilt_ans - mean(model1_estilt_ans)).^2);

R2_ilt_techn = 1 - MSE_ilt_techn/MSD
max_error_ilt_techn = max(abs(model1_estilt - model1_estilt_ans))
   

RMSE_iltx_techn = ((mean((model1_estiltx - model1_estiltx_ans).^2))^.5)
NRMSE_iltx_techn = 100*(RMSE_iltx_techn)/(mean((model1_estiltx_ans).^2))^.5
max_error_iltx_techn = max(abs(model1_estiltx - model1_estiltx_ans))

MSE_iltx_techn = (mean((model1_estiltx - model1_estiltx_ans).^2));
MSD_iltx = mean((model1_estiltx_ans - mean(model1_estiltx_ans)).^2);
R2_iltx_techn = 1 - MSE_iltx_techn/MSD_iltx






figure(1)
subplot(2,1,2)
hold on
plot(T2,model1_true)
p2 = plot(T2,model1_estilt,'r-');
p2.LineWidth = 2;
p3 = plot(T2,model1_estiltx, 'g-')
p3.LineWidth = 2;
hold off
set(gca, 'XScale', 'log') 
xlim([1e-4 10])
ylim([0 0.025])
legend('True','Estimated ILT','Estimated ILT+')

figure(2)
hold on

p2 = plot(T2,model1_estilt,'r--');
p2.LineWidth = 2;
hold off
set(gca, 'XScale', 'log') 
xlim([1e-4 10])
ylim([0 0.025])
legend('Actual ILT Result','Own ILT Result')


figure(3)
clf
hold on
p2 = plot(T2,model1_estiltx_ans,'b-');
p2.LineWidth = 2;
p2 = plot(T2,model1_estiltx,'r-');

p2.LineWidth = 2;
hold off
set(gca, 'XScale', 'log') 
xlim([1e-4 10])
ylim([0 0.025])
legend('Actual ILT+ Result','Own ILT+ Result')



function plot_actual_omega(f, T2, omega)


    momentKern = ones(size(omega,2), length(f));
    poro = trapz(f); 
    indx = 1;
    for w = omega
        kern = (T2.^(w))/poro;
        momentKern(indx,:) = kern;    
        indx = indx +1;
    end
    
    actual_mom = momentKern * f;
    
    figure(55)
    clf
    subplot(2,1,1)
    hold on
    p3 = plot(omega, actual_mom, 'b--')
    p3.LineWidth = 2;
    subplot(2,1,2)
    xlabel('$\omega$')
    
    
end

function plot_overall_est_omega(f, T2, omega)


    momentKern = ones(size(omega,2), length(f));
    poro = trapz(f); 
    indx = 1;
    for w = omega
        kern = (T2.^(w))/poro;
        momentKern(indx,:) = kern;    
        indx = indx +1;
    end
    
    actual_mom = momentKern * f;
    
    figure(55)
    subplot(2,1,1)
    hold on
    p3 = plot(omega, actual_mom, 'r--')
    p3.LineWidth = 2;
    subplot(2,1,2)
    xlabel('$\omega$')
    
    
end


function actual_tapered_area = plot_actual_tapered_area(f, T2, T_cutoff)


    taperedKern = ones(size(T_cutoff,2), length(f));
    indx = 1;
    for Tc = T_cutoff
        C = 0.7213 / Tc;
        alpha = 1.572*Tc;
        beta = 0.4087 / Tc;
        gamma = 1./T2 + beta;
        kern = (C./gamma).*tanh(alpha*gamma);
        taperedKern(indx,:) = kern;
        indx = indx +1;       
    end
    
    actual_area = taperedKern * f;
    actual_tapered_area = actual_area;
    
    
    figure(56)
    clf
    subplot(2,1,1)
    hold on
    p3 = plot(T_cutoff, actual_area, 'b--')
    p3.LineWidth = 2;
    subplot(2,1,2)
    xlabel('$T_c$')
    
    
end


function est_tapered_area = plot_est_tapered_area(f, T2, T_cutoff)


    taperedKern = ones(size(T_cutoff,2), length(f));
    indx = 1;
    for Tc = T_cutoff
        C = 0.7213 / Tc;
        alpha = 1.572*Tc;
        beta = 0.4087 / Tc;
        gamma = 1./T2 + beta;
        kern = (C./gamma).*tanh(alpha*gamma);
        taperedKern(indx,:) = kern;
        indx = indx +1;       
    end
    
    actual_area = taperedKern * f;
    est_tapered_area = actual_area;
    
    
    figure(56)
    subplot(2,1,1)
    hold on
    p3 = plot(T_cutoff, actual_area, 'r--')
    p3.LineWidth = 2;
    subplot(2,1,2)
    xlabel('$T_c$')
    
    
end
