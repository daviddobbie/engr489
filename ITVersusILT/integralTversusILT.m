%% David Dobbie
% Victoria University of Wellington
% Comparion of EHT, MT versus ILT
%
% 

%Aim: Comparing the performance of linear functional transforms with the
%original ILT method in Venk 2002.



clc
clear

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',18)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.2)


%% init variables
% number of data points in each dimension
N2 = 1000;
% number of bins in relaxation time grids
Ny = 50;      
%sets how many singular values we compress to
sing_val=10; %no singular values
tE = 200e-6; % sample interval
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel

% RECREATE MODEL 1 from Gruber et al pap 4



%generate the density fn
T2_mean1 = 0.03
T2_var1 = 0.03

T2_mean2 = 0.006
T2_var2 = 0.04


f_answer = .25*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + .75*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';


f_answer = f_answer./trapz(f_answer); % normalise to unity porosity
porosity = trapz(f_answer);


%{

%generate the density fn
T2_mean1 = 0.01
T2_var1 = 0.04

T2_mean2 = 0.15
T2_var2 = 0.08


% formation of distribution
f_answer = .1*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + .25*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';



f_answer = f_answer./trapz(f_answer); % normalise to unity porosity
porosity = trapz(f_answer);
%}

%delta distribut
%f_answer = zeros(Ny,1);
%f_answer(500) = 1;


figure(1)
plot(T2, f_answer);
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('$f(T_2)$')
%title('Correct Density Function of $T_2$');


% generate the noise
noise_mean = 0;
n_std_dev = 0.1;

%calc the logarithmic mean
actualMean = exp((log(T2))*f_answer)

%--------------- running simulations and results

results_leng = 10000;
results_tpvIT = zeros(3,results_leng);
results_momIT = zeros(4,results_leng);

results_tpvILT = zeros(3,results_leng);
results_momILT = zeros(4,results_leng);

results_tpvTrue = zeros(3,results_leng);
results_momTrue = zeros(4,results_leng);


areaErrorOld = 0;
areaErrorNew = 0;


tic
for i = 1:results_leng
    i
    [tpvIT momIT tpvILT momILT tpvTrue momTrue] =  ...
        evaluateIntegralTransforms(n_std_dev, noise_mean, f_answer, ...
        K2, N2, Ny, tE, T2, tau2, porosity); 
    
    results_tpvIT(:,i) = tpvIT(:,1);
    results_momIT(:,i) = momIT(:,1);    
    results_tpvILT(:,i) = tpvILT(:,1);
    results_momILT(:,i) = momILT(:,1);    
    results_tpvTrue(:,i) = tpvTrue(:,1);
    results_momTrue(:,i) = momTrue(:,1);     
    
    
    
    
end    
toc




% tapered Areas comparison

figure(9)
clf
x_axis = 0:0.01:1;
compareTechniques(results_tpvILT(1,:), results_tpvIT(1,:), ...
        results_tpvTrue(1,:), x_axis)
                subplot(2,1,1)
        title('A')
                subplot(2,1,2)
        title('B')
figure(10)
clf
x_axis = 0:0.01:1;
compareTechniques(results_tpvILT(2,:), results_tpvIT(2,:), ...
        results_tpvTrue(2,:), x_axis)
                subplot(2,1,1)
        title('A')
                subplot(2,1,2)
        title('B')  
figure(11)
clf
x_axis = 0:0.01:1;
compareTechniques(results_tpvILT(3,:), results_tpvIT(3,:), ...
        results_tpvTrue(3,:), x_axis)
                subplot(2,1,1)
        title('C')
                subplot(2,1,2)
        title('D')
 
% moments comparison
figure(12)
clf
x_axis = 0:0.01:1;
compareTechniques(results_momILT(1,:), results_momIT(1,:), ...
        results_momTrue(1,:), x_axis)
            subplot(2,1,1)
    title('MT $\omega$ = -0.5')
           subplot(2,1,1)
        title('A')
                subplot(2,1,2)
        title('B')
    
figure(13)
clf
x_axis = 0:0.01:1;
compareTechniques(results_momILT(2,:), results_momIT(2,:), ...
        results_momTrue(2,:), x_axis)
            subplot(2,1,1)
        title('MT $\omega$ = 0.1')
figure(14)
clf
x_axis = 0:0.01:1;
compareTechniques(results_momILT(3,:), results_momIT(3,:), ...
        results_momTrue(3,:), x_axis) 
        subplot(2,1,1)
        title('MT $\omega$ = 0.5')
        
        subplot(2,1,1)
        title('C')
                subplot(2,1,2)
        title('D')
figure(15)
clf
x_axis = 0:0.01:1;
compareTechniques(results_momILT(3,:), results_momIT(3,:), ...
        results_momTrue(4,:), x_axis) 
            subplot(2,1,1)
        title('MT $\omega$ = 1')

    
    
    
    

%% function definitions:

%
%
% INPUTS:
% IT
% ILT
% True
% xaxis to plot on
% OUTPUTS:
% creates plot of estimations for each technique
function compareTechniques(ILT, IT, true, x_axis)

    min_val = min([ILT IT]);
    max_val = max([ILT IT]);

    subplot(2,1,1)
    
    hold on
    h1 = histogram(IT, 30, 'BinLimits', [min_val max_val]);
    p = plot([mean(true) mean(true)], [0 max(h1.Values)]);
    p.LineWidth = 2;
    p.Color = 'k';
    p = plot([mean(IT) mean(IT)], [0 max(h1.Values)]);
    p.LineWidth = 2;
    p.Color = 'r';
    

    title('IT')
    xlabel('$T2$ [s]')
    ylabel('Frequency')
    
    ITError = normalisedRootMeanSquareError(mean(IT), mean(true));
    ITStdDev = (std(IT)/mean(IT))*100;
    
    str = sprintf('NRMSE = %4.2f %% \nNSD = %4.2f %%',ITError, ITStdDev);
    dim = [0.832196669275554 0.900822892171908 0.172638431822438 0.0906800988038782];
%     annotation('textbox',...
%     dim,...
%     'String',[str],...
%     'FontSize',8,...
%     'FontName','Times New Roman',...
%     'LineStyle','none',...
%     'FitBoxToText','on');

     hold off  
    
    subplot(2,1,2)
    
    hold on
    
    h2 = histogram(ILT, 30, 'BinLimits', [min_val max_val]);
    p = plot([mean(true) mean(true)], [0 max(h2.Values)]);
    p.LineWidth = 2;
    p.Color = 'k';
    p = plot([mean(ILT) mean(ILT)], [0 max(h2.Values)]);
    p.LineWidth = 2;
    p.Color = 'r';
    title('ILT')
    xlabel('$T2$ [s]')
    ylabel('Frequency')

    ILTError = normalisedRootMeanSquareError(mean(ILT), mean(true));
    ILTStdDev = (std(ILT)/mean(ILT))*100;
    
    str = sprintf('NRMSE = %4.2f %% \nNSD = %4.2f %%',ILTError, ILTStdDev);
    dim = [0.826806127531373 0.438287153502255 0.172638431822438 0.0906800988038782];
%     annotation('textbox',...
%     dim,...
%     'String',[str],...
%     'FontSize',8,...
%     'FontName','Times New Roman',...
%     'LineStyle','none',...
%     'FitBoxToText','on');

    hold off

end



% Estimation of the density function from measured data. Returns two
% results, the ILT method in Venk. 2002 (old) and the ILT+ method in Gruber
% 2013 (new).
% INPUTS: 
%    m
% OUTPUTS:
%    area
function [tpdAreasVectIntegralTransform, momentVectIntegralTransform, ...
    tpdAreasVectILT, momentVectILT, tpdAreasTrue, momentVectTrue] ...
    = evaluateIntegralTransforms(n_std_dev, noise_mean, f_answer, K2, ...
    N2, Ny, tE, T2, tau2, porosity)


    %omega = linspace(0.1,1,4);
    omega = [-0.5 0.1 0.5 1];
    T_cutoff = [0.01 0.1 1];

    noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
    m = K2*f_answer + noise;  
    %{
        figure(2)
        clf
        hold on
        plot(tau2, m);
        hold off
        %title('Measured Signal')
        xlabel('t [s]')
        ylabel('$m(t)$')
    %}
    
    %---------- ESTIMATION with only measured data
    
    % estimate tapered areas for different cut off times
    tpdAreasVectIntegralTransform = ones(size(T_cutoff,2), 2);
    
    indx = 1;
    for Tc = T_cutoff
        [tpd, tpd_var] = exponentialHaarTransform(tE, m, n_std_dev, Tc, T2,tau2);
        
        k = exponentialHaarTransform(tE, ones(size(m)), n_std_dev, Tc, T2,tau2);
        tpd_var = n_std_dev^2 *(norm(k))^2;
        tpdAreasVectIntegralTransform(indx,:) = [tpd, sqrt(tpd_var)];
        
        indx = indx + 1;
    end
    tpdAreasVectIntegralTransform;

    % estimate different moments of the T2 distribution
    momentVectIntegralTransform = ones(size(omega,2), 2);

    indx = 1;
    for w = omega

        [mom, mv] = mellinTransform(m, w, tE, porosity, 0, n_std_dev);
        
        % to estimate the uncertainty of the moment (eq 11)
        k = mellinTransform(ones(size(m)), w, tE, porosity, 0, n_std_dev);
        mom_var = n_std_dev^2 *(norm(k))^2;
        
        momentVectIntegralTransform(indx,:) = [mom, sqrt(mom_var)]; %%since we use uncertainty
        indx = indx + 1;
    end

    momentVectIntegralTransform;

    % create compressed m vector of values for optimisation
    [m_comp, k_comp] = compressData(m,K2,10);
    m_comp = m_comp'; % compressed m vector

    
    f_est_old = optimisationInverseTransform(m_comp, k_comp, eye(size(m_comp,2)), n_std_dev);

    %----------- ESTIMATION with the ILT
    
    
    tpdAreasVectILT = ones(size(T_cutoff,2), 1);
    indx = 1;
    for Tc = T_cutoff
        [tpd] = actualTaperedArea(f_est_old, Tc,T2);
        tpdAreasVectILT(indx,:) = [tpd];
        indx = indx + 1;
    end    
    
    momentVectILT = ones(size(omega,2), 1);
    indx = 1;
    for w = omega
        mom = actualMoment(f_est_old, T2, w);
        momentVectILT(indx) = mom;
        indx = indx + 1;
    end

    %----------- TRUE DISTIRBUTION ANSWERS
    tpdAreasTrue = ones(size(T_cutoff,2), 1);
    indx = 1;
    for Tc = T_cutoff
        [tpd] = actualTaperedArea(f_answer, Tc,T2);
        tpdAreasTrue(indx) = [tpd];
        indx = indx + 1;
    end    
    
    momentVectTrue = ones(size(omega,2), 1);
    indx = 1;
    for w = omega
        mom = actualMoment(f_answer, T2, w);
        momentVectTrue(indx) = mom;
        indx = indx + 1;
    end    
    
    
end


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

% Calculate the kernel used to compute the tapered area in the T2 domain.
% This is the exponetial Haar transform
% INPUTS: 
%    Tc = the bound and fluid fraction point (in time)
%    T2 = T2 relaxation axis
% OUTPUTS:
%    kern = the kernel of the tapered heaviside function used.
%  
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

% Calculate the actual moment of the distribution (fully polarised).
% uses the definition of moment to compute it.
% INPUTS: 
%    f = Ny x 1 desnity function vector
%    T2 = T2 relaxation axis
%    omega = the wth moment being calculated
% OUTPUTS:
%    moment = the density function moment
%  

function moment = actualMoment(f, T2, omega)
    porosity = trapz(f);
    top_arg = (T2.^omega).*f';  
    moment = trapz(top_arg) / porosity;
end


% Calculates the normalised root mean square error of two different
% functions. Quantifies correctness.
% INPUTS: 
%    est = estimate function
%    true = the true function
% OUTPUTS:
%    nrmse = normalised root mean square error of the dist.
% 
function nrmse = normalisedRootMeanSquareError(est, true)
    s =(est - true).^2;
    ms = s;
    rms = sqrt(ms);
    nrmse = (rms/true)*100;
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

    %svd results sorted by magnitude singular values. i.e we only have to
    %truncate to s1 rowsxcols.
    
    [U2, S2, V2] = svd(K);    
    %only leave trunc number of largest values. This removes small weighted
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
    
    %{
    figure(3)
    plot(M_compressed);
    title('Compressed Measured Signal')
    xlabel('Time 2  $ \tau_2 $ [s]')
    ylabel('$M_c(t)$')
    %}
    
end


% Optimisation function. This function uses regularisation and optimisation
% technqiues with a weighted matrix and the cost matrices to apply an ILT
% (or ILT+). This is adapted from Venk 2002. solving Fredholm integrals
% INPUTS: 
%    G = N x 1 linear functional and compressed measurement data 
%    K = compressed input kernel
%    L = mapping matrix conversion from T2 to t space
%    W = weighting matrix
% OUTPUTS:
%    f_est = the estimated density function
% 
function f_est = optimisationInverseTransform(G, L, W, n_std_dev)
    alpha = 1000;
    Ny = size(L,2);
    c = ones([length(G)  1]);

    % lexiographical sorting of the m matrix
    %Glex = sortrows(reshape(G,1,[]).');

    G = W*G;

    L = W*L;
    %K2 = sortrows(reshape(K2,Ny,N2).');

    alpha_hist = [];
    f_est = c;

    %this is the method that prevents it being divergent
    for i=1:20
        
        %L_square = L*L'; 
        stepFnMatrix = (heaviside(L'*c))'.* eye(Ny);
        L_square = L *stepFnMatrix * L';       %recreate eq 30
        %made symmetric and semi-positive definite
        c = inv(L_square + alpha*eye(length(G))); %eq 29
        c = c*G;
        alpha =  n_std_dev * sqrt(length(G))/ norm(c); %implement eq 17 BRD paper  
        %plot(c) 
    end
    hold off

    f_est = L'*c;
    under = min(f_est);
    f_est = f_est - under;
    f_est = f_est ./ trapz(f_est); %normalise to unity


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
        tau_min = (tE)^omega; %eq 19a
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
              
        %moment = 1/(gamma(omega)*poro) * tE*(((tE:tE:(N)*tE).^(omega - 1))*m)
        %eq 23
        var = (sum(delta.^2)/gamma(omega+1)^2)*(sigma_n/poro)^2;
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
        
       % moment =   (1/(gamma(omega + 1)*poro)) * (const + delta*m);       
        
        
       if moment < 5e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            %disp('delta * m')
            (delta*m);
            %disp('1/gamma(omg + 1)')
            1/(gamma(omega + 1)*poro);
            moment = NaN;
        end
        
        
        
        var = (sum(delta.^2)/(gamma(omega+1))^2)*(sigma_n/poro)^2 + var;
        var= var + (moment - k)^2*(sigma_p/poro)^2;
        var = var + ((omega*tau_min^((omega+1)/omega))/gamma(omega + 2)) * (a1_stddev / poro)^2;
        return;
    else %allows outside case, shouldn't be called
        moment = 0;
        var = 0;
    end
end


