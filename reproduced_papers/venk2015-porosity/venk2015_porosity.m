%% David Dobbie
% Victoria University of Wellington
% Recreating paper 5 
% New Method to Estimate Porosity More Accurately from NMR Data with Short
% Relaxation Times
% 
% L. Venkataramanan et al/ Petrophysics Vol 56 no 2 April 2015 Pg 147-157

%Aim: Estimating Porosity by using a correcting factor with a known
%porosity sensitivity curve

% algorithm goes as follows:
% 1) compute ILT for several noise realisations
% 2) Compute the porosity sensitivity curve for all of these realisations
% 3) Add the correction to all of the ILT data points for short relaxation
% times - so we get better estimate of amount of bound fluid in the rock

clc
clf
clear

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesTitleFontSizeMultiplier', 1)
set(0,'defaultAxesFontSize',14)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.1)


%loading M4 dist from paper 2015 porosity estimation
density_funcload = load('datasets\m4.csv');
%density_funcload(:,2) = density_funcload(:,2)
[C,ia,ic]  = unique(density_funcload(:,1)),'stable';
density_funcload = density_funcload(ia,:);
%density_funcload = density_funcload(:~non_unique,:);
%{
figure(3)
clf
plot(density_funcload(:,1), density_funcload(:,2))
set(gca, 'XScale', 'log')
%}
%% Step 0: intialise variables

% number of data points in each dimension
N2 = 1000;
% number of bins in relaxation time grids
Ny = 30;      
%sets how many singular values we compress to
sing_val=5; %no singular values
tE = 10e-6;
%tE = 200e-6; % sample interval
T2 = logspace(log10(300e-6),log10(3),Ny); %form T2 domain, use log since will be small
%T2 = logspace(-5,1,Ny);
%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*tE;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation kernel




% RECREATE MODEL 2 from Gruber et al pap 4

%{
%generate the density fn
T2_mean1 = 0.05
T2_var1 = 0.03

T2_mean2 = 0.006
T2_var2 = 0.04

% formation of distribution
f_answer = .25*normpdf(log10(T2), log10(T2_mean1), sqrt(T2_var1))';
f_answer = f_answer + .75*normpdf(log10(T2), log10(T2_mean2), sqrt(T2_var2))';
%f_answer = ones(Ny,1)


porosity = 20
f_answer = porosity*f_answer./trapz(f_answer); % normalise to unity porosity
%}

f_answer = density_funcload
f_answer = interp1(density_funcload(:,1),density_funcload(:,2),T2,'pchip')'
f_answer = 0.20*f_answer./trapz(f_answer)

porosity = trapz(f_answer);


figure(3)
clf
plot(T2, f_answer)
set(gca, 'XScale', 'log')



noise_mean = 0;
n_std_dev = 0.2;
f_calibrate = eye(Ny);
%f_calibrate = f_calibrate./trapz(f_calibrate);



%% Step 1: ILT (use BRD)
results_leng = 10;
bins_ILTold = zeros(Ny,results_leng*Ny);




for idx = 1:results_leng
    
    for eachDelta = 1:Ny
        f_cal_row = f_calibrate(:,eachDelta);
        [f_est_ilt] = estimateDensityFunction(n_std_dev, noise_mean,  ... 
        f_cal_row, K2, N2, Ny, tE, T2, tau2, porosity,sing_val, 10); 

        eachDelta + Ny*(idx-1)
        
        mask = zeros(Ny,1);
        mask(eachDelta)  =1;
        %bins_ILTold(:,eachDelta + Ny*(idx-1)) = mask.*f_est_ilt;
        bins_ILTold(:,eachDelta + Ny*(idx-1)) = f_est_ilt;
        
            figure(50)
            
            clf
            hold on
            stem(T2, f_cal_row,'-b');
            stem(T2, f_est_ilt,'-r');
            hold off
            set(gca, 'XScale', 'log')
            xlabel('$T_2(s)$')
            ylabel('$f(T_2)$')
            title('Density Function of $T_2$');
            legend('True','Estimated ILT')

            %pause(0.005)
        
    end
    
end  





%% Step 2: Calc porosity curve

bias_T2 = (sum(bins_ILTold,2)/results_leng)'-1;


%bias_T2 = (mean(bins_ILTold)) - mean(f_calibrate);

%bias_T2 = bias_T2 ./bias_T2(8);


%bias_T2_augment = [bias_T2(1:9)  bias_T2(10)*ones(1,21)];



figure(1)
clf
hold on
plot(T2, bias_T2 + 1)
%plot(T2, bias_T2_augment + 1)
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('Sensitivity')
legend('Normal','Augmented')


%bias_T2 = bias_T2_augment;




%% Step 3: Compute and add correction
 
%bias_T2 = [bias_T2(1:10)  bias_T2(11)*ones(1,20)];

r_t2 = (bias_T2 +1)./n_std_dev;
correction_T2 = 1./(1 + bias_T2 .* (r_t2 ./ (mean(r_t2) + r_t2)    ));

correction_T2_simple = 1./(1+bias_T2);

figure(3)
clf
hold on
plot(T2, correction_T2)
%plot(T2, correction_T2_simple)
hold off
set(gca, 'XScale', 'log')
xlabel('$T_2(s)$')
ylabel('Sensitivity')
legend('correction', 'simple correction')


N_p_est = 500;

overall_corrected_p = zeros(1,N_p_est);
overall_old_p = zeros(1,N_p_est);
overall_answer_p = trapz(f_answer)


n_std_dev = 0.2.*trapz(f_answer);

for el = 1:N_p_est
    
    [f_est_ilt] = estimateDensityFunction(n_std_dev, noise_mean,  ... 
    f_answer, K2, N2, Ny, tE, T2, tau2, porosity,sing_val, -1);

    r_t2 = (f_est_ilt')./n_std_dev;
    correction_T2 = 1./(1 + bias_T2 .* (r_t2 ./ (mean(r_t2) + r_t2)    ));

    old =  f_est_ilt;
    corrected_simple = correction_T2_simple' .* f_est_ilt;
    corrected = correction_T2' .* f_est_ilt;    

    corrected_porosity = trapz(corrected);
    old_porosity = trapz(old);

    overall_corrected_p(el) = corrected_porosity;
    overall_old_p(el) = old_porosity;
 
    %corrected = corrected / trapz(corrected);
    
    figure(4)
    clf
    hold on
    plot(T2, f_answer,'-b');
    plot(T2, old,'-r');
    %plot(T2, corrected_simple,'-k');
    plot(T2, corrected,'-g');
    hold off
    set(gca, 'XScale', 'log')
    xlabel('$T_2(s)$')
    ylabel('$f(T_2)$')
    %legend('True','ILT', 'Simple Correction', 'Correction')
    legend('True','ILT', 'Correction')
    


end

std_corrected = 100*std(overall_corrected_p)/ overall_answer_p;
bias_corrected = 100*abs(abs(overall_answer_p - mean(overall_corrected_p))/overall_answer_p);


std_old = 100*std(overall_old_p)/ overall_answer_p;
bias_old = 100*abs(abs(overall_answer_p - mean(overall_old_p))/overall_answer_p);


figure(5)
clf
hold on
plot(bias_corrected, std_corrected,'.b', 'MarkerSize', 20)
plot(bias_old, std_old, '.r', 'MarkerSize', 20)
hold off
xlabel('Bias $\frac{B_\phi}{\phi_T} \times 100$');
ylabel('Imprecision $\frac{\sigma_\phi}{\phi_T} \times 100$');
legend('corrected','old')
xlim([0 60])
ylim([0 50])




%% FUNCTION
% Estimation of the density function from measured data. Returns 
% result, the ILT method in Venk. 2002
% INPUTS: 
%    noise standard deviation
%    noise mean
%    density function answer
%    K2 kernel
%    N2 size of t axis
%    Ny size of T2 axis
%    tE time between samples
%    T2 relaxation time axis
%    tau2 time axis
%    porosity of density function
% OUTPUTS:
%    f_est_old estimation of density function with ILT old method
function [f_est_old] = estimateDensityFunction(n_std_dev, ...
    noise_mean, f_answer, K2, N2, Ny, tE, T2, tau2, porosity, sing_val, alpha)

    noise = n_std_dev*normrnd(noise_mean, 1, [N2 ,1]);
    m = K2*f_answer + noise;  
 


    % create compressed m vector of values for optimisation
    [m_comp, k_comp] = compressData(m,K2,sing_val);

    m_comp = m_comp'; % compressed m vector 
    f_est_old = optimisationInverseTransform(m_comp, k_comp, eye(size(m_comp,1)), n_std_dev, alpha);

    
    figure(66)
    clf
    hold on
    plot(tau2, m)
    c = plot(tau2, K2*f_est_old) ;
    c.Color = 'k';
    c.LineStyle = '--';
    c.LineWidth = 2.5;
    xlabel('$T_2(s)$')
    ylabel('$f(T_2)$')
    ylim([-0.5 1.5])
    
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
function f_est = optimisationInverseTransform(G, L, W, n_std_dev, alpha)
    manualRegularisation = 1;
    if(alpha <= 0)
        manualRegularisation = 0; % invalid alpha value
    else
        manualRegularisation = 1; % sets manual regularisation
    end
        
    G = W*G;
    %{
    if(length(G) > 11)
        figure(34)
        clf
        plot(abs(G))   
    end
    %}
    L = W*L;    
    
    N = nnz(W);
    
    Ny = size(L,2);
    c = ones([length(G)  1]);

    % lexiographical sorting of the m matrix
    %Glex = sortrows(reshape(G,1,[]).');


    %K2 = sortrows(reshape(K2,Ny,N2).');

    alpha_hist = [];
    f_est = c;

    %this is the method that prevents it being divergent
    for i=1:30
        
        %L_square = L*L'; 
        %stepFnMatrix = (heaviside(L'*c))'.* eye(Ny);
        %stepFnMatrix = heaviside(L'*L);
        
        %L_square = L *stepFnMatrix * L';       %recreate eq 30
        L_square = L* L';
        L_square = L_square .* heaviside(L_square);
        all(eig(L_square) >= eps); %verify is positive semi-definite
        %made symmetric and semi-positive definite
        c = inv(L_square + alpha*eye(length(G))); %eq 29
        c = c*G;
        
        if(manualRegularisation == 0)
            alpha =  n_std_dev * sqrt(N)/ norm(c);
        end
        %alpha_hist = [alpha_hist alpha];
        %alpha =  n_std_dev * sqrt(size(nnz(W),1))/ norm(c); %implement eq 17 BRD paper  
        %plot(c) 
    end
    hold off
    %{
    figure(77)
    plot(alpha_hist)
    %}
    f_est = L'*c;
    %f_est = abs(f_est);
    %under = min(f_est);
    %f_est = f_est - under;
    %f_est = f_est ./ trapz(f_est); %normalise to unity


end
