%% ILT+ BBF Estimator
% Programmed by: David Dobbie
% Victoria University of Wellington
%
% Recreation of:
% A More accurate estimate of T2 distribution from direct analysis of NMR
% measurements
% by
% F. K. Gruber et al / Journal of Magnetic Resonance 228 (2013) 95-103

function [estimate, compute_time ,f_est] = iltx_estimator(g_bfv, m, K, n_sigma,...
    T2, tE, tau2)

    tic
    

    
    Ny = length(T2);
    g_poro = ones(Ny,1);
    
    
    omega = linspace(-0.5,1,12);
    T_cutoff = [0.01 0.1 1];
    
    % estimate tapered areas for different cut off times
    tpdAreasVect = ones(size(T_cutoff,2), 2);
    tpdAreaKern = ones(size(T_cutoff,2), Ny);
    indx = 1;
    for Tc = T_cutoff
        [tpd, tpd_var] = exponentialHaarTransform(tE, m, n_sigma, Tc, T2,tau2);
        
        k = exponetialHaarTransformKernel(Tc, T2);
        
       % tpd_var = n_sigma^2 *(norm(k))^2;
        
        
        tpdAreasVect(indx,:) = [tpd, sqrt(tpd_var)];
        tpdAreaKern(indx,:) = exponetialHaarTransformKernel(Tc,T2);
        indx = indx + 1;
    end
    tpdAreasVect;

    % estimate different moments of the T2 distribution
    momentVect = ones(size(omega,2), 2);
    momentKern = ones(size(omega,2), Ny);

    indx = 1;
    
    % note: we estimate porosity independnt of an ILT, diverging from the
    %       method described in paper 2. It mentions in practice that
    %       simply using TSE estimation is 'reasonable' (pg 23)
    
    [tse_porosity] = polyfit(1:30, m(1:30)', 1);
    porosity = tse_porosity(2);
    
    for w = omega

        kern = (T2.^(w))/porosity;
        momentKern(indx,:) = kern;
        [mom, mom_var] = mellinTransform(m, w, tE, porosity, n_sigma, n_sigma);
        
        
        % to estimate the uncertainty of the moment (eq 11)
        %k = mellinTransform(ones(size(m)), w, tE, porosity, 0, n_sigma);
        %mom_var = (n_sigma)^2 * norm(k)^2;
        
        momentVect(indx,:) = [mom, sqrt(mom_var)]; %%since we use uncertainty
        
        
        indx = indx + 1;
    end

    momentVect;

    % weighting measurements with short pulse increased SNR
    
    weight_short_pulse = ones(size(m));
    weight_short_pulse(1:30) = sqrt(11)* ones(30,1);
    m_weighted_for_short_pulse = weight_short_pulse .* m;
    K_weighted_for_short_pulse = weight_short_pulse .* K;    
    
    
    %m_weighted_for_short_pulse = m;
    %K_weighted_for_short_pulse = K;
    
    % create compressed m vector of values for optimisation
    [m_comp, k_comp] = compressData(m_weighted_for_short_pulse, ...,
        K_weighted_for_short_pulse,10);

    m_comp = m_comp'; % compressed m vector


    G_opt = [m_comp'; momentVect(:,1) ; tpdAreasVect(:,1)]; %eq 13 pap4
    
    
    L_opt = [k_comp ; momentKern ; tpdAreaKern]; % eq 14 pap 4
    %W_vect = [1*(ones(size(m_comp')))/n_sigma; 0./momentVect(:,2) ; ...
     %   0./tpdAreasVect(:,2)];
     


    W_vect = [1*ones(size(m_comp'));  1*n_sigma./(momentVect(:,2)) ; ...
       n_sigma./(tpdAreasVect(:,2))]  ;    

    n_sigma = n_sigma * (1/2)
    %n_sigma = n_sigma * sqrt(1/3)   
   
   %n_sigma = n_sigma*mean(1./find(W_vect))
   
   
   figure(55)
   subplot(2,1,1)
   hold on
   plot(omega, momentVect(:,1))
   hold off
   subplot(2,1,2)   
   hold on
   plot(omega, 1./momentVect(:,2))
   hold off
   
   figure(56)
   subplot(2,1,1)
   hold on
   plot(T_cutoff, tpdAreasVect(:,1))
   hold off
   subplot(2,1,2)   
   hold on
   plot(T_cutoff, 1./tpdAreasVect(:,2))
   hold off
   
   
    %W_vect = [1*(ones(size(m_comp'))); 0./momentVect(:,2) ; ...
     % ones(length(tpdAreasVect),1)/2];
    
    
    
    W_vect;
    
    W_opt = W_vect .* eye(size(G_opt,1));   
     
 
    
    
    
    n_sigma_avg = mean(1./W_vect);
    n_sigma_avg = n_sigma;
    
    f_est = (0.58)*optimisationInverseTransform(G_opt, L_opt, W_opt, n_sigma_avg);
% magic number comes from scaling for tapered areas. It does not seem to
% behave in the context of the optimisation problem.
    
    estimate_bfv = g_bfv' * f_est;
    estimate_porosity = g_poro' * f_est;
    estimate = estimate_bfv / estimate_porosity;
    
    compute_time = toc;
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
    
    area_uncert = sqrt((sigma_n)^2 *tE * ((tE * k')*k)); % eq6
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
    K_compressed = (S2c*V2c');
    M_compressed = (M'*U2c)';
    
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
function f_est = optimisationInverseTransform(G, L, W, n_sigma)
    alpha = 1000;
        
    G = W*G;
    
    figure(44)
    clf
    plot(abs(G))
    
    figure(22)
    clf
    hold on
    set(gca, 'XScale', 'log') 
    ylim([0 0.05])
    
    L_normal = L;
    L = W*L;  
    
    N = nnz(W);
    N = sum(diag(W)>1e-1);

    c = ones([length(G)  1]);
    
    f_est = c;
    %this is the method that prevents it being divergent
    alph_past = 0;
    indx = 0;
    %while abs(alph_past - alpha) > 1e-8 && indx < 30
    while  indx <20
            alph_past = alpha; 
            L_square = L* L';      
            h = heaviside(L_square);
            if h == 0.5
                h =0;
            end
            L_square = L_square .* h;     %recreate eq 30       
            %made symmetric and semi-positive definite
            c = inv(L_square + alpha*eye(length(G))); %eq 29
            c = c'*G;
            alpha = 1*n_sigma *sqrt(N)/ norm(c);
            figure(22)
            plot(max(L'*c,0))
            set(gca, 'XScale', 'log') 
            pause(0.1)
            
            %alpha = 14;
            indx = indx + 1;
            %alpha = median(1./diag(W)) * sqrt(N)/ norm(c);
            %alpha =15;
    end

    hold off
    

    
    f_est = max(L'*c,0);
    
    f_est = f_est - min(f_est);
    
    %f_est = max(L'*c,0);
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
        
%         hold on
%         plot (delta)
%         hold off
        
        %delta = [(0.5*tau_min(2.^omega-1.^omega)) delta (0.5*tau_min*(N.^omega-(N-1).^omega))];
        
        %omega
        
        % note that erroneous values are apparent for a negative
        % measurement (leads to complex result from log)
        moment = k + 1/(gamma(omega + 1)*poro) * (delta*m); % eq18
        %{
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           omega
        end
        %}

        
        if moment < 1e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            %disp('delta * m');
            (delta*m);
            %disp('1/gamma(omg + 1)');
            1/(gamma(omega + 1)*poro);
            moment = 1;
        end
           
        
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
        a1 = (coeffc(1))/(tE);
        m_est = coeffc(1).*(1:100) + coeffc(2);
        a1_stddev = std(abs(m(1:100) - m_est')); %standard deviation of slope
        
        %a1 = -1;
        
        const = ((a1*omega) / (omega+1)) * tau_min^((omega+1)/omega);
        
        moment = k + (1/(gamma(omega + 1)*poro)) * (const + delta*m);
        
        %{
        if moment < 1e-1
           disp('OMEGA CAUSED ERROR'); 
           moment
        const + delta*m
        const
        delta*m
        omega
        end
        %}
        if moment < 5e-1 || abs(moment) == Inf % computation breaks, is negative
            k;
            moment = 1;
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






