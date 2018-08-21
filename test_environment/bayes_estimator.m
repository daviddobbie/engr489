%% Bayesian BBF Estimator
% David Dobbie
% Victoria University of Wellington


function [estimate compute_time] = bayes_estimator(g_bfv, m, K, n_sigma,...
    T2, tE, prior_matrix)
    
    tic;


    N2 = length(m);
    Ny = length(T2);
    
    g_poro = ones(Ny,1);
    

    Cf = cov(prior_matrix');
    Cn = (n_sigma)^2*eye(N2);
    
    mu_prior_f = mean(prior_matrix')';
    
    R = Cf * K' * inv(K* Cf * K' + Cn);
    
    est_f = (  R * (m - K*mu_prior_f) + mu_prior_f);
    
    estimate_bfv = g_bfv' * est_f;
    estimate_porosity = g_poro' * est_f;

    estimate = estimate_bfv / estimate_porosity;
    
    %est_uncertainty = sqrt(g_bfv'  *   (Cf - R*K*Cf')  *  g_bfv);
    
    
    
    
    
    compute_time = toc;
end

