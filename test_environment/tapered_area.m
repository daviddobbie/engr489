%% Tapered Area BFF Estimator
% Programmed by: David Dobbie
% Victoria University of Wellington
%
% Recreation of:
% Estimation of petrophysical and fluid properties using integral
% transforms in nuclear magnetic resonance
% by
% F. K. Gruber et al / Journal of Magnetic Resonance 228 (2013) 104-115

function [estimate compute_time] = tapered_area(Tc, m, n_sigma,...
    T2, tE, tau2)
    
    tic;

    area = 0;
    area_uncert = 0;

    % set up in table 1 of paper
    C = 0.7213 / Tc;
    alpha = 1.572*Tc;
    beta = 0.4087 / Tc;
    gamma = 1./T2 + beta;
    N = length(m);

    % implementing (-1)^n as vector operator  
    n = floor(tau2/(2*alpha));
    n_term = ((-1).^n)';    
    k = C.*n_term'.*exp(-beta * tau2);
    estimate_tapered_area = tE * (k' * m); % eq5
    
    
   % note: we estimate porosity independnt of an ILT, diverging from the
    %       method described in paper 2. It mentions in practice that
    %       simply using TSE estimation is 'reasonable' (pg 23)
    
    [tse_porosity] = polyfit(1:100, m(1:100)', 1);
    porosity = tse_porosity(2);
    estimate = estimate_tapered_area/porosity;
    
    
    
    estimate_uncert = (n_sigma)^2 *tE * ((tE * k')*k); % eq6
    %kernel = (C./gamma).*tanh(alpha*gamma);
    compute_time = toc;
end


