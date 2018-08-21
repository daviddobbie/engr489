%% IKT BBF Estimator
% Programmed by: David Dobbie
% Victoria University of Wellington
%
% Recreation of:
% Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions
% by
% Venkataramanan et al


function [estimate compute_time] = ilt_estimator(g_bfv, m, K, n_sigma,...
    T2, tE)
    tic;

    Ny = length(T2);
    g_poro = ones(Ny,1);

    [m_comp k_comp] = compressData(m,K);

    N2 = length(m_comp);
    alpha = 1000;
    c = ones(N2, 1);
    f_est = c;

    %this is the method that prevents it being divergent
    for i=1:20
             
        K_square = k_comp* k_comp'; %recreate eq 30
        K_square = K_square .* heaviside(K_square);
        %made symmetric and semi-positive definite
        c = inv(K_square + alpha*eye(length(m_comp))); %eq 29
        c = c'*m_comp;
        alpha = n_sigma * sqrt(N2)/ norm(c);
    end
    hold off
    f_est = k_comp'*c;
    
    estimate_bfv = g_bfv' * f_est;
    estimate_porosity = g_poro' * f_est;
    estimate = estimate_bfv / estimate_porosity;
    
    compute_time = toc;
    
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

function [M_compressed K_compressed] = compressData(M,K)
    N = length(M);

    
    sing_val = 10;
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
