%% ILT BBF Estimator
% Programmed by: David Dobbie
% Victoria University of Wellington
%
% Recreation of:
% Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions
% by
% Venkataramanan et al


function [estimate, compute_time, f_est] = ilt_estimator(g_bfv, m, K, n_sigma,...
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
    alph_past = 0;
    indx = 0;
    while indx < 15
        alph_past = alpha;
        abs(alph_past - alpha);
        K_square = k_comp* k_comp'; %recreate eq 30
        h = heaviside(K_square);
        if h == 0.5
            h =0;
        end
        K_square = K_square .* h;
        %made symmetric and semi-positive definite
        c = inv(K_square + alpha*eye(length(m_comp))); %eq 29
        c = c'*m_comp;
        alpha = n_sigma * sqrt(N2)/ norm(c);
        %alpha = 10;
        indx = indx + 1;
    end

    hold off
    f_est = max(0, k_comp'*c);
    
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

    
    sing_val = 20;
    
    


    %svd results sorted by magnitude singular values. i.e we only have to
    %truncate to s1 rowsxcols.
    
    [U2, S2, V2] = svd(K); 
    
    sing_val = find_condition_threshold(S2, 1000); % to set condition required
    
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


function cond_indx = find_condition_threshold(S, cond_barrier)
    cond_indx = min(size(S));
    
    while cond(S(1:cond_indx,:)) > cond_barrier
        cond_indx = cond_indx - 1;
    end
    cond_indx = cond_indx +1; 
end

