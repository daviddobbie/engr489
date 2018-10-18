%% ILT BBF Estimator
% Programmed by: David Dobbie
% Victoria University of Wellington
%
% Recreation of:
% Solving Fredholm Integrals of the first Kind With
% Tensor Product Structure in 2 and 2.5 Dimensions
% by
% Venkataramanan et al

function [estimate, compute_time, f_est, estimate_bfv] = ilt_estimator(g_bfv, m, K, n_sigma,...
    T2, tE)
    tic;

    
    %ilt_estimator_multiple(g_bfv, m, K, n_sigma, T2, tE)
    
    
    Ny = length(T2);
    g_poro = ones(Ny,1);

    [m_comp k_comp] = compressData(m,K);

    N2 = length(m_comp);
    alpha = 1000;
    c = 1000*ones(N2, 1);
    f_est = c;
    
    %this is the method that prevents it being divergent
    alph_past = 0;
    indx = 0;
    %alpha = n_sigma * sqrt(N2)/ norm(c)
    while indx < 20
    %while abs(alph_past - alpha) > 1e-8 && indx < 30
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
              
        %c = newton_search(K_square, alpha, c, m_comp);
        
        %alpha = 0.65*n_sigma * sqrt(N2)/ norm(c);
        alpha = n_sigma * sqrt(N2)/ norm(c);
        %alpha = 14;
        indx = indx + 1;
    end
    alpha;
    hold off
    f_est = max(0, k_comp'*c);
    
    estimate_bfv = g_bfv' * f_est;
    estimate_porosity = g_poro' * f_est;
    estimate = estimate_bfv / estimate_porosity;
    
    compute_time = toc;
    
end

function c_new = newton_search(M, alpha, c_old, m) %method for finding approx inverse 
    c = 100*c_old;
    lambda = 1000;
    indx = 1;
    c_func_diff = (M + alpha*eye(length(m))) * c - m;
    while norm(c_func_diff) / norm(m) > 1e-6
        c_func_diff = (M + alpha*eye(length(m))) * c - m;
        c_func_2nd_diff = M + alpha*eye(length(m))  ;

        grad = inv(M + alpha*eye(length(m))) * c_func_diff;
        s = (grad'*c_func_diff) / ...
            ( pinv(grad) * c_func_2nd_diff * grad);
        
        s = 1;
        
        indx = 1;
        gamma = 0.5;
        
        cond_bool = 0;
        indx = 1;
        while (cond_bool == 0)
                    
            c_next = c - 0.5^indx * s * grad;
            
            c_func = .5 * pinv(c) * (M + alpha*eye(length(m))) * c - c' * m;
            c_func_next = .5 * pinv(c_next) * (M + alpha*eye(length(m))) * c_next - c_next' * m;
            
            
            c_func_diff_next = (M + alpha*eye(length(m))) * c_next - m;
            
            
            % setting the next c
            c_func_is_next_good = c_func_next - c_func;
            
            c_func_grad_increase = c_func_diff' * c_func_diff_next;
            
            
            if(c_func_next < c_func) 
                cond_bool = 1;
            end
            if (c_func_diff' * c_func_diff_next > 0)
                cond_bool = 1;
            end
            
            indx = indx + 1;
        end
        c = c_next;
      norm(c_func_diff -m) / norm(m); 
    end
    c_new = c;

end




function [estimate, compute_time, f_est] = ilt_estimator_multiple(g_bfv, m, K, n_sigma,...
    T2, tE)
    tic;

    Ny = length(T2);
    g_poro = ones(Ny,1);

    [m_comp k_comp] = compressData(m,K);

    N2 = length(m_comp);

    c = ones(N2, 1);
    f_est = c;
    %alpha = n_sigma * sqrt(N2)/ norm(c)   
    
    alpha = 1e4;
    
    start_alpha = [0.01 1 10 1000]
    alpha_indx = 1;
    
    
    f_est = zeros(Ny,length(start_alpha))    
    
    for alpha = start_alpha
        indx = 0;   
        while indx < 30
            K_square = k_comp* k_comp'; %recreate eq 30
            h = heaviside(K_square);
            if h == 0.5 
                h =0;
            end
            K_square = K_square .* h;
            %made symmetric and semi-positive definite

            c = inv(K_square + alpha*eye(length(m_comp))); %eq 29
            c = c'*m_comp;
            
            %f_est = max(0, k_comp'*c);      
            %alpha = 0.75*n_sigma * sqrt(N2)/ (norm(c));
            alpha = n_sigma * sqrt(N2)/ norm(c);
            indx = indx + 1;
        end
        hold off
        f_est(:,alpha_indx) = max(0, k_comp'*c); 
        
        alpha_indx = alpha_indx + 1;
    end
    
    figure(33)
    clf
    hold on
    p2 = plot(T2, f_est, '-s')
    p2(1).LineWidth = 1.5
    p2(1).LineStyle = '-.' 
    p2(1).Marker = 's'     
    p2(2).LineWidth = 1.5
    p2(2).LineStyle = '--'
    p2(2).Marker = 'o'     
    p2(3).LineWidth = 1.5
    p2(3).LineStyle = ':'
    p2(3).Marker = 'x'     
    p2(4).LineWidth = 1.5
    p2(4).LineStyle = '-.'
    p2(4).Marker = '+'     
    hold off
    xlabel('$T_2$')
    ylabel('$f_{est}$')
    set(gca, 'XScale', 'log') 
    lgnd = legend('$\alpha_\textrm{start} = 0.01$',...
        '$\alpha_\textrm{start} = 0.1$', '$\alpha_\textrm{start} = 1$',...
        '$\alpha_\textrm{start} = 10$', ...
        '$\alpha_\textrm{start} = 100$', '$\alpha_\textrm{start} = 1000$')
    lgnd.Interpreter = 'latex'
    
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

