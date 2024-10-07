%% subproblem 1    
% 3.1 p_d optimziation
acc = 0.1;

Zsk = Sk * 1000;

p_d_ori = p_d * sigma_n;
beta_n = beta / sigma_n;
lambda_n = lambda / sigma_n;

omega_k = Zsk ./ (1 - eta) / Bandwidth;
R_k_req = omega_k ./ ((Delay_req - delay_1 - delay_3) * 1000);
gamma_k_req = 2.^(R_k_req) - 1;

for m0 = 1:M

    varphi_m = (sum(sqrt(p_d_ori .* lambda_n)) - sqrt(p_d_ori(m0,:) .* lambda_n(m0,:)))';
    nu_m = (sum(p_d_ori' * beta_n) - sum(p_d_ori(m0,:)' * beta_n(m0,:)) + 1)';

    p_dm0_fixed = p_d_ori(m0,:)';
    
    OB_ori1 = sum((omega_k .* (p_dm0_fixed + (sum(p_d_ori) - p_d_ori(m0,:))')) ./ (log(N * p_dm0_fixed .* lambda_n(m0,:)' + sum(p_dm0_fixed * beta_n(m0,:))' + 2*N*sqrt(p_dm0_fixed) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2) - (log(sum(p_dm0_fixed * beta_n(m0,:))' + nu_m) / log(2))));
    
    if m0 == 1
        ori_outer = OB_ori1;
    end
    
    for i2 = 1:10
        Delta_g2_vector = beta_n(m0,:)' ./ (log(2) * (sum(p_dm0_fixed * beta_n(m0,:))' + nu_m));
        Delta_g2_matrix = repmat(Delta_g2_vector, 1, K); % 梯度矩阵，每一行代表一个用户，这一行的所有列代表这个用户对应的K个梯度矩阵（每一行的元素都是一样的）；
        g2k_const = log(sum(p_dm0_fixed * beta_n(m0,:))' + nu_m) / log(2);
        
        p_dm0_inver = p_dm0_fixed;
        
        AK_inver = omega_k .* (p_dm0_inver + (sum(p_d_ori) - p_d_ori(m0,:))');
        rk_inver = log(N * p_dm0_inver .* lambda_n(m0,:)' + sum(p_dm0_inver * beta_n(m0,:))' + 2*N*sqrt(p_dm0_inver) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2) - (g2k_const + Delta_g2_matrix * (p_dm0_inver - p_dm0_fixed));
        
        z_k_opt = sqrt(rk_inver) ./ AK_inver;
        OB_ori3 = sum(AK_inver ./ rk_inver);
        OB_ori2 = OB_ori3;
%         fprintf('m0: %d, i2: %d, inverse Q before OB_ori3 ***: %f\n', m0, i2, OB_ori3);
        
        for i3 = 1:10000
            AK_inver = omega_k .* (p_dm0_inver + (sum(p_d_ori) - p_d_ori(m0,:))');
            rk_inver = log(N * p_dm0_inver .* lambda_n(m0,:)' + sum(p_dm0_inver * beta_n(m0,:))' + 2*N*sqrt(p_dm0_inver) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2) - (g2k_const + Delta_g2_matrix * (p_dm0_inver - p_dm0_fixed));
            z_k_opt = sqrt(rk_inver) ./ AK_inver;
            
            cvx_begin quiet
            variable p_dm0_var(K,1)
            minimize sum( inv_pos(2 * ((log(N * p_dm0_var .* lambda_n(m0,:)' + sum(p_dm0_var * beta_n(m0,:))' + 2*N*sqrt(p_dm0_var) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2)) - (g2k_const + Delta_g2_matrix * (p_dm0_inver - p_dm0_fixed))) - z_k_opt .*  (omega_k .* (p_dm0_var + (sum(p_d_ori) - p_d_ori(m0,:))')) ))
            subject to
            2 * ((log(N * p_dm0_var .* lambda_n(m0,:)' + sum(p_dm0_var * beta_n(m0,:))' + 2*N*sqrt(p_dm0_var) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2)) - (g2k_const + Delta_g2_matrix * (p_dm0_inver - p_dm0_fixed))) - z_k_opt .*  (omega_k .* (p_dm0_var + (sum(p_d_ori) - p_d_ori(m0,:))')) >= 0;
            sum(p_dm0_var) - rho_d * sigma_n <= 0;
            p_dm0_var >= 0;
            N * p_dm0_var .* lambda_n(m0,:)' + 2*N*sqrt(p_dm0_var) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 >= (sum(p_dm0_var * beta_n(m0,:))' + nu_m) .* gamma_k_req;   % SINR requirement
            cvx_end
            
            Ak_var_after = omega_k .* (p_dm0_var + (sum(p_d_ori) - p_d_ori(m0,:))');
            yk1_var = log(N * p_dm0_var .* lambda_n(m0,:)' + sum(p_dm0_var * beta_n(m0,:))' + 2*N*sqrt(p_dm0_var) .* sqrt(lambda_n(m0,:)') .* varphi_m + N * varphi_m.^2 + nu_m) / log(2);
            yk2_var = (g2k_const + Delta_g2_matrix * (p_dm0_inver - p_dm0_fixed));
            
            rk_var_after = yk1_var - yk2_var;
            OB_after3 = sum(Ak_var_after ./ rk_var_after);
            
            if (OB_after3 - OB_ori3)/(OB_ori3) > 0
                break;
            elseif (OB_after3 - OB_ori3)/(OB_ori3) > - acc
                p_dm0_inver = p_dm0_var;
                OB_ori3 = OB_after3;
                break;
            else
                p_dm0_inver = p_dm0_var;
                OB_ori3 = OB_after3;
            end
%             fprintf('m0: %d, i2: %d, i3: %d,inverse Q OB_after3: %f\n', m0, i2, i3, OB_after3);
        end
        
        p_dm0_fixed = p_dm0_inver;
        
        OB_after2 = OB_after3;
        if (OB_after2 - OB_ori2)/(OB_ori2) > - acc
            break;
        end
        
    end
    
    OB_after1 = OB_after2;
%     fprintf('m0: %d, OB_after1: %f ********* \n', m0, OB_after1);
    if OB_after1 < OB_ori1
        OB_ori1 = OB_after1;
        p_d_ori(m0,:) = p_dm0_fixed;
    end
    
    if m0 == M && ori_outer >= OB_ori1
        break;
    end
end

p_d = p_d_ori / sigma_n;

[gamma_k, Rate_k, delay_2, total_delay, energy_2, energy_initial] = systemupdate_pd(rho_k, p_d, GLk, sigma_n, N, Bandwidth, lambda, beta, eta, delay_1, delay_3, energy_1, energy_3);









