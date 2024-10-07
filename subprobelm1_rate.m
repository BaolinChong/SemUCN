   %% subproblem 1
    % 1.1 rho_k optimization
    rho_k_start = rho_k;
    
    weight_y1 = 0.55;
    coe1_1 = 1.1;
    coe1_2 = 0.9;
    
    y1_1_Ori = weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k - C2).^C3 + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3);
    y2_1_Ori = C4 .* pow_p(rho_k, -C5);
    OB_Ori = kappa * sum(y1_1_Ori .* f.^2) + sum(GLk .* rho_k ./ Rate_k .* (sum(p_d)' * sigma_n)) + kappa * sum(y2_1_Ori .* g.^2);
    
%     fprintf('首次OB：%f\n',OB_Ori)
    
    y1_1_Iter = weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k_start - C2).^C3 - C1 .* C3 .* (rho_k_start - C2).^(C3 - 1) .* (rho_k_start - rho_k) + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3);    
    
    for i_sub1 = 1:10000        
        cvx_begin quiet
        variable rho_k_var(K, 1);
        minimize kappa * sum((weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k_start - C2).^C3 - C1 .* C3 .* (rho_k_start - C2).^(C3 - 1) .* (rho_k_var - rho_k) + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3)) .* f.^2) + sum(GLk .* rho_k_var ./ Rate_k .* (sum(p_d)' * sigma_n)) + kappa * sum((C4 .* pow_p(rho_k_var, -C5)) .* g.^2);
        subject to
        (weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k_start - C2).^C3 - C1 .* C3 .* (rho_k_start - C2).^(C3 - 1) .* (rho_k_var - rho_k) + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3)) ./ f + GLk .* rho_k_var ./ Rate_k + (C4 .* pow_p(rho_k_var, -C5)) ./ g <= Delay_req;
        rho_req <= rho_k_var <= ones(K,1); % semantic similarity requrement
        cvx_end
        
        if all((weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k_start - C2).^C3 - C1 .* C3 .* (rho_k_start - C2).^(C3 - 1) .* (rho_k_var - rho_k) + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3)) ./ f + GLk .* rho_k_var ./ Rate_k + (C4 .* pow_p(rho_k_var, -C5)) ./ g <= Delay_req) == 0
            break;
        end
        
        if isnan(rho_k_var(1))
            break;
        end
        
        rho_k_start = rho_k_var;
        
        y1_1_after = weight_y1 * Lk.^(coe1_1) + (1 - weight_y1) * Lk.^(coe1_2) - C1 .* (rho_k_var - C2).^C3 + max(C1 .* (0 - C2).^C3, C1 .* (1 - C2).^C3);
        y2_1_after = C4 .* pow_p(rho_k_var, -C5);
        OB_after = kappa * sum(y1_1_after .* f.^2) + sum(GLk .* rho_k_var ./ Rate_k .* (sum(p_d)' * sigma_n)) + kappa * sum(y2_1_after .* g.^2);
        
%         rho_k = rho_k_var;
        if (OB_Ori-OB_after) / OB_Ori < 0.01
            rho_k = rho_k_var;
            break;
        else
            OB_Ori = OB_after;
        end
    end
    
    % 1.2 rho_k and system update
    [Sk, y1, y2, delay_1, delay_2, delay_3, total_delay, energy_1, energy_2, energy_3, energy_initial] = systemupdate_rhok(rho_k, f, g, p_d, Rate_k, Lk, GLk, kappa, sigma_n, C1, C2, C3, C4, C5);
  