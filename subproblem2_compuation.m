    %% subproblem 2
    % 2.1 f, g optimization
    cvx_begin quiet
    % cvx_solver
    variables f_val(K, 1) g_val(K, 1);
    minimize kappa * 1e20 * sum(y1 .* f_val.^2) + kappa * 1e20 * sum(y2 .* g_val.^2);
    subject to
    sum(f_val) <= f_max * 1e-9;
    f_val >= 0;
    0 <= g_val <= ones(K, 1) * g_max * 1e-9;
    y1 * 1e-9 .* inv_pos(f_val) + GLk .* rho_k ./ Rate_k + y2 * 1e-9 .* inv_pos(g_val) <= Delay_req;
    cvx_end
    
    % 2.2 f, g, and system update
    f = f_val * 1e9; g = g_val * 1e9;
    [delay_1, delay_3, total_delay, energy_1, energy_3, energy_initial] = systemupdate_fg(f, g, y1, y2, delay_2, energy_2, kappa);
%     fprintf('The total energy consumption at %d-th iteration after computation capacity control is: %f\n', outer_iter_times, energy_initial);