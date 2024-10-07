system_init;

cvx_solver SDPT3
% joint optimization algorithm
outer_iter_times = 1;
max_iter = 1000;
max_iter_inner = 1000;
tol = 1 * 10^(-1);

outer_obj_now = energy_initial;
outer_obj_old = outer_obj_now;
% fprintf('Begining: the total energy consumption is: %f\n', outer_obj_now);

subproblem3_power;
% fprintf('3.1 iteration is: %f\n', energy_initial);

while  outer_iter_times < max_iter && ( (outer_obj_old - outer_obj_now) / (outer_obj_now) > tol || outer_iter_times == 1)
    outer_obj_old = outer_obj_now;
%     subproblem2_compuation;

    %% subproblem 1
    % 2.1 rho_k optimization
    subprobelm1_rate;
%     fprintf('3.2 iteration is: %f\n', energy_initial)
    
   %% subproblem 2
    % 3.1 f, g optimization
    subproblem2_compuation;
%     fprintf('3.3 iteration is: %f\n', energy_initial)
    
    %%
    outer_obj_now = energy_initial;
%     fprintf('The total energy consumption at %d-th iteration is: %f\n', outer_iter_times, outer_obj_now);
    outer_iter_times = outer_iter_times + 1;
end
% energy_initial
