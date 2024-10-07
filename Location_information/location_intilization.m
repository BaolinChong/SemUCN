rng(9);
% initialization
average_times = 200; % iteration times 
M = 100; % maximum number of APs
K = 50;  % maximum number of UEs


location_APs = rand(average_times,M,2)*1000; % 2-D location of APs
location_UEs = rand(average_times,K,2)*1000; % 2-D location of UEs

a_shadow = rand(average_times, M);
b_shadow = rand(average_times, K);
delta_shadow = rand(average_times,1);

save('location_information_rng10', 'location_APs', 'location_UEs','a_shadow','b_shadow','delta_shadow')