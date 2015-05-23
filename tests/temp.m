clear
clc

load mbefore
bc_normal_velocity = compute_Womersley_profiles(pptp(:, :, 1:2), bctp(:, 1:2), time, T, mu, input_file_list{i});

clear
clc
load tbefore
bc_normal_velocity = compute_Womersley_profiles(vx, bctp(:, 1:2), t, T, mu, '/home/jo/ptc/tests/fig18.csv');