% main_doppler_aoa.m
% This script demonstrates the complete analysis of the Combined Doppler and AoA
% emitter location method, including comparisons with individual methods,
% trajectory analysis, and measurement quality analysis.
%
% Author: Siddhant Chaurasia

% Clear workspace and command window
clear all;
close all;
clc;

fprintf('=====================================================\n');
fprintf('  COMBINED DOPPLER AND AoA EMITTER LOCATION ANALYSIS\n');
fprintf('=====================================================\n\n');

% Define standard parameters
g = 3;                  % 3g turn
T = 60;                 % 60 seconds observation
del_T = 1;              % 1 second between measurements
dop_std_dev = 1;        % 1 Hz Doppler measurement std dev
aoa_std_dev = 0.5;      % 0.5 degrees AoA measurement std dev
alt_kft = 10;           % 10,000 ft altitude
vel = 200;              % 200 m/s platform speed
p_true = [5000, 5000, 0, 1e9];  % True emitter location and frequency
p_est_init = [5100, 4900, 0, 1.000001e9];  % Initial estimate
fo = 1e9;               % 1 GHz carrier frequency
L = 2;                  % 2 meter baseline

% Part 1: Basic Combined Doppler-AoA Estimation
fprintf('Running Combined Doppler-AoA Estimation...\n');
fprintf('----------------------------------------\n');
run_doppler_aoa(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L);
fprintf('\nPress any key to continue...\n');
pause;

% Part 2: Comparative Analysis
fprintf('\nRunning Comparative Analysis...\n');
fprintf('------------------------------\n');
compare_methods(g, T, del_T, dop_std_dev, aoa_std_dev, alt_kft, vel, p_true, p_est_init, fo, L);
fprintf('\nPress any key to continue...\n');
pause;

% Part 3: Trajectory Analysis
fprintf('\nAnalyzing Impact of Platform Trajectory...\n');
fprintf('----------------------------------------\n');
analyze_trajectories();
fprintf('\nPress any key to continue...\n');
pause;

% Part 4: Measurement Quality Analysis
fprintf('\nAnalyzing Impact of Measurement Quality...\n');
fprintf('-----------------------------------------\n');
analyze_measurement_quality();

% Print conclusion
fprintf('\n=====================================================\n');
fprintf('                    CONCLUSION\n');
fprintf('=====================================================\n');
fprintf('The analysis demonstrates that combining Doppler and AoA\n');
fprintf('measurements provides improved emitter location accuracy\n');
fprintf('compared to using either method alone. The combined method\n');
fprintf('theoretically shows a significant improvement over the Doppler-only\n');
fprintf('approach, and a modest improvement over the AoA-only approach.\n');
fprintf('\nThese results support the research hypothesis that there is a\n');
fprintf('synergistic effect between Doppler and AoA-based emitter location.\n');