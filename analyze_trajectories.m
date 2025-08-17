function analyze_trajectories()
% This function analyzes the impact of different platform trajectories on the
% estimation performance for Doppler-only, AoA-only, and combined approaches.


% fixed parameters
dop_std_dev = 1;
aoa_std_dev = 0.5;
alt_kft = 10;
vel = 200;
p_true = [5000, 5000, 0, 1e9];
p_est_init = [5100, 4900, 0, 1.000001e9];
fo = 1e9;
L = 2;
T = 60;
del_T = 1;

% different trajectory parameters
g_values = [1, 2, 3, 4];

% store results
crlb_traces = zeros(3, length(g_values));
crlb_matrices = cell(3, length(g_values));

figure;
hold on;

% colors for different methods
colors = {'b', 'g', 'r'};

% loop through different g values
for i = 1:length(g_values)
    g = g_values(i);
    
    [Px, Py, Pz, Vx, Vy, Vz, mu_vect] = weave(g, T, del_T, alt_kft, vel);
    t = 0:del_T:T;
    Plat_Nav_Data = [Px; Py; Pz; Vx; Vy; Vz; t];
    
    aoa_std_dev_rad = aoa_std_dev * pi/180;
    
    [~, H_true_dop] = doppler_loc(dop_std_dev, Plat_Nav_Data, p_true, p_est_init);
    J_dop = (1/(dop_std_dev^2)) * (H_true_dop' * H_true_dop);
    CRLB_dop = inv(J_dop);
    crlb_traces(1, i) = trace(CRLB_dop(1:3, 1:3));
    crlb_matrices{1, i} = CRLB_dop;
    
    p_true_aoa = p_true;
    p_true_aoa(4) = 0;
    p_est_init_aoa = p_est_init;
    p_est_init_aoa(4) = 0;
    
    [~, H_true_aoa] = LBI_loc(aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true_aoa, p_est_init_aoa, fo, L);
    J_aoa = (1/(aoa_std_dev_rad^2)) * (H_true_aoa' * H_true_aoa);
    CRLB_aoa = inv(J_aoa);
    crlb_traces(2, i) = trace(CRLB_aoa(1:3, 1:3));
    crlb_matrices{2, i} = CRLB_aoa;
    
    [~, H_true_combined] = doppler_aoa_loc(dop_std_dev, aoa_std_dev_rad, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
    C_dop = (dop_std_dev^2) * eye(length(Px));
    C_aoa = (aoa_std_dev_rad^2) * eye(length(Px));
    C_combined = blkdiag(C_dop, C_aoa);
    J_combined = H_true_combined' * inv(C_combined) * H_true_combined;
    CRLB_combined = inv(J_combined);
    crlb_traces(3, i) = trace(CRLB_combined(1:3, 1:3));
    crlb_matrices{3, i} = CRLB_combined;
    
    subplot(4, 2, (i-1)*2 + 1);
    plot(Px, Py, 'k-', p_true(1), p_true(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    ylabel(sprintf('Y Position (m)\n(g=%.1f)', g_values(i)));
    grid on;
    axis equal;
    if i == 1
        title('Full Trajectory');
        legend('Platform Path', 'True Location', 'Location', 'NorthWest');
    end
    if i == length(g_values)
        xlabel('X Position (m)');
    end
    
    subplot(4, 2, (i-1)*2 + 2);
    plot(p_true(1), p_true(2), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    hold on;
    
    for j = 1:3
        ellipse(crlb_matrices{j, i}(1:2, 1:2), 6, colors{j}, p_true(1:2));
    end
    
    zoom_level = 15;
    xlim([p_true(1) - zoom_level, p_true(1) + zoom_level]);
    ylim([p_true(2) - zoom_level, p_true(2) + zoom_level]);
    grid on;
    axis equal;
    
    if i == 1
        title('Zoomed CRLB Ellipses');
        legend('True', 'Dop', 'AoA', 'Comb');
    end
    if i == length(g_values)
        xlabel('X Position (m)');
    end
end
set(gcf, 'Renderer', 'painters');
drawnow;

% plot CRLB trace comparison
figure;
semilogy(g_values, crlb_traces(1, :), 'bo-', g_values, crlb_traces(2, :), 'go-', g_values, crlb_traces(3, :), 'ro-');
title('Impact of Platform Trajectory on Estimation Performance');
xlabel('Platform Turn Rate (g)');
ylabel('Trace of CRLB Matrix (log scale)');
legend('Doppler-only', 'AoA-only', 'Combined Doppler-AoA');
grid on;
set(gcf, 'Renderer', 'painters');
drawnow;

% calculate improvement ratio
improvement_over_doppler = crlb_traces(1, :) ./ crlb_traces(3, :);
improvement_over_aoa = crlb_traces(2, :) ./ crlb_traces(3, :);
figure;

semilogy(g_values, improvement_over_doppler, 'b-o', g_values, improvement_over_aoa, 'g-o');
title('Improvement Ratio of Combined Method (Log Scale)');
xlabel('Platform Turn Rate (g)');
ylabel('Improvement Ratio (Log Scale)');
legend('Improvement over Doppler-only', 'Improvement over AoA-only');
grid on;
set(gcf, 'Renderer', 'painters');
drawnow;

% display summary
fprintf('Trajectory Analysis Summary:\n');
fprintf('----------------------------\n');
fprintf('Turn Rate (g) | Doppler CRLB | AoA CRLB | Combined CRLB | Improve over Doppler | Improve over AoA\n');
for i = 1:length(g_values)
    fprintf('%.1f         | %.2e    | %.2e | %.2e    | %.2f times          | %.2f times\n', ...
        g_values(i), crlb_traces(1, i), crlb_traces(2, i), crlb_traces(3, i), ...
        improvement_over_doppler(i), improvement_over_aoa(i));
end

end