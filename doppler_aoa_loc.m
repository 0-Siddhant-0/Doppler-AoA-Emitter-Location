function [p_est_iter, H_true_combined] = doppler_aoa_loc(freq_accuracy, phi_accuracy, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L)
%
% Usage: [p_est_iter, H_true_combined] = doppler_aoa_loc(freq_accuracy, phi_accuracy, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L);
% 
% Inputs:   freq_accuracy = standard deviation of frequency measurements in Hz
%           phi_accuracy = standard deviation of angle measurements in radians
%           Plat_Nav_Data = Platform Navigation Data  (7xT Matrix)
%                   Each row is a time history of a Nav Quantity at T instants
%                   Order of Rows: Px, Py, Pz, Vx, Vy, Vz, t
%           mu_vect = baseline unit vectors at T instants (3xT Matrix)
%           p_true = True Emitter Location & Frequency (1x4 Vector)
%                   [xe ye ze fo] with units [m m m Hz]
%           p_est_init = Initial Estimate of Emitter (1x4 Vector)
%                   [xe_init ye_init ze_init fo_init] with units [m m m Hz]
%           fo = signal's center frequency  (in Hz)
%           L = baseline length in meters
%
% Outputs: p_est_iter = matrix of estimate iterations
%                   Each column is 4x1 estimate of p_true.'
%          H_true_combined = True Jacobian computed at true emitter location for CRLB computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Combined Doppler+AoA Processing Parameters Here %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
c = 2.998e8;             %% Speed of Light in m/s
N = 25;                  %% Max. # of Iterations
lambda = c/fo;           %% Wavelength
L_scaled = L*mu_vect;    %% Gives unit baseline vector length of L
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Nav Data %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Px = Plat_Nav_Data(1,:);
Py = Plat_Nav_Data(2,:);
Pz = Plat_Nav_Data(3,:);
Vx = Plat_Nav_Data(4,:);
Vy = Plat_Nav_Data(5,:);
Vz = Plat_Nav_Data(6,:);
t = Plat_Nav_Data(7,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% True Emitter Location/Frequency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xe = p_true(1);    %% meters
ye = p_true(2);    %% meters
ze = p_true(3);    %% meters
true_fo = p_true(4);    %% Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate Perfect Doppler Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = sqrt((Px-xe).^2 + (Py-ye).^2 + (Pz-ze).^2);
f_nf = true_fo - (true_fo/c)*(Vx.*(Px-xe) + Vy.*(Py-ye) + Vz.*(Pz-ze))./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make Noise-Corrupted Doppler Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_dop = freq_accuracy.^2;  % Convert Std. Dev. into Variance
f = f_nf + sqrt(var_dop)*randn(size(f_nf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate Perfect AoA Phase Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_nf = 0 - (2*pi/lambda)*(mu_vect(1,:).*(Px-xe) + mu_vect(2,:).*(Py-ye) + mu_vect(3,:).*(Pz-ze))./R;
% Note: Using a phase offset of 0 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make Noise-Corrupted AoA Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_aoa = phi_accuracy.^2;  % Convert Std. Dev. into Variance
phi = phi_nf + sqrt(var_aoa)*randn(size(phi_nf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now Start the Nonlinear LS Estimation Processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pre-Allocate space for the iterated estimates
p_est_iter = zeros(4, N+1);  %% One column for each iteration

%%% Get Initial Estimate and put into Estimate Iteration Matrix%%%
xe_hat = p_est_init(1);
ye_hat = p_est_init(2);
ze_hat = p_est_init(3);
fo_hat = p_est_init(4);

p_est_iter(:, 1) = p_est_init.';

%%% Iterations:
for n = 1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%% Generate "Predicted" measurements based on Current Estimate:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    R_hat = sqrt((Px-xe_hat).^2 + (Py-ye_hat).^2 + (Pz-ze_hat).^2);
    
    % Doppler predictions
    f_hat = fo_hat - (fo_hat/c)*(Vx.*(Px-xe_hat) + Vy.*(Py-ye_hat) + Vz.*(Pz-ze_hat))./R_hat;
    
    % AoA predictions
    phi_hat = 0 - (2*pi/lambda)*(mu_vect(1,:).*(Px-xe_hat) + mu_vect(2,:).*(Py-ye_hat) + mu_vect(3,:).*(Pz-ze_hat))./R_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Generate the LS residual:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    res_f = f - f_hat;
    res_phi = phi - phi_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute Jacobian Matrix needed for LS Update Computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% First need some preliminary things to make things easier: 
    diff_X = Px-xe_hat;  
    diff_Y = Py-ye_hat;   
    diff_Z = Pz-ze_hat;
        
    %% For Doppler Jacobian
    Numer = (Vx.*diff_X) + (Vy.*diff_Y) + (Vz.*diff_Z);
    Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
    Ratio = Numer./Denom;
    
    H_dop_1 = (fo_hat/c)*(Vx.*Denom - Ratio.*diff_X)./(Denom.^2);  % df/dxe
    H_dop_2 = (fo_hat/c)*(Vy.*Denom - Ratio.*diff_Y)./(Denom.^2);  % df/dye
    H_dop_3 = (fo_hat/c)*(Vz.*Denom - Ratio.*diff_Z)./(Denom.^2);  % df/dze
    H_dop_4 = 1 - (1/c)*Ratio;                                    % df/dfo
    
    H_dop = [H_dop_1; H_dop_2; H_dop_3; H_dop_4]';
    
    %% For AoA Jacobian
    Numer_aoa = (L_scaled(1,:).*diff_X) + (L_scaled(2,:).*diff_Y) + (L_scaled(3,:).*diff_Z);
    Ratio_aoa = Numer_aoa./Denom;
    
    H_aoa_1 = (2*pi/lambda)*(L_scaled(1,:).*Denom - (Ratio_aoa).*diff_X)./(Denom.^2);  % dphi/dxe
    H_aoa_2 = (2*pi/lambda)*(L_scaled(2,:).*Denom - (Ratio_aoa).*diff_Y)./(Denom.^2);  % dphi/dye
    H_aoa_3 = (2*pi/lambda)*(L_scaled(3,:).*Denom - (Ratio_aoa).*diff_Z)./(Denom.^2);  % dphi/dze
    H_aoa_4 = zeros(size(Px));  % dphi/dfo - phase doesn't depend on frequency directly
    
    H_aoa = [H_aoa_1; H_aoa_2; H_aoa_3; H_aoa_4]';
    
    %% Set up weighted least squares problem
    
    % Weight the measurements by their standard deviations
    W_dop = eye(length(f)) / var_dop;
    W_aoa = eye(length(phi)) / var_aoa;
    
    % Compute the normal equations individually
    A_dop = H_dop' * W_dop * H_dop;
    b_dop = H_dop' * W_dop * res_f';
    
    A_aoa = H_aoa' * W_aoa * H_aoa;
    b_aoa = H_aoa' * W_aoa * res_phi';
    
    % Combine the normal equations
    A_combined = A_dop + A_aoa;
    b_combined = b_dop + b_aoa;
    
    % Add strong regularization to stabilize the solution
    lambda_reg = 1;  % Strong regularization parameter
    A_reg = A_combined + lambda_reg * eye(size(A_combined));
    
    % Solve the system using SVD for better numerical stability
    [U, S, V] = svd(A_reg);
    s = diag(S);
    tol = max(size(A_reg)) * eps(max(s));
    r = sum(s > tol);
    s_inv = zeros(size(s));
    s_inv(1:r) = 1./s(1:r);
    
    % Compute the update
    del_p_est = V * diag(s_inv) * U' * b_combined;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Update the Estimate, with damping to prevent divergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    damping = 0.5;  % Damping factor to prevent overshooting
    p_est_iter(:, n+1) = p_est_iter(:, n) + damping * del_p_est;
    
    %% Stuff new estimate back into names used in iteration calculations:
    xe_hat = p_est_iter(1, n+1);
    ye_hat = p_est_iter(2, n+1);
    ze_hat = p_est_iter(3, n+1);
    fo_hat = p_est_iter(4, n+1);
    
    % Check for convergence
    if norm(del_p_est) < 1e-6
        p_est_iter = p_est_iter(:, 1:n+1);  % Truncate unused iterations
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute Jacobian Matrix needed for CRLB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the true emitter location
diff_X = Px-xe;
diff_Y = Py-ye;
diff_Z = Pz-ze;

% For Doppler CRLB
Numer = (Vx.*diff_X) + (Vy.*diff_Y) + (Vz.*diff_Z);
Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
Ratio = Numer./Denom;

H_dop_true_1 = (true_fo/c)*(Vx.*Denom - Ratio.*diff_X)./(Denom.^2);  % df/dxe
H_dop_true_2 = (true_fo/c)*(Vy.*Denom - Ratio.*diff_Y)./(Denom.^2);  % df/dye
H_dop_true_3 = (true_fo/c)*(Vz.*Denom - Ratio.*diff_Z)./(Denom.^2);  % df/dze
H_dop_true_4 = 1 - (1/c)*Ratio;                                     % df/dfo

H_dop_true = [H_dop_true_1; H_dop_true_2; H_dop_true_3; H_dop_true_4]';

% For AoA CRLB
Numer_aoa = (L_scaled(1,:).*diff_X) + (L_scaled(2,:).*diff_Y) + (L_scaled(3,:).*diff_Z);
Ratio_aoa = Numer_aoa./Denom;

H_aoa_true_1 = (2*pi/lambda)*(L_scaled(1,:).*Denom - (Ratio_aoa).*diff_X)./(Denom.^2);  % dphi/dxe
H_aoa_true_2 = (2*pi/lambda)*(L_scaled(2,:).*Denom - (Ratio_aoa).*diff_Y)./(Denom.^2);  % dphi/dye
H_aoa_true_3 = (2*pi/lambda)*(L_scaled(3,:).*Denom - (Ratio_aoa).*diff_Z)./(Denom.^2);  % dphi/dze
H_aoa_true_4 = zeros(size(Px));                                                      % dphi/dfo

H_aoa_true = [H_aoa_true_1; H_aoa_true_2; H_aoa_true_3; H_aoa_true_4]';

% Combined Jacobian for CRLB
H_true_combined = [H_dop_true; H_aoa_true];

end
