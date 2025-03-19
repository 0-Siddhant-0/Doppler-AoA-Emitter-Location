function [p_est_iter, H_true] = LBI_loc(phi_accuracy, Plat_Nav_Data, mu_vect, p_true, p_est_init, fo, L)
%
% Usage: p_est_iter=LBI_loc(Plat_Nav_Data,mu_vect,p_true,p_est_init,fo,L);
% 
% Inputs:   phi_accuracy = standard deviation of LBI phase measurements in radians
%           Plat_Nav_Data = Platform Navigation Data  (7xT Matrix)
%                   Each row is a time history of a Nav Quantity at T instants
%                   Order of Rows: Px, Py, Pz, Vx, Vy, Vz, t
%           mu_vect = baseline unit vectors at T instants (3xT Matrix)
%           p_true = True Emitter Location & Frequency (1x4 Vector)
%                   [xe ye ze phi_o] with units [m m m rad]
%           p_est_init = Initial Estimate of Emitter (1x4 Vector)
%                   [xe_init ye_init ze_init phi_o] with units [m m m rad]
%           fo = signal's center frequency  (in Hz)
%           L = baseline length in meters
%
% Outputs: p_est_iter = matrix of estimate iterations
%                   Each column is 4x1 estimate of p_true.'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Put LBI Processing Parameters Here %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 2.998e8;             %% Speed of Light in m/s
lambda = c/fo;
N = 25;                  %% Max. # of Iterations

L_scaled = L*mu_vect;   %% Gives unit baseline vector length of L
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get Nav Data %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Px = Plat_Nav_Data(1,:);
Py = Plat_Nav_Data(2,:);
Pz = Plat_Nav_Data(3,:);
Vx = Plat_Nav_Data(4,:);
Vy = Plat_Nav_Data(5,:);
Vz = Plat_Nav_Data(6,:);
t = Plat_Nav_Data(7,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get True Emitter Location/Frequency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xe = p_true(1);    %% meters
ye = p_true(2);    %% meters
ze = p_true(3);    %% meters
phi_o = p_true(4);    %% rad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate Perfect LBI Phase Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = sqrt((Px-xe).^2 + (Py-ye).^2 + (Pz-ze).^2);
phi_nf = phi_o - (2*pi/lambda)*(mu_vect(1,:).*(Px-xe) + mu_vect(2,:).*(Py-ye) + mu_vect(3,:).*(Pz-ze))./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make Noise-Corrupted LBI Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_LBI = phi_accuracy.^2;  % Convert Std. Dev. into Variance
phi = phi_nf + sqrt(var_LBI)*randn(size(phi_nf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now Start the Nonlinear LS Estimation Processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pre-Allocate space for the iterated estimates
p_est_iter = zeros(4, N+1);  %% One column for each iteration


%%% Get Initial Estimate and put into Estimate Iteration Matrix%%%
xe_hat = p_est_init(1);
ye_hat = p_est_init(2);
ze_hat = p_est_init(3);
phi_o_hat = p_est_init(4);

p_est_iter(:, 1) = p_est_init.';

%%% Iterations:

for n = 1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%% Generate "Predicted" measurements based on Current Estimate:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    R_hat = sqrt((Px-xe_hat).^2 + (Py-ye_hat).^2 + (Pz-ze_hat).^2);
    phi_hat = phi_o_hat - (2*pi/lambda)*(mu_vect(1,:).*(Px-xe_hat) + mu_vect(2,:).*(Py-ye_hat) + mu_vect(3,:).*(Pz-ze_hat))./R_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Generate the LS residual:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    res_phi = phi - phi_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute Jacobian Matrix needed for LS Update Computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	
    diff_X = Px-xe_hat;
    diff_Y = Py-ye_hat;
    diff_Z = Pz-ze_hat;
    
    Numer = (L_scaled(1,:).*diff_X) + (L_scaled(2,:).*diff_Y) + (L_scaled(3,:).*diff_Z);
    Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
    Ratio = Numer./Denom;
    
    H1 = (2*pi/lambda)*(L_scaled(1,:).*Denom - (Ratio).*diff_X)./(Denom.^2);
    H2 = (2*pi/lambda)*(L_scaled(2,:).*Denom - (Ratio).*diff_Y)./(Denom.^2);
    H3 = (2*pi/lambda)*(L_scaled(3,:).*Denom - (Ratio).*diff_Z)./(Denom.^2);
    H4 = ones(size(Px));
    
    H = [H1; H2; H3; H4]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute Estimate-Update for this Iteration:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add small regularization to avoid singularity issues
    lambda_reg = 1e-8;
    H_reg = (H')*H + lambda_reg*eye(size(H, 2));
    
    del_p_est = inv(H_reg)*(H')*(res_phi');  
    %% NOTE: This assumes that measurement covariance = (sigma^2)*Identity
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Update the Estimate: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Update and stuff into storage matrix
    p_est_iter(:, n+1) = p_est_iter(:, n) + del_p_est;  
    
    %% Stuff new estimate back into names used in iteration calculations:
    xe_hat = p_est_iter(1, n+1);
    ye_hat = p_est_iter(2, n+1);
    ze_hat = p_est_iter(3, n+1);
    phi_o_hat = p_est_iter(4, n+1);

end   %%% End "LS Iterations" loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute Jacobian Matrix needed for CRLB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	
diff_X = Px-xe;
diff_Y = Py-ye;
diff_Z = Pz-ze;

Numer = (L_scaled(1,:).*diff_X) + (L_scaled(2,:).*diff_Y) + (L_scaled(3,:).*diff_Z);
Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
Ratio = Numer./Denom;

H1 = (2*pi/lambda)*(L_scaled(1,:).*Denom - (Ratio).*diff_X)./(Denom.^2);
H2 = (2*pi/lambda)*(L_scaled(2,:).*Denom - (Ratio).*diff_Y)./(Denom.^2);
H3 = (2*pi/lambda)*(L_scaled(3,:).*Denom - (Ratio).*diff_Z)./(Denom.^2);
H4 = ones(size(Px));

H_true = [H1; H2; H3; H4]';
end