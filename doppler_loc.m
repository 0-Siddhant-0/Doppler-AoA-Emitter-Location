function [p_est_iter, H_true] = doppler_loc(freq_accuracy, Plat_Nav_Data, p_true, p_est_init)
%
% Usage: p_est_iter=doppler_loc(Plat_Nav_Data,p_true,p_est_init);
% 
% Inputs:   freq_accuracy = standard deviation of frequency measurements in Hz
%           Plat_Nav_Data = Platform Navigation Data  (7xT Matrix)
%                   Each row is a time history of a Nav Quantity at T instants
%                   Order of Rows: Px, Py, Pz, Vx, Vy, Vz, t
%           p_true = True Emitter Location & Frequency (1x4 Vector)
%                   [xe ye ze fo] with units [m m m Hz]
%           p_est_init = Initial Estimate of Emitter (1x4 Vector)
%                   [xe_init ye_init ze_init fo_init] with units [m m m Hz]
%
% Outputs: p_est_iter = matrix of estimate iterations
%                   Each column is 4x1 estimate of p_true.'
%          H_true = True Jacobian computed at true emitter location for use in CRLB computation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Put Doppler Processing Parameters Here %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
c = 2.998e8;             %% Speed of Light in m/s
N = 25;                  %% Max. # of Iterations
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
fo = p_true(4);    %% Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate Perfect Doppler Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = sqrt((Px-xe).^2 + (Py-ye).^2 + (Pz-ze).^2);
f_nf = fo - (fo/c)*(Vx.*(Px-xe) + Vy.*(Py-ye) + Vz.*(Pz-ze))./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make Noise-Corrupted Doppler Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_dop = freq_accuracy.^2;  % Convert Std. Dev. into Variance
f = f_nf + sqrt(var_dop)*randn(size(f_nf));

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
fo_hat = p_est_init(4);

p_est_iter(:, 1) = p_est_init.';

%%% Iterations:

for n = 1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%% Generate "Predicted" measurements based on Current Estimate:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_hat = sqrt((Px-xe_hat).^2 + (Py-ye_hat).^2 + (Pz-ze_hat).^2);
    f_hat = fo_hat - (fo_hat/c)*(Vx.*(Px-xe_hat) + Vy.*(Py-ye_hat) + Vz.*(Pz-ze_hat))./R_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Generate the LS residual:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    res_f = f - f_hat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute Jacobian Matrix needed for LS Update Computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% First need some preliminary things to make things easier: 
    diff_X = Px-xe_hat;  
    diff_Y = Py-ye_hat;   
    diff_Z = Pz-ze_hat;
    
    Numer = (Vx.*diff_X) + (Vy.*diff_Y) + (Vz.*diff_Z);
    Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
    Ratio = Numer./Denom;
    
    %% Then Compute Columns of Jacobian Matrix:
    H1 = (fo_hat/c)*(Vx.*Denom - Ratio.*diff_X)./(Denom.^2);    %% df/dxe
    H2 = (fo_hat/c)*(Vy.*Denom - Ratio.*diff_Y)./(Denom.^2);    %% df/dye
    H3 = (fo_hat/c)*(Vz.*Denom - Ratio.*diff_Z)./(Denom.^2);    %% df/dze
    H4 = 1 - (1/c)*Ratio;                                     %% df/dfo
    
    H = [H1; H2; H3; H4]';   %% Put columns into matrix...  Jacobian Matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute Estimate-Update for this Iteration:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add small regularization to avoid singularity issues
    lambda_reg = 1e-8;
    H_reg = (H')*H + lambda_reg*eye(size(H, 2));
    
    del_p_est = inv(H_reg)*(H')*(res_f'); 
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
    fo_hat = p_est_iter(4, n+1);

end   %%% End "LS Iterations" loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute Jacobian Matrix needed for CRLB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	
diff_X = Px-xe;  
diff_Y = Py-ye;   
diff_Z = Pz-ze;

Numer = (Vx.*diff_X) + (Vy.*diff_Y) + (Vz.*diff_Z);
Denom = sqrt(diff_X.^2 + diff_Y.^2 + diff_Z.^2);
Ratio = Numer./Denom;

H1 = (fo/c)*(Vx.*Denom - Ratio.*diff_X)./(Denom.^2);    %% df/dxe
H2 = (fo/c)*(Vy.*Denom - Ratio.*diff_Y)./(Denom.^2);    %% df/dye
H3 = (fo/c)*(Vz.*Denom - Ratio.*diff_Z)./(Denom.^2);    %% df/dze
H4 = 1 - (1/c)*Ratio;                                  %% df/dfo

H_true = [H1; H2; H3; H4]';   %% Put columns into matrix...  Jacobian Matrix
end