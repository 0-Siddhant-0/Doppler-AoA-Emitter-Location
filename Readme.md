# Combined Doppler and AoA Emitter Location

MATLAB implementation and analysis of Combined Doppler and Angle-of-Arrival (AoA) methods for passive emitter localization.

## Mathematical Models

### Doppler Measurement Model

The Doppler frequency shift measurement model is given by:

![Equation](Equations/equation_0_f_t_mathbfx_f_c_.png)

where:
- $\mathbf{x} = [X, Y, Z, f_c]^T$ is the parameter vector
- $f_c$ is the carrier frequency
- $\mathbf{P}(t) = [X_p(t), Y_p(t), Z_p(t)]^T$ is the platform position
- $\mathbf{V}(t) = [V_x(t), V_y(t), V_z(t)]^T$ is the platform velocity
- $\mathbf{X} = [X, Y, Z]^T$ is the emitter position
- $c$ is the speed of light

The observed Doppler measurements include additive noise:

![Equation](Equations/equation_1_tildef_t_i_mathbfx_.png)

where $\nu(t_i) \sim \mathcal{N}(0, \sigma_f^2)$ is Gaussian noise.

### Angle-of-Arrival Measurement Model

The AoA phase measurement model is given by:

![Equation](Equations/equation_2_phi_t_i_phi_0_f.png)

where:
- $\phi_0$ is the phase offset
- $\lambda = \frac{c}{f_c}$ is the wavelength
- $\mathbf{L}(t_i)$ is the scaled baseline vector at time $t_i$

The noisy AoA measurements are modeled as:

![Equation](Equations/equation_3_tildephi_t_i_phi_.png)

where $w_\phi(t_i) \sim \mathcal{N}(0, \sigma_\phi^2)$.

## Jacobian Matrices

### Doppler Jacobian

The Jacobian matrix for Doppler measurements consists of partial derivatives:

$$\mathbf{H}_{Dop} = \begin{bmatrix} \frac{\partial f}{\partial X} & \frac{\partial f}{\partial Y} & \frac{\partial f}{\partial Z} & \frac{\partial f}{\partial f_c} \end{bmatrix}$$

where:

![Equation](Equations/equation_5_fracpartial_fpartial.png)

![Equation](Equations/equation_6_fracpartial_fpartial.png)

![Equation](Equations/equation_7_fracpartial_fpartial.png)

![Equation](Equations/equation_8_fracpartial_fpartial.png)

with:
- $r = \|\mathbf{P} - \mathbf{X}\|$ is the range from platform to emitter
- $\text{Ratio} = \frac{\mathbf{V} \cdot (\mathbf{P} - \mathbf{X})}{r}$ is the projection of velocity onto the line-of-sight vector

### AoA Jacobian

The Jacobian matrix for AoA measurements is:

$$\mathbf{H}_{AoA} = \begin{bmatrix} \frac{\partial \phi}{\partial X} & \frac{\partial \phi}{\partial Y} & \frac{\partial \phi}{\partial Z} & \frac{\partial \phi}{\partial \phi_0} \end{bmatrix}$$

where:

![Equation](Equations/equation_10_fracpartial_phiparti.png)

![Equation](Equations/equation_11_fracpartial_phiparti.png)

![Equation](Equations/equation_12_fracpartial_phiparti.png)

![Equation](Equations/equation_13_fracpartial_phiparti.png)

with:
- $\text{Ratio}_{AoA} = \frac{\mathbf{L} \cdot (\mathbf{P} - \mathbf{X})}{r}$ is the projection of baseline vector onto the line-of-sight

## Combined Nonlinear Least Squares Estimator

The combined estimator minimizes the weighted sum of squared residuals:

![Equation](Equations/equation_14_hatmathbfx_argmin_.png)

where:
- $\mathbf{C}_f = \sigma_f^2 \mathbf{I}$ is the Doppler measurement covariance matrix
- ![C_phi](https://latex.codecogs.com/png.latex?\bg_white%20\mathbf{C}_\phi%20=%20\sigma_\phi^2%20\mathbf{I}) is the AoA measurement covariance matrix

This is solved iteratively using the Gauss-Newton method:

![Equation](Equations/equation_15_mathbfx_k_1_mathbf.png)

where the update step $\delta \mathbf{x}_k$ is:

![Equation](Equations/equation_16_delta_mathbfx_k_le.png)

with:
* The combined Jacobian matrix:
  
  ![H_combined](https://latex.codecogs.com/png.latex?\bg_white%20\mathbf{H}_{combined}%20=%20\begin{bmatrix}%20\mathbf{H}_{Dop}%20\\%20\mathbf{H}_{AoA}%20\end{bmatrix})

* The block-diagonal combined covariance:
  
  ![C_combined](https://latex.codecogs.com/png.latex?\bg_white%20\mathbf{C}_{combined}%20=%20\begin{bmatrix}%20\mathbf{C}_f%20&%20\mathbf{0}%20\\%20\mathbf{0}%20&%20\mathbf{C}_\phi%20\end{bmatrix})

* The residual vector:
  
  ![r_combined](https://latex.codecogs.com/png.latex?\bg_white%20\mathbf{r}_{combined}%20=%20\begin{bmatrix}%20\tilde{\mathbf{f}}%20-%20\mathbf{f}(\mathbf{x}_k)%20\\%20\tilde{\boldsymbol{\phi}}%20-%20\boldsymbol{\phi}(\mathbf{x}_k)%20\end{bmatrix})
  
- $\lambda$ is the regularization parameter for numerical stability

## Cramer-Rao Lower Bound (CRLB) Analysis

The CRLB provides the theoretical lower bound on the covariance matrix of any unbiased estimator:

![Equation](Equations/equation_20_mathbfC_CRLB_mathbfx.png)

where $\mathbf{J}$ is the Fisher Information Matrix (FIM):

### Doppler-only FIM

![Equation](Equations/equation_21_mathbfJ_Dop_frac1s.png)

### AoA-only FIM

![Equation](Equations/equation_22_mathbfJ_AoA_frac1s.png)

### Combined FIM

![Equation](Equations/equation_23_mathbfJ_combined_m.png)

The CRLB trace is used as a scalar measure of estimation uncertainty:

![Equation](Equations/equation_24_textCRLB_texttrace_.png)

## Platform Trajectory Model

The platform trajectory is generated using the weave model, which provides position and velocity as functions of time. The motion is characterized by alternating turns with specified g-force.

The turn angle is modeled as:

![Equation](Equations/equation_25_theta_t_theta_max.png)

where:
- $\theta_{max}$ is the maximum turn angle (30 degrees)
- $v_{horiz}$ is the horizontal velocity
- $r_{turn} = \frac{v_{horiz}^2}{a_{horiz}}$ is the turn radius
- $a_{horiz} = g \cdot 9.81$ m/s² is the horizontal acceleration

The platform velocities are computed as:

![Equation](Equations/equation_26_v_x_v_horiz_cdot_c.png)

![Equation](Equations/equation_27_v_y_v_horiz_cdot_.png)

![Equation](Equations/equation_28_v_z_v_z_a_z_cdot.png)

where $a_z$ alternates between positive and negative values to create vertical oscillation.

## Monte Carlo Simulation

The Monte Carlo simulation generates multiple realizations of noisy measurements and applies the estimator to each:

![Equation](Equations/equation_29_tildef_i_t_j_f_t_.png)

![Equation](Equations/equation_30_tildephi_i_t_j_ph.png)

The RMS error is calculated as:

![Equation](Equations/equation_31_textRMS_x_sqrtfrac.png)

![Equation](Equations/equation_32_textRMS_y_sqrtfrac.png)

![Equation](Equations/equation_33_textRMS_z_sqrtfrac.png)

![Equation](Equations/equation_34_textRMS_3D_sqrttex.png)

## Improvement Ratio Calculation

The improvement ratio of the combined method over Doppler-only is:

![Equation](Equations/equation_35_R_Dop_fractexttrac.png)

The improvement ratio over AoA-only is:

![Equation](Equations/equation_36_R_AoA_fractexttrac.png)

## Results

### Emitter Location Estimation

![Emitter Location Estimation](Results/1.jpg)

### Error Analysis

![Error Ellipse Zoomed View](Results/2.jpg)

### Method Comparison

![Comparison of Location Estimates](Results/3.jpg)
![Zoomed Error Ellipses](Results/4.jpg)

### Trajectory Analysis

![Platform Trajectories](Results/5.jpg)
![Estimation Performance vs Trajectory](Results/6.jpg)
![Improvement Ratio vs Trajectory](Results/7.jpg)

### Measurement Quality Impact

![Measurement Quality Analysis](Results/8.jpg)
![Improvement Ratio Surface](Results/9.jpg)

## Key Findings

1. **Position Estimation Accuracy**:
   - Doppler-only: 3D RMS error = 9.18 m
   - AoA-only: 3D RMS error = 1.29 m
   - Combined: 3D RMS error = 1.05 m

2. **CRLB Analysis**:
   - Doppler-only: $\text{trace}(\mathbf{C}_{CRLB,Dop}) = 8.47 \times 10^1$
   - AoA-only: $\text{trace}(\mathbf{C}_{CRLB,AoA}) = 1.64 \times 10^0$
   - Combined: $\text{trace}(\mathbf{C}_{CRLB,combined}) = 1.13 \times 10^0$

3. **Improvement Ratios**:
   - Over Doppler-only: $R_{Dop} = 75.05$
   - Over AoA-only: $R_{AoA} = 1.45$

4. **Trajectory Impact**:
   - 1.0g: $R_{Dop} = 103.64$, $R_{AoA} = 2.04$
   - 2.0g: $R_{Dop} = 81.36$, $R_{AoA} = 1.58$
   - 3.0g: $R_{Dop} = 75.05$, $R_{AoA} = 1.45$
   - 4.0g: $R_{Dop} = 74.63$, $R_{AoA} = 1.44$

5. **Measurement Quality Impact**:
   - Optimal improvement (2.46×) occurs at $\sigma_f = 3.16$ Hz, $\sigma_\phi = 10.00$ degrees
   - Improvement ratio varies systematically with measurement quality, following a pattern where the combined method provides greatest advantage when compensating for poor quality in one measurement type

## References

[1] N. O'Donoughue, Emitter Detection and Geolocation for Electronic Warfare, ser. Artech House electronic warfare library. Artech House, 2019. [Online]. Available: https://books.google.com/books?id=TbjEDwAAQBAJ

[2] S. Frattasi and F. Rosa, Mobile Positioning and Tracking: From Conventional to Cooperative Techniques, ser. IEEE Press. Wiley, 2017. [Online]. Available: https://books.google.com/books?id=rvImDwAAQBAJ

[3] R. G. Stansfield, "Statistical theory of d.f. fixing," Journal of the Institution of Electrical Engineers - Part IIIA: Radiocommunication, vol. 94, no. 15, pp. 762–770, April 1947.

[4] D. J. Torrieri, "Statistical theory of passive location systems," IEEE Transactions on Aerospace and Electronic Systems, vol. AES-20, pp. 183–198, 1984. [Online]. Available: https://api.semanticscholar.org/CorpusID:8983806

[5] K. Becker, "An efficient method of passive emitter location," IEEE Transactions on Aerospace Electronic Systems, vol. 28, no. 4, pp. 1091–1104, Oct. 1992.

[6] J. Foutz, A. Spanias, and M. Banavar, Narrowband Direction of Arrival Estimation for Antenna Arrays, 01 2008, vol. 3.

[7] M. Fowler, "Eece 522 estimation theory lecture notes," http://www.ws.binghamton.edu/fowler/fowler%20personal%20page/EE522.htm, Jan. 29 2014.

[8] M. Khalaf-Allah, "Emitter location with azimuth and elevation measurements using a single aerial platform for electronic support missions," Sensors, vol. 21, no. 12, 2021. [Online]. Available: https://www.mdpi.com/1424-8220/21/12/3946

[9] Fundamentals of Statistical Signal Processing, Volume 1: Estimation Theory. Pearson Education. [Online]. Available: https://books.google.com/books?id=pDnV5qf1f6IC

[10] E. J. Bailey, "Single platform geolocation of radio frequency emitters," Air Force Institute of Technology, 2015. [Online]. Available: https://scholar.afit.edu/etd/22

[11] M. Fowler, "Analysis of single-platform passive emitter location with terrain data," IEEE Transactions on Aerospace and Electronic Systems, vol. 37, no. 2, pp. 495–507, 2001.

