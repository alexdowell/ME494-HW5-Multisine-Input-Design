# ME494 HW5: Multisine Input Design  

## Description  
This repository contains the fifth homework assignment for **ME494**, focusing on system identification through optimized input design. The homework explores the use of multisine and frequency sweep signals for estimating the dynamics of a pendulum system without a controller. The repository includes MATLAB scripts, datasets, and a PDF containing problem descriptions.  

## Files Included  

### **Part 1: Multisine and Frequency Sweep Input Generation**  
- **File:** SID_HW5.m  
- **Topics Covered:**  
  - Generating frequency sweep and multisine signals  
  - Frequency domain analysis using FFT  
  - Evaluating signal characteristics and peak factors  

### **Part 2: Pendulum Simulation and System Identification**  
- **File:** pend.m  
- **Topics Covered:**  
  - Simulating pendulum motion with applied inputs  
  - Recording system response for identification  
  - Least squares parameter estimation  
  - Confidence interval calculation for estimated parameters  

### **Part 3: Residual Analysis and Normalized Regressors**  
- **File:** cvec.m  
- **File:** peakfactor.m  
- **File:** pf_cost.m  
- **Topics Covered:**  
  - Residual analysis for model performance evaluation  
  - Normalization of regressor variables  
  - Coefficient significance and relative importance assessment  

### **Part 4: Input Signal Processing and Optimization**  
- **File:** mkmsswp.m  
- **File:** mksswp.m  
- **Topics Covered:**  
  - Generating multisine and frequency sweep signals  
  - Optimizing input signals for system identification  
  - Phase optimization to minimize peak factors  

### **Homework Assignment Document**  
- **File:** SID_HW5_2022.pdf  
- **Contents:**  
  - Problem descriptions and equations  
  - MATLAB implementation steps  
  - Expected results and discussion points  

## Installation  
Ensure MATLAB is installed before running the scripts. No additional toolboxes are required.  

## Usage  

### **Generating and Analyzing Input Signals**  
1. Open MATLAB.  
2. Run the script:  
   ```SID_HW5```  
3. View generated multisine and frequency sweep signals.  
4. Compare frequency responses using FFT analysis.  

### **Running the Pendulum Simulation and System Identification**  
1. Open MATLAB.  
2. Run the script:  
   ```SID_HW5```  
3. Analyze pendulum motion response and estimated parameters.  
4. Compute confidence intervals and residual plots.  

### **Performing Residual Analysis and Regression Normalization**  
1. Open MATLAB.  
2. Run the script:  
   ```SID_HW5```  
3. Observe residuals over time and normalized regressor coefficients.  

## Example Output  

- **Input Signal Comparison**  
  - Frequency sweep vs. multisine signal comparison in time and frequency domains  
  - Relative peak factors for input signals  

- **System Identification Results**  
  - Estimated model coefficients `[A, B, C, D]` from least squares estimation  
  - Confidence intervals for identified parameters  

- **Residual and Normalization Analysis**  
  - Residual distribution plots for both input designs  
  - Normalized regressor importance and ranking  

## Contributions  
This repository is intended for academic research and educational use. Contributions and modifications are welcome.  

## License  
This project is open for research and educational purposes.  

---  
**Author:** Alexander Dowell  

