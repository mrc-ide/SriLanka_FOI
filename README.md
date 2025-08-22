# ShinyApp: Dengue FOI and Model Analysis

The **Dengue FOI and Model Analysis** Shiny application is designed to analyze the Force of Infection (FOI) and epidemiological models for dengue cases across selected provinces and districts in Sri Lanka.  

It allows users to:  
- Upload datasets  
- Filter by province and year  
- Run statistical models (either time-constant or time-varying analysis)
- Generate interactive plots and metrics for analysis  

---

## Getting Started

1. Open **RStudio** and the file `ShinyApp.R`. Ensure the following packages have been installed and otherwise run: install.packages("shiny", "ggplot2", "ggpubr", "rstan", "dplyr", "tidyr", "reshape2", "writexl"))
2. Set the working directory to the **ShinyApp folder**.  
3. Highlight the entire code and run it.  
4. The ShinyApp window will appear automatically.  

---

## 1. Data Upload and Configuration

- **File Inputs**: Requires two CSV files:  
  - Dengue cases (simulated)  
  - Population data  
  Columns should include *province, district, year, and age group*.  

- **Separator Selection**: Choose between comma, semicolon, or tab delimiters.  

- **Province and Year Selection**:  
  - Dropdown for province selection  
  - Slider for year range  

---

## 2. Model Selection and Parameters

- **Models Available**:  
  - `time_constant.stan`: FOI assumed constant across time  
  - `time_varying.stan`: Year-specific FOI estimates  

- **Iterations**:  
  - Slider to select number of iterations  
  - Increase iterations if the model does not converge (e.g., Rhat > 1.1 or poor chain convergence).  

- **Run Model**:  
  - Click **Run Model** to start.  
  - Progress bar shows model fitting stages.  
  - Results appear once complete.  
  - If an error (`incorrect number of dimensions`) occurs, it will resolve automatically after the process finishes.  

---

## 3. Plots and Outputs

- **Fit Plot**:  
  - Observed vs predicted incidence by year, age group, and location  
  - Median (blue line) with 95% credible intervals (blue shaded area)  

- **FOI (Lambda) Plot**:  
  - Historical FOI (Lam_H)  
  - Time-varying FOI (for `time_varying.stan`)  
  - 95% credible intervals  

- **Parameter (Pars) Plot**:  
  - Province-specific parameters (rho, gamma, chi1, chi2)  
  - Median with 95% credible intervals  

- **Traceplot**:  
  - MCMC chain diagnostics and convergence  

- **Rhat Values Table**:  
  - Summary of Rhat values to assess convergence  

---

## 4. Model Output Management

- **Model Feedback Tab**:  
  - Real-time updates on model progress and potential issues  

- **Excel Output Files**:  
  - For `time_constant.stan`: **one Excel file** with posterior parameter estimates  
  - For `time_varying.stan`: **two Excel files** with posterior parameter estimates  
  - Files are automatically saved in the working folder and include summaries of the posterior distributions  

---
