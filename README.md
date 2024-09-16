<a href="https://doi.org/10.5281/zenodo.11473051"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.11473051.svg"></a>

# InflammationODE
MATLAB implementation of model code from "Mathematical Model of the Inflammatory Response to Acute and Prolonged Lipopolysaccharide Exposure in Humans"


## Repository Structure

Below is an overview of each file and its function:

### Files

1. **Main.m**
    - This is the main script to run the model. It includes code for parameter estimation and basic visualization of the results. This script integrates the various components of the model and handles the execution of the simulation.

2. **model_code.m**
    - This file contains the implementation of the Ordinary Differential Equation (ODE) model. It includes the model equations and parameters necessary for the simulation. All mathematical formulations and dynamic behaviors of the system are defined here.

3. **input_func.m**
    - This script specifies the different Lipopolysaccharide (LPS) input profiles. It defines how the input stimuli are applied to the model over time, allowing for the simulation of different experimental conditions or scenarios.

4. **data_model.m**
    - This file assigns the various datasets used for parameter estimation and plotting purposes. It links the model to the experimental data, facilitating the comparison between model predictions and observed data.

5. **cost_func.m**
    - This script defines the cost function used for parameter estimation. The cost function quantifies the difference between the model predictions and the experimental data, guiding the optimization process to find the best-fitting parameters.

6. **data.mat**
    - This file contains the population averages extracted from literature figures using WebPlotDigitizer. The data in this file serves as the reference for parameter estimation and model validation.

### Folders

In addition to the main files, there are three folders that contain scripts for different types of parameter sensitivity and profile analysis:

1. **LPSA/** (Local Parameter Sensitivity Analysis)
    - This folder contains the code for local parameter sensitivity analysis, which explores how small changes in individual parameters affect model outputs. To run this analysis, simply execute the script starting with `Main_` in this folder (e.g., `Main_LPSA.m`).

2. **MPSA/** (Multiple Parameter Sensitivity Analysis)
    - This folder includes code for multiple parameter sensitivity analysis, where multiple parameters are varied simultaneously to assess their collective influence on the model's behavior. Run the script starting with `Main_` to perform this analysis (e.g., `Main_MPSA.m`).

3. **PLA/** (Profile Likelihood Analysis)
    - This folder contains the code for profile likelihood analysis, which is used to evaluate parameter identifiability and confidence intervals. To execute this analysis, run the `Main_` script in this folder (e.g., `Main_PLA.m`).

Each folder is self-contained, and you only need to run the respective `Main_` script in each folder to perform the associated analysis.


## Usage

To run the model and perform basic visualizations, execute the `Main.m` script. This script integrates all other components and provides a comprehensive simulation environment.

### Running the Model

1. Open MATLAB.
2. Navigate to the directory containing the repository files.
3. Run the `Main.m` script by typing `Main` in the MATLAB command window.

### Modifying Input Profiles

To change the LPS input profiles, edit the `input_func.m` file. In `Main` different pre-defined datasets can be selected. Define your desired input profile within this script to simulate different experimental conditions.

### Parameter Estimation

The parameter estimation process is guided by the `cost_func.m` script. Ensure that your data in `data.mat` is correctly loaded in `Main.m` before running the parameter estimation.

## Dependencies

- MATLAB (version R2023b used for development)
- Optimization Toolbox

## Contributing

Please follow the standard GitHub workflow for contributing:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature/YourFeature`).
3. Commit your changes (`git commit -am 'Add some feature'`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Create a new Pull Request.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License. See the [LICENSE](LICENSE.md) file for details.
