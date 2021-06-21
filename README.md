# Probabilistic-Robot-Localization
## Running:
`study=;`

where study can equal 1 (loop-closure study) or 2 (Monte Carlo simualtions) depending on the study of interest.

If study=2 (Monte Carlo Simulation), choose number of runs

`M_run=;`

## Evaluation:
The two localization algorithms differ only in what the robot knows about the environment.

The beacon-based localization (BBL) algorithm has prior knowledge of the landmark locations.

The SLAM algorithm does not have prior knowledge of landmark locations and allows the robot to map out unknown environments. 
