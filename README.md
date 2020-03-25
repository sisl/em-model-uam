# em-model-uam

Repository for generating UAM trajectories and write them to output CSV files.

## Quick Start Guide
Ensure that the following packages are installed in your current version of Julia: `Distributions`, `Convex`, `ECOS`, `DataFrames`, `CSV`.

The following lines will generate a trajectory file titled "test.csv" with data for a nominal landing trajectory.

```
include("UAMTrajectoryGenerator.jl")
τ = nominal_landing()
generate_trajectory_file(τ, "test.csv")
```

If no trajectory type is specified, a random trajectory type will be chosen and generated (uniform among all possible types):

```
generate_trajectory_file("test.csv")
```

If you want to generate more than 1 trajectory per file, you can add this as the last required argument. The following line will generate a file called "test.csv" with 3 trajectories.

```
generate_trajectory_file("test.csv", 3)
```

The first column of the resulting CSV file will be the time in seconds and the following three columns will be the x-, y-, and z-position respectively in feet. To test on your system, include the `run_uam.jl` script. It should generate a file that matches the "sample_uam_traj.csv" file.

## Type Descriptions
`SAMPLER` - abstract type of which specific trajectory samplers are subtypes. A sampler contains the distributions of random variables associated with a particular trajectory type (e.g. initial altitude, glideslope, etc.)

`UAM_TRAJECTORY` - abstract type of which specific trajectory types such as `NOMINAL_LANDING` and `VERTICAL_ASCENT` are subtypes. All trajectory types contain a corresponding sampler, a time step, acceleration constraints, position, velocity, and a place to store each variable in the sampler. 

`TRAJECTORY` - contains position (`p`), velocity (`v`), and acceleration (`a`) for an entire trajectory. Each row is the x, y, and z values of the variable at a particular time step with time increasing with the rows.

## Key Functions
`sample_features!(τ::UAM_TRAJECTORY)` - replaces feature values in the trajectory by sampling new values according to the distributions specified in the sampler.

`generate_trajectory!(τ::UAM_TRAJECTORY)` - runs `sample_features!(τ)` on the trajectory and uses the sampled features to generate the trajectory by solving a constrained convex optimization problem. When the `solve_trajectory!(τ::UAM_TRAJECTORY)` function is called, the position `p` and velocity `v` in the trajectory object will be replaced with the values for the newly generated trajectory.

If you want to hold any features constant during the trajectory, you will need to add a line after `sample_features!(τ)` gets called. For example, if you know you want a nominal landing with an starting altitude of 500 ft, using the following lines of code to generate the trajectory:

```
τ = nominal_landing()
sample_features!(τ)
τ.z_init = 500
solve_trajectory!(τ)
```

## File Descriptions
`UAMTrajectoryGenerator.jl` - Main file to include for UAM trajectory generation that sets up constants and types; includes all necessary files, and defines general functions for all types of trajectories.

`HighReconaissance.jl` - defines functions and types for high reconaissance trajectories.

`Nominal Landing.jl` - defines optimization problem and types for nominal landing trajectories.

`NominalTakeoff.jl` - defines optimization problem and types for nominal takeoff trajectories.

`VerticalAscent.jl` - defines optimization problem and types for vertical ascent takeoff trajectories.

`VerticalDescent.jl` - defines optimization problem and types for vertical descent landing trajectories.

`run_uam.jl` - sample script for generating uam trajectory files

`sample_uam_traj.csv` - sample model output for one uam trajectory
