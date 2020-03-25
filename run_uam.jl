# Basic test script - generates the file sample_uam_traj.csv containing one uam trajectory

# Include the files
include("UAMTrajectoryGenerator.jl")
# Set the random seed
Random.seed!(1)
# Generate the trajectory file
generate_trajectory_file("sample_uam_traj.csv")