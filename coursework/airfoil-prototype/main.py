import yaml
from airfoil_solver import AirfoilSolver
import time

def read_config(file_path="config.yaml"):
    with open(file_path, "r") as file:
        config = yaml.safe_load(file)
    return config

if __name__ == "__main__":

     # Record the start time
    start_time = time.time()


    # Read input configuration from YAML file
    config = read_config()

    # Create an instance of AirfoilSolver using the configuration
    airfoil_solver = AirfoilSolver(**config)

    # visualize   
    airfoil_solver.plot_airfoil()
    airfoil_solver.visualize(save=True)
 
    # Record the end time
    end_time = time.time()

    print(f"Total Runnung time: {end_time - start_time} seconds")
