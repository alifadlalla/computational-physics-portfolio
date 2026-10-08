import json
import numpy as np
from modules import Body

def load_input(file_path):
    """Reads the input JSON file and returns structured data."""
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Convert bodies data into Body objects
    bodies = [
        Body(
            name=body['name'],
            gm=body['gm'],
            position=np.array(body['position']),
            velocity=np.array(body['velocity']),
            color=body['color']
        )
        for body in data['bodies']
    ]

    # Extract simulation parameters
    simulation_params = data['simulation']

    

    return bodies, simulation_params
