import numpy as np
from sklearn.model_selection import train_test_split
import os
from PIL import Image
import paths


# Data preprocessing

# Define the paths to the input and output folders
input_folder = paths.input_folder
output_folder = paths.output_folder

# dict of patients and their arc
# The private patient mapping is omitted from this public archive.
patient = {}
if not patient:
    raise RuntimeError("The private dataset and patient mapping are not distributed.")

# Initialize empty lists to store input and output data
input_data = []
output_data = []

# Loop over the patients
for p, arc_list in patient.items(): 
    patient = p     # Numéro de la phase de traitement
    for arc in arc_list: # Numéro de l'arc
        
        # Define the path to the input and output folders for the current patient and arc
        input_path = os.path.join(input_folder, f"patient {patient}", f"arc_{arc}")
        output_path = os.path.join(output_folder, f"patient {patient}", f"arc_{arc}")
        
        # Loop over the images in the input folder
        for filename in os.listdir(input_path):
            # Load the input image
            input_image = Image.open(os.path.join(input_path, filename)).convert('L')
            # Resize the input image to (320, 320)
            input_image = input_image.resize((320, 320))
            # Convert PIL image to numpy array
            input_array = np.array(input_image)
            # Normalize the input array
            input_array = input_array / 255.0
            # Append the input array to the input data list
            input_data.append(input_array)

        # Loop over the images in the output folder
        for filename in os.listdir(output_path):
            # Load the output image
            output_image = Image.open(os.path.join(output_path, filename)).convert('L')
            # Resize the output image to (320, 320)
            output_image = output_image.resize((320, 320))
            # Convert PIL image to numpy array
            output_array = np.array(output_image)
            # Normalize the output array
            output_array = output_array / 255.0
            # Append the output array to the output data list
            output_data.append(output_array)

print("Done with loading the data") # Just to check if everything is working

# Convert the input and output data lists to numpy arrays
input_data = np.array(input_data)
output_data = np.array(output_data)

# Reshape the input and output data to include the channel dimension
input_data = np.expand_dims(input_data, axis=-1)
output_data = np.expand_dims(output_data, axis=-1)

# Split the dataset into training and validation sets
# 80% training, 20% validtion
x_train, x_val, y_train, y_val = train_test_split(input_data, output_data, test_size=0.2, random_state=42)

# Normalize the data
normalize_factor = 255.0
x_train = x_train / normalize_factor
x_val = x_val / normalize_factor

print("Training set - Input shape:", x_train.shape, " Output shape:", y_train.shape)
print("Validation set - Input shape:", x_val.shape, " Output shape:", y_val.shape)
