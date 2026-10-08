import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models
from tensorflow.keras.callbacks import ModelCheckpoint, TensorBoard
from sklearn.model_selection import train_test_split
import os
from PIL import Image
import paths


# personal packages
from unet import unet
from data_utils import x_train, x_val, y_train, y_val

# making sure that the dataset is correct
print("Training set - Input shape:", x_train.shape, " Output shape:", y_train.shape)
print("Validation set - Input shape:", x_val.shape, " Output shape:", y_val.shape)


# ----------------------- Training ----------------------- #

# Define input and output shapes
input_shape = (320, 320, 1)
output_shape = (320, 320, 1)

# Create the U-Net model
model = unet(input_shape)

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mae'])

# change these paths
# Define the paths to save the model and logs
model_path = paths.model_path
log_dir = paths.log_dir

os.makedirs(model_path, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# Create the callbacks
checkpoint = ModelCheckpoint(model_path, monitor='val_loss', save_best_only=True)
tensorboard = TensorBoard(log_dir=log_dir)

callbacks = [checkpoint, tensorboard]

# Train the model with callbacks
batch_size = 8
epochs = 10

model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, validation_data=(x_val, y_val), callbacks=callbacks)


# Evaluate the model
loss, accuracy = model.evaluate(x_val, y_val)
print("Validation Loss:", loss)
print("Validation Accuracy:", accuracy)

#change these paths
# Create directories to save the images
input_images_dir = paths.input_images_dir
reconstructed_images_dir = paths.reconstructed_images_dir

os.makedirs(input_images_dir, exist_ok=True)
os.makedirs(reconstructed_images_dir, exist_ok=True)

# Save the input and reconstructed images
for i in range(len(x_val)):
    # Get the input image
    input_image = (x_val[i] * 255).astype(np.uint8)
    input_image_path = os.path.join(input_images_dir, f'input_{i+1}.png')
    Image.fromarray(input_image.squeeze()).save(input_image_path)
    
    # Get the reconstructed image
    reconstructed_image = model.predict(np.expand_dims(x_val[i], axis=0))
    reconstructed_image = (reconstructed_image[0] * 255).astype(np.uint8)
    reconstructed_image_path = os.path.join(reconstructed_images_dir, f'reconstructed_{i+1}.png')
    Image.fromarray(reconstructed_image.squeeze()).save(reconstructed_image_path)

