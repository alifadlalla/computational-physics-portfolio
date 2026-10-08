import numpy as np

# Function to define the OR logic gate
def or_gate(x1, x2):
    inputs = np.array([x1, x2, 1])  # Adding bias node with value 1
    weights = np.random.rand(3)      # Initialize weights with small random values
    output = np.dot(inputs, weights) > 0
    return int(output)

# Training data for OR gate
input_data = np.array([
    [0, 0],
    [0, 1],
    [1, 0],
    [1, 1]
])

output_data = np.array([0, 1, 1, 1])  # OR gate outputs

# Training the Perceptron
learning_rate = 0.1
epochs = 10

print("Initial weights:", weights)

for epoch in range(epochs):
    for i in range(len(input_data)):
        x1, x2 = input_data[i]
        target_output = output_data[i]
        
        inputs = np.array([x1, x2, 1])
        predicted_output = np.dot(inputs, weights) > 0
        error = target_output - predicted_output
        
        weights += learning_rate * error * inputs

    print(f"Iteration {epoch + 1}: Weights = {weights}")

# Testing the trained Perceptron
print("\nTesting the trained Perceptron:")
for i in range(len(input_data)):
    x1, x2 = input_data[i]
    predicted_output = or_gate(x1, x2)
    print(f"Input: ({x1}, {x2}) | Predicted Output: {predicted_output}")
