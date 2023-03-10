import numpy as np
import matplotlib.pyplot as plt

# Create sample image data
image = np.random.rand(10, 10)

# Set some values to NaN
image[image > 0.5] = np.nan

# Get the indices of the NaN values
nan_indices = np.argwhere(np.isnan(image))

# Create a figure and axis object
fig, ax = plt.subplots()

# Plot the image
ax.imshow(image, cmap='gray')

# Add a red cross at the location of each NaN value
for i in range(nan_indices.shape[0]):
    x, y = nan_indices[i]
    ax.scatter(y, x, marker='x', color='red', s=50)

# Display the plot
plt.show()