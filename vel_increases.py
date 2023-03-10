import matplotlib.pyplot as plt
import numpy as np

# Create some sample data
data1 = np.random.rand(10, 10)
data2 = np.random.rand(10, 10)
data3 = np.random.rand(10, 10)
data4 = np.random.rand(10, 10)
data5 = np.random.rand(10, 10)
data6 = np.random.rand(10, 10)

# Create a figure with three rows and two columns of subplots
fig, axes = plt.subplots(3, 2)

# Plot the data on each subplot
im1 = axes[0, 0].imshow(data1, cmap='cool')
im2 = axes[0, 1].imshow(data2, cmap='hot')
im3 = axes[1, 0].imshow(data3, cmap='cool')
im4 = axes[1, 1].imshow(data4, cmap='hot')
im5 = axes[2, 0].imshow(data5, cmap='cool')
im6 = axes[2, 1].imshow(data6, cmap='hot')

# Adjust the figure size and subplot spacing
fig.set_size_inches(8, 8)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.5, hspace=0.03)

# Display the plot
plt.show()

