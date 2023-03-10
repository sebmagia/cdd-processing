import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Load the image
img = Image.open('Paper_Figures/imagen.png')
width, height = img.size
scale_x = 80.42/width
scale_y = 105.2/height
pixel_x=np.arange(width+1)
pixel_y=np.arange(height+1)

x = (pixel_x -width/2)*scale_x
y = (pixel_y -height/2)*scale_y

# Convert the image to a numpy array
img_array = np.array(img)

# Plot the image
plt.imshow(img_array)

# Overlay your data on top of the image
#x = [10, 20, 30, 40, 50]
#y = [50, 40, 30, 20, 10]

#plt.plot(x, y, 'ro')

# Show the plot
plt.show()