#!/usr/bin/python3


import os

import matplotlib.pyplot    as plt
import numpy                as np

file = os.getenv('file_image')
img = os.getenv('final_img')
data = np.fromfile(file, dtype=np.float32)
nx = int(os.getenv('nx'))
ny = int(os.getenv('ny'))
image = data.reshape((nx, ny))

plt.imshow(image, cmap='gray')
plt.axis('off')
plt.savefig(img, bbox_inches='tight', pad_inches=0)
