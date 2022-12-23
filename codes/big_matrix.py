import numpy as np

nrows, ncols = 1000000, 100

f = np.memmap('memmapped.dat', dtype=np.float32,
              mode='w+', shape=(nrows, ncols))
