from itertools import product
import numpy as np

def my_func(a,b,c,d): #Function that I have to evaluate
    #do some matrix computations and return some scalar q#
    return q

Range = [x for x in np.arange(0,3.60,0.5)] # A list that contains 8 elements
it = product(Range, repeat=4) #Generate all possiblities
Max = 0
for i in it:
        Sum = my_func(*i)
        if Sum > Max:
            Max = Sum
#return Max

#Result = Evaluate() #Output


# with Pool(processes=4) as pool:
#   pool.map(my_func, itertools.product(Range, repeat=8))
