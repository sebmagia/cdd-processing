import itertools
from multiprocessing import Pool

# attributes is always accessible as a global, so worker processes can directly access it
attributes = ('Age', 'Workclass', 'Fnlwgt', 'Education', 'Education-num',
              'marital-status', 'Occupation', 'Relationship', 'Race', 'Sex',
              'Capital-gain', 'Capital-loss', 'Hours-per-week', 'Native country',
              'Probability', 'Id')

def comb(n): # the argument n is the number of items to select
    res = list(itertools.combinations(attributes, n)) # create a list from the iterator
    return res

def main():
    p = Pool(4)
    times = range(0, len(attributes)+1)
    values = p.map(comb, times) # pass the range as the sequence of arguments!
    p.close()
    p.join()
    print(values)

if __name__ == '__main__':
    main()
