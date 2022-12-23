from multiprocessing import Pool

class Test:

    def __init__(self):
        self.main()

    @classmethod
    def methodForMultiprocessing(cls, x):
        print(x*x)

    def main(self):
        if __name__ == "__main__":
            p = Pool()
            p.map(Test.methodForMultiprocessing, list(range(1, 11)))
            p.close()

TestObject = Test()

from functools import partial

def multiply(x, y,z,w):
    print(x*y + z*w)


g=partial(multiply,2,3)