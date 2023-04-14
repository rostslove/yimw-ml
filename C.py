import hashlib
import numpy as np

class Hyperloglog:
    def __init__(self, b=10):
        self.m = 1 << b
        self.M = np.zeros(self.m, dtype=np.uint8)
        self.alpha = self.get_alpha(b)
        
    def add(self, item):
        x = int(hashlib.md5(item.encode('utf-8')).hexdigest(), 16)
        j = x & (self.m - 1) # вычисляем младшие b бит хеша x
        w = x >> self.m.bit_length() # вычисляем оставшиеся (64 - b) бит хеша x
        self.M[j] = max(self.M[j], self.rho(w) + 1)
        
    def count(self, iterations=5):
        results = []
        for i in range(iterations):
            E = self.alpha * self.m * self.m / np.sum(np.power(2, -self.M))
            if E <= self.m / 2:
                V = np.count_nonzero(self.M == 0)
                if V != 0:
                    results.append(int(-self.m * np.log(V / self.m)))
                else:
                    results.append(E)
            elif E <= 1/30 * np.power(2, 64):
                results.append(E)
            else:
                results.append(-np.power(2, 64) * np.log(1 - E / np.power(2, 64)))
            self.M.fill(0)
        return int(np.exp(np.sum(np.log(results)) / iterations))
        
    def rho(self, x):
        i = 1
        while (x & 1) == 0 and i <= 64:
            x = x >> 1
            i += 1
        return i
        
    def get_alpha(self, b):
        if b == 4:
            return 0.673
        elif b == 5:
            return 0.697
        elif b == 6:
            return 0.709
        else:
            return 0.7213 / (1 + 1.079 / (1 << b))
n = int(input())
hyper = Hyperloglog(n)
for i in range(n):
    hyper.add(input())
print(hyper.count())
