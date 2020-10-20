"""
    author:  xdanielsb
    created: 20.10.2020
"""

import matplotlib.pyplot as plt
import math

class RungeKutta:
    def __init__(self, iters=1200, h = 0.1):
        self.iters = iters
        self.lx = []
        self.ly = []
        self.lt = []
        self.avg = []
        self.h = h

    def runge(self, f1, f2, x, y):
        for it in range(self.iters):
            self.lt.append(it*self.h)
            self.lx.append(x)
            self.ly.append(y)
            self.avg.append( f2(it+5, self.h))
            k1, m1 = f1( x, y, self.h )
            k2, m2 = f1( x + k1/2, y + m1/2, self.h)
            k3, m3 = f1( x + k2/2, y + m2/2, self.h)
            k4, m4 = f1( x + k3,  y + m3, self.h)
            x = x + ( k1 + 2*k2 + 2*k3 + k4) / 6.0
            y = y + ( m1 + 2*m2 + 2*m3 + m4) / 6.0
        fr = 200
        plt.figure(1)
        plt.plot(self.lx[fr:], self.ly[fr:], 'g:', self.lx[fr:], self.avg[fr:], 'b-')
        plt.figure(2)
        plt.plot(self.lt[fr:], self.ly[fr:], 'g:', self.lt[fr:], self.avg[fr:], 'b-')
        plt.show()


def exact_solution(t, h):
    # this is the solution given by the average equations
    # of the vanderpol system 
    ep = 0.1
    t = (t-1)*h
    ans = 2.0 / math.sqrt ( 3* (math.e ** (-t*ep)) + 1)
    ans = ans * math.cos(t)
    return ans


def vanderpol(x, y, h):
    # this is the system of vanderpol
    # expresed in diferential equations
    u = 0.1
    yp = -u*( x*x - 1) * y - x
    return y*h, yp*h


if __name__ == "__main__":
    obj = RungeKutta()
    obj.runge(vanderpol, exact_solution, x=1, y=2)
    
