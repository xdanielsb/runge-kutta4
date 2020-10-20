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
        self.avg = []
        self.h = h

    def exact_solution(self, t):
        # this is the solution given by the average equations
        # of the vanderpol system 
        ep = 0.1
        t = (t-1)*self.h
        ans = 2.0 / math.sqrt ( 3* (math.e ** (-t*ep)) + 1)
        ans = ans * math.cos(t)
        return ans


    def runge(self, f, x, y):
        for it in range(self.iters):
            self.lx.append(x)
            self.ly.append(y)
            self.avg.append( self.exact_solution(it+5))
            k1, m1 = f( x, y, self.h )
            k2, m2 = f( x + k1/2, y + m1/2, self.h)
            k3, m3 = f( x + k2/2, y + m2/2, self.h)
            k4, m4 = f( x + k3,  y + m3, self.h)
            x = x + ( k1 + 2*k2 + 2*k3 + k4) / 6.0
            y = y + ( m1 + 2*m2 + 2*m3 + m4) / 6.0
        fr = 0
        plt.plot(self.lx[fr:], self.ly[fr:], 'g:', self.lx[fr:], self.avg[fr:], 'b-')
        plt.show()



def vanderpol(x, y, h):
    # this is the system of vanderpol
    # expresed in diferential equations
    u = 0.1
    yp = -u*( x*x - 1) * y - x
    return y*h, yp*h


if __name__ == "__main__":
    obj = RungeKutta()
    obj.runge(vanderpol, x=1, y=2)
    
