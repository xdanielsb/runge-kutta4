"""
    author:  xdanielsb
    created: 20.10.2020
"""

import matplotlib.pyplot as plt
import math


class RungeKutta:
    def __init__(self, iters=1200, h=0.1):
        self.iters = iters
        self.lx = []
        self.ly = []
        self.lt = []
        self.avg = []
        self.h = h

    def runge(self, f1, f2, x, y):
        for it in range(self.iters):
            self.lt.append(it * self.h)
            self.lx.append(x)
            self.ly.append(y)
            self.avg.append(f2(it  , self.h))
            k1, m1 = f1(x, y, self.h)
            k2, m2 = f1(x + k1 / 2, y + m1 / 2, self.h)
            k3, m3 = f1(x + k2 / 2, y + m2 / 2, self.h)
            k4, m4 = f1(x + k3, y + m3, self.h)
            x = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            y = y + (m1 + 2 * m2 + 2 * m3 + m4) / 6.0
        fr = 0
        tr = 800
        plt.figure(1)
        plt.plot(self.lx[fr:], self.ly[fr:], color="g", label="f_aprox")
        plt.plot([ x*1 for x in self.lx[fr:tr]],[y*1 for y in self.avg[fr:tr]], color="b", label="Avg")
        plt.gca().set_aspect("equal", adjustable="box")
        plt.legend()
        plt.figure(2)
        plt.plot(self.lt[fr:], self.ly[fr:], color="g", label="f_aprox")
        plt.plot(self.lt[fr:], self.avg[fr:], color="b", label="Avg")
        leg = plt.legend(loc="best", ncol=2, mode="expand", shadow=True, fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.show()


class VanderPol:
    # exact solution
    def exact(self, t, h):
        # this is the solution given by the average equations
        # of the vanderpol system
        ep = 0.1
        t = (t-1)*h
        ans = 2.0 / math.sqrt(3 * (math.e **t) + 1)
        ans = ans * math.cos(t)
        return ans

    # dynamic system
    def system(self, x, y, h):
        # this is the system of vanderpol
        # expresed in diferential equations
        u = 0.1
        yp = -u * (x * x - 1) * y - x
        return y * h, yp * h


class SystemDT:
    # exact solution
    def exact(self, t, h):
        ep = 0.1
        # changing this t, we rotate the system ?
        t = (t + 6)*h
        ans  = math.sqrt(6/ ( 1 + math.e ** (-(t/2)*ep)))
        ans = ans * math.cos(t)
        return ans
    # dynamic system
    def system(self, x, y, h):
        # this is the system of vanderpol
        # expresed in diferential equations
        u = 0.1
        yp = -u * (x * x - 1) * (y ** 3) - x
        return y * h, yp * h


if __name__ == "__main__":
    # vd = VanderPol()
    vd = SystemDT()
    RungeKutta().runge(vd.system, vd.exact, x=0.5, y=0)
