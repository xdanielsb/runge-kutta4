import matplotlib.pyplot as plt

class RungeKutta:
    def __init__(self, iters=1000, h = 0.1):
        self.iters = iters
        self.lx = []
        self.ly = []
        self.h = h

    def runge(self, f, x, y):
        for _ in range(self.iters):
            self.lx.append(x)
            self.ly.append(y)
            k1, m1 = f( x, y, self.h )
            k2, m2 = f( x + k1/2, y + m1/2, self.h)
            k3, m3 = f( x + k2/2, y + m2/2, self.h)
            k4, m4 = f( x + k3,  y + m3, self.h)
            x = x + ( k1 + 2*k2 + 2*k3 + k4) / 6.0
            y = y + ( m1 + 2*m2 + 2*m3 + m4) / 6.0

        plt.plot(self.lx, self.ly)
        plt.show()

def vanderpol(x, y, h):
    # this is the system of vanderpol
    u = 0.1 
    yp = -u*( x*x - 1) * y - x
    return y*h, yp*h

if __name__ == "__main__":
    obj = RungeKutta()
    obj.runge(vanderpol, x=2, y=2)
    
