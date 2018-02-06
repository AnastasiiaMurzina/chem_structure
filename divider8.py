'''
u from -pi to pi
v from -pi/2 to pi/2
u   |  0..pi/2 |  pi/2..pi|-pi..-pi/2|  -pi/2..0|
    |  1  |  2 |  3  |  4 |  5  |  6 |  7  |  8 |
v   |>0   | <0 |>0   | <0 |>0   | <0 |>0   | <0 |
'''

from math import pi
# return num of zone with divided 8
def divider8(u, v):
    if 0 < u < pi/2:
        if v < 0:
            return 1
        else:
            return 2
    elif pi/2 < u < pi:
        if v < 0:
            return 3
        else:
            return 4
    elif -pi < u < -pi/2:
        if v < 0:
            return 5
        else:
            return 6
    else:
        if v < 0:
            return 7
        else:
            return 8

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import cmath
    from cmath import pi
    from math import sin, cos, tan
    from mpl_toolkits.mplot3d import Axes3D


    def anglesToCoords(u, v):
        return np.array([np.cos(u) * np.cos(v), np.sin(u) * np.cos(v), np.sin(v)])
    def graphic_test():
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        u, v = np.mgrid[-pi:pi:20j, -pi/2:pi/2:10j]
        x, y, z = anglesToCoords(u, v)
        ax.plot_wireframe(x, y, z, color="r", alpha=0.5)
        import random
        color_dict = {1: 'r', 2: 'b', 3: 'g', 4: 'y', 5: 'pink', 6: 'sienna', 7: 'darkred', 8: 'lime'}
        for i in range(1000):
            alpha = 2*pi*random.random()-pi
            betta = pi*random.random()-pi/2
            point = anglesToCoords(alpha, betta)
            ax.scatter(point[0], point[1], point[2], color=color_dict[divider(alpha, betta)])

            # if divider(alpha, betta) == 1:
            #     point = anglesToCoords(alpha, betta)
            #     ax.scatter(point[0], point[1], point[2], color='r')
        plt.show()

    def square_test(n):
        import  random
        arr = [0 for i in range(8)]
        for i in range(n):
            alpha = 2 * pi * random.random() - pi
            betta = pi * random.random() - pi / 2
            arr[divider8(alpha,betta)-1]+=1
        print([i/n for i in arr])

    square_test(100000)
