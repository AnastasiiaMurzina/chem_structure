from math import tan, sin, pi, sqrt

'''
u from -pi to pi
v from -pi/2 to pi/2
u              |  sqrt(2/3)sin u < tg v |  sqrt(2/3)sin u > tg v | sqrt(2/3)sin u > tg v | sqrt(2/3)sin u < tg v |
               | -sqrt(2/3)sin u < tg v | -sqrt(2/3)sin u < tg v |-sqrt(2/3)sin u > tg v |-sqrt(2/3)sin u > tg v |
                
abs(u)>pi/2    |            1           |          2  |        3 |               4       |      5     |      6   |
abs(u)<pi/2    |            7           |          8  |        9 |               10      |     11     |     12   |

               |                        | v>0         | v<0      |                       |v>0         |      v<0 |
it's smth hard!!!!     
'''

def divider(u,v):
    left_sin = (2/3)*sin(u)
    right_tan = tan(v)
    if left_sin < right_tan:
        if -left_sin < right_tan:
            if abs(u)>pi/2:
                return 1
            else:
                return 7
        else:
            if abs(u)>pi/2:
                if v>0:
                    return 5
                else:
                    return 6
            else:
                if v>0:
                    return 11
                else:
                    return 12
    else:
        if -left_sin < right_tan:
            if abs(u) > pi / 2:
                if v > 0:
                    return 2
                else:
                    return 3
            else:
                if v > 0:
                    return 8
                else:
                    return 9
        else:
            if abs(u)>pi/2:
                return 4
            else:
                return 10

def divider2(u,v):
    left_sin = (2/3)*sin(u)
    right_tan = tan(v)
    if left_sin < right_tan:
        if -left_sin < right_tan:
            if abs(u)>pi/2:
                return 1
            else:
                return 7
        else:
            if abs(u)>pi/2:
                if v>0:
                    return 5
                else:
                    return 6
            else:
                if v>0:
                    return 11
                else:
                    return 12
    else:
        if -left_sin < right_tan:
            if abs(u) > pi / 2:
                if v > 0:
                    return 2
                else:
                    return 3
            else:
                if v > 0:
                    return 8
                else:
                    return 9
        else:
            if abs(u)>pi/2:
                return 4
            else:
                return 10

def divider12(u, v):
    local = 0
    if -pi<u<-2*pi/3:
        local = 1
    elif -2*pi/3< u < -pi/3:
        local = 3
    elif -pi/3<u<0:
        local = 5
    elif 0<u<pi/3:
        local = 7
    elif pi/3<u<2*pi/3:
        local = 9
    elif 2*pi/3< u < pi:
        local = 11
    if v<0:
        local+=1
    return local

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    from cmath import pi
    from math import sin, cos, tan
    from mpl_toolkits.mplot3d import Axes3D


    def anglesToCoords(u, v):
        return np.array([np.cos(u) * np.cos(v), np.sin(u) * np.cos(v), np.sin(v)])

    def graphic_test():
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        u, v = np.mgrid[-pi:pi:20j, -pi / 2:pi / 2:10j]
        x, y, z = anglesToCoords(u, v)
        ax.plot_wireframe(x, y, z, color="r", alpha=0.5)
        import random
        color_dict = {1: 'r', 2: 'b', 3: 'g', 4: 'y', 5: 'pink', 6: 'sienna', 7: 'darkred', 8: 'lime',
                      9: 'orange', 10: 'crimson', 11: 'purple', 12: 'olive'}
        for i in range(1000):
            alpha = 2 * pi * random.random() - pi
            betta = pi * random.random() - pi / 2
            point = anglesToCoords(alpha, betta)
            ax.scatter(point[0], point[1], point[2], color=color_dict[divider(alpha, betta)])

            # if divider(alpha, betta) == 1:
            #     point = anglesToCoords(alpha, betta)
            #     ax.scatter(point[0], point[1], point[2], color='r')
        plt.show()

    def square_test(n):
        import  random
        arr = [0 for i in range(12)]
        for i in range(n):
            alpha = 2 * pi * random.random() - pi
            betta = pi * random.random() - pi / 2
            arr[divider12(alpha,betta)-1]+=1
        print([i/n for i in arr])

    square_test(10000)
