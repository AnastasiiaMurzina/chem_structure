import numpy as np
from numpy import arctan2


def cartesian_to_spherical(vector):
    r_xy = vector[0] ** 2 + vector[1] ** 2
    theta = arctan2(vector[1], vector[0])
    phi = arctan2(vector[2], r_xy ** 0.5)
    return theta, phi

class Spherical_divider():
    def __init__(self, n):
        self.n = n
        side = 2. / n
        divs = [np.array([1, -1, 1]), np.array([1, 1, 1])]
        divs += [divs[0] + (divs[1] - divs[0]) * i / n for i in range(1, n)]
        for ix in range(n + 1):
            divs += [np.array([divs[ix][0] - side * j, divs[ix][1], divs[ix][2]]) for j in range(1, n + 1)]
        for ix in range((n + 1) ** 2):
            divs += [np.array([divs[ix][0], divs[ix][1], divs[ix][2] - 2])]
        for ix in range(3 * n + 1):
            divs += [np.array([divs[ix][0], divs[ix][1], divs[ix][2] - side * j]) for j in range(1, n)]
        for ix in (range(2 * (n + 1) ** 2 + 2 * (n - 1), 2 * (n + 1) ** 2 + 2 * (n - 1) + (n - 1) ** 2)):
            divs += [np.array([divs[ix][0] - 2, divs[ix][1], divs[ix][2]])]
        for i in range(len(divs)):
            divs[i] = divs[i] / np.linalg.norm(divs[i])
        self.scube = divs
        self.d_min = np.linalg.norm(divs[0]-divs[2])
        self.anti_scube = {}
        self.set_anti_scube()


    def find_section(self, p1, p0=np.zeros(3), get_error=False):
        """
            :param p0: this point has already basis
            :param p1: not important basis of p1
            :return: section of p0 atom in which there's p1
        """
        vec = p0 - p1
        vec = vec / np.linalg.norm(vec)
        ds = []
        for i, j in enumerate(self.scube):
            r = np.linalg.norm(vec - j)
            if r < 2 * self.d_min:
               ds.append([r, i])
        d, ix = sorted(ds)[0]
        if get_error: return ix, d
        return ix

    def nearest_from_to(self, s1, s2):
        v1 = -self.scube[s1] # Warning: with "-"
        v2 = -self.scube[s2]
        dv = v2-v1
        k = self.n
        for i in range(k):
            current = v1+i*dv/k
            current = current/np.linalg.norm(current)
            fs = self.find_section(current)
            if fs != s1:
                return fs
        return s2

    def set_anti_scube(self):
        for i, j in enumerate(self.scube):
            self.anti_scube.update({i: self.find_section(j)})

if __name__ == '__main__':
    # pass
    n = 5
    s = Spherical_divider(n)
    k1, k2 = 5, 10
    ss = s.nearest_from_to(k1, k2)
    print(np.linalg.norm(s.scube[k1]-s.scube[k2]))
    print(np.linalg.norm(s.scube[k1]-s.scube[ss]))
    # k = len(s)
    # d = [np.linalg.norm(s1-s2) for s1, s2 in (combinations(s, 2))]
    # print(min(d))
    # print(len(np.where(d == min(d))[0]))

