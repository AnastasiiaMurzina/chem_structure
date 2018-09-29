from math import acos, sin, pi, cos, atan, asin
import matplotlib.pyplot as plt
import numpy as np
from numpy import arctan2
from itertools import product
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from colors_from_plots import cnames


def cartesian_to_spherical(vector):
    r_xy = vector[0] ** 2 + vector[1] ** 2
    theta = arctan2(vector[1], vector[0])
    phi = arctan2(vector[2], r_xy ** 0.5)
    return theta, phi


def cube(n=3):
    side = 2./(n)
    divs = [np.array([1, -1, 1]), np.array([1,1,1])]
    divs += [divs[0]+(divs[1]-divs[0])*i/(n) for i in range(1,n)]
    for ix in range(n+1):
        divs += [np.array([divs[ix][0]-side*j, divs[ix][1], divs[ix][2]]) for j in range(1,n+1)]
    for ix in range((n+1)**2):
        divs += [np.array([divs[ix][0], divs[ix][1], divs[ix][2]-2])]
    for ix in range(3*n+1):
        divs += [np.array([divs[ix][0], divs[ix][1], divs[ix][2]-side*j]) for j in range(1,n)]
    for ix in (range(2 * (n + 1) ** 2 + 2 * (n - 1), 2 * (n + 1) ** 2 + 2 * (n - 1) + (n - 1) ** 2)):
        divs += [np.array([divs[ix][0] - 2, divs[ix][1], divs[ix][2]])]
    return divs

def show_cube(n=3):
    ds = cube(n=n)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for num, i in enumerate(ds):
        ax.scatter(i[0], i[1], i[2])
        ax.text(i[0], i[1], i[2], str(num))
    plt.show()

def plot_cube(cube_definition, n=3):
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ds = cube(n=n)
    for num, i in enumerate(ds):
        ax.scatter(i[0], i[1], i[2])
        # ax.text(i[0], i[1], i[2], str(num))
    faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))
    ax.add_collection3d(faces)
    ax.set_aspect('equal')
    plt.show()

def show_with_surface(n=3):
    cube_definition = [
        (-1,-1,-1), (-1,1,-1), (1,-1,-1), (-1,-1,1)
    ]
    plot_cube(cube_definition, n=n)

def spherical_cube(n=3):
    original = cube(n=n)
    for i in range(len(original)):
        original[i] = original[i]/np.linalg.norm(original[i])
    return original


scube = spherical_cube(3)
colors_sc = {}
name = [key for key in cnames.keys()]
for i in range(len(scube)):
    colors_sc.update({i: name[i]})
d_min = np.linalg.norm(scube[0]-scube[2])


def find_section_old(p1, p0=np.zeros(3), n=3, eps=0.05, get_error=False):
    '''
        :param p0: this point has already basis
        :param p1: not important basis of p1
        :return: section of p0 atom in which there's p1
        '''
    vec = p1 - p0
    vec = vec / np.linalg.norm(vec)
    ds = 1
    i = -1
    scube_l = spherical_cube(n=n)
    while ds > d_min/2+eps:
        if i >= len(scube_l)-1:
            i=-1
            # eps*=1.1
        i += 1
        ds = np.linalg.norm(vec - scube_l[i]) # here warning for use parameter n
    if get_error: return i, ds
    return i

def find_section(p1, p0=np.zeros(3), n=3, eps=0.05, get_error=False):
    '''
        :param p0: this point has already basis
        :param p1: not important basis of p1
        :return: section of p0 atom in which there's p1
        '''
    vec = p1 - p0
    vec = vec / np.linalg.norm(vec)
    ds = []
    scube_l = spherical_cube(n=n)
    for i, j in enumerate(scube_l):
        r = np.linalg.norm(vec - j)
        if r < d_min*2:
           ds.append([r, i])
    d, ix = sorted(ds)[0]
    if get_error: return ix, d
    return ix



def anti_scube(n=3):
    scube = spherical_cube(n)
    antis = {}
    for i, j in enumerate(scube):
        antis.update({i: find_section(-j)})
    return antis


def show_spherical_cube(n=3):
    ds = spherical_cube(n=n)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for num, i in enumerate(ds):
        ax.scatter(i[0], i[1], i[2])
        ax.text(i[0], i[1], i[2], str(num))
    plt.show()


def get_spherical_cube_in_2d(n=3):
    planar = []
    for i in spherical_cube(n=n):
        planar.append(cartesian_to_spherical(i))
    return planar


def show_planar_cube(n=3):
    for i in get_spherical_cube_in_2d(n=n):
        plt.scatter(i[0], i[1])
    plt.xlim(-pi, pi)
    plt.ylim(-pi/2, pi/2)
    plt.show()


def sc_err(n=6):
    fig = plt.figure()
    ax1 = plt.subplot(111)
    color_map = plt.get_cmap('cool')
    average_error, mx_error = 0,0
    nit = 8000
    for _ in range(nit):
        x, y, z = np.random.normal(size=3)
        r = (x**2+y**2+z**2)**0.5
        x /= r
        y /= r
        z /= r
        point = np.array([x, y, z])
        theta, phi = cartesian_to_spherical(np.array([x, y, z]))
        cur_err = np.linalg.norm(scube[find_section(point, n=n)]-point) # may be create cubes with parameter n
        d_err = mx_error - average_error
        sc = ax1.scatter(theta, phi, c=cur_err, cmap=color_map, vmin=average_error-d_err/2, vmax=average_error+d_err/2)
        average_error += cur_err / nit
        mx_error = max(mx_error, cur_err)
    ax1.set_title('quadro, n=' + str(n))
    cb = fig.colorbar(sc, orientation='horizontal', drawedges=False)
    cb.set_label(
        'Discrete errors, average=' + str(format(average_error, '.5g')) + ', Maximal = ' + str(format(mx_error, '.5g')),
        fontsize=14)
    plt.show()

def diagram_divide(n=3):
    for _ in range(10000):
        x, y, z = np.random.rand(3)
        r = np.linalg.norm(np.array([x, y, z]))
        x /= r
        y /= r
        z /= r
        x = -x if np.random.rand() < 0.5 else x # one fourth without these lines
        y = -y if np.random.rand() < 0.5 else y
        z = -z if np.random.rand() < 0.5 else z
        vec = np.array([x, y, z])
        ix = find_section(vec, n=n)
        angles = cartesian_to_spherical(vec)
        plt.scatter(angles[0], angles[1], c=colors_sc[ix])
    plt.show()


def spherical_diagram_divide(n=3):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for _ in range(10000):
        x, y, z = np.random.rand(3)
        r = np.linalg.norm(np.array([x, y, z]))
        x /= r
        y /= r
        z /= r
        x = -x if np.random.rand() < 0.5 else x # one fourth without these lines
        y = -y if np.random.rand() < 0.5 else y
        z = -z if np.random.rand() < 0.5 else z
        vec = np.array([x, y, z])
        ix = find_section(vec, n=n)
        ax.scatter(x,y,z, c=colors_sc[ix])
    plt.show()


#################Checkers#############

# def get_reversed_section_and_basis(s, b0, method='first', **kwargs):
#     if method == 'first':
#         return find_basis(np.zeros(3), rotate_by_basis(icos[s], b0[0], b0[1], **kwargs), **kwargs)
#     elif method == 'ten':
#         return find_basis(np.zeros(3), rotate_ten_vars(icos[s], b0, **kwargs), **kwargs)
#     return find_basis(np.zeros(3), rotate_non_perpendicular(icos[s], b0[0], b0[1], **kwargs), **kwargs)
#

if __name__ == '__main__':
    # anti_scube()
    # diagram_divide()
    # spherical_diagram_divide()
    pass