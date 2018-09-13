from math import acos, sin, pi, cos, atan, asin
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D

'''
Describe division: first parameter of bond is section [0..11]
second parameter is rotate by two numbers [0..n_y-1, 0..n_z-1]
P.S. if section is 0 or 6, rotate's described by [0..n_y-1, 0]
'''

r = 1
a = 4 * r * (10 + 22 / 5 ** 0.5) ** (-0.5)  # length of the side
d = a * sin(3 * pi / 10) / sin(2 * pi / 5)  # distance between the center and the pentagon angle
h = d * sin(3 * pi / 10)  # hight to the middle of the side
angle = acos(- 5 ** (-0.5))  # between two surfaces
center_up = np.array([0, 0, 1])
center_down = np.array([0, 0, -1])
center = np.array([0, 0, 0])
min_distance = h / 2  # for decoding position
# Rotation consts
eps_length = 0.001
d_hor_angle = 2 * pi / 5
d_ver_angle = pi / 10
n_y = 4
n_z = 4


##########################Build dodecahedron##############################
def get_side_centers_hotizontal(point):
    '''
    :param point: center_point of dodecahedron
    :return: five points - centers of horizontal pentagon sides, Radious = 1 (global)
    '''
    initial_angle = 0
    centers = []
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) + point[0], h * sin(initial_angle) + point[1], point[2]]))
        initial_angle += 2 * pi / 5
    return centers


def get_up_surface_centers(point_from):
    '''
    :param point_from: the highest point of dodecahedron
    :return: points in the higher semispace of dodecahedron
    '''
    initial_angle = pi/5
    centers = [point_from]
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) * (1 - cos(angle)) + point_from[0],
                                 h * sin(initial_angle) * (1 - cos(angle)) + point_from[1],
                                 h * sin(angle) + point_from[2]]))
        initial_angle += 2 * pi / 5
    return centers


def get_down_surface_centers(point_from):
    '''
    :param point_from: the lowest point of dodecahedron
    :return: points in thr lower semispace of dodecahedron
    '''
    initial_angle = 0
    centers = [point_from]
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) * (1 - cos(angle)) + point_from[0],
                                 h * sin(initial_angle) * (1 - cos(angle)) + point_from[1],
                                 -h * sin(angle) + point_from[2]]))
        initial_angle += 2 * pi / 5
    return centers


def get_penta_points(point_center=np.zeros(3)):
    '''
    :param point_center:
    :return: centers of dodecahedron (oriented with (0, 0, 0) degrees) sides
    '''
    center_up = np.array([i for i in point_center])
    center_up[2] += 1
    center_down = np.array([i for i in point_center])
    center_down[2] -= 1
    return np.concatenate((get_down_surface_centers(center_up), get_up_surface_centers(center_down)), axis=0)


pp = get_penta_points()

def ico_points_from_penta():
    icos = []
    for i in range(1,5):
        icos.append(np.array([(pp[0][j]+pp[i][j]+pp[i+1][j])/3. for j in range(3)]))
        icos.append(np.array([(pp[6+i][j]+pp[i][j]+pp[i+1][j])/3. for j in range(3)]))
    icos.append(np.array([(pp[0][j]+pp[1][j]+pp[5][j])/3. for j in range(3)]))
    icos.append(np.array([(pp[11][j]+pp[1][j]+pp[5][j])/3. for j in range(3)]))
    for i in range(7,11):
        icos.append(np.array([(pp[6][j]+pp[i][j]+pp[i+1][j])/3. for j in range(3)]))
        icos.append(np.array([(pp[i-5][j]+pp[i][j]+pp[i+1][j])/3. for j in range(3)]))
    icos.append(np.array([(pp[6][j]+pp[7][j]+pp[11][j])/3. for j in range(3)]))
    icos.append(np.array([(pp[1][j]+pp[7][j]+pp[11][j])/3. for j in range(3)]))
    return [i/np.linalg.norm(i) for i in icos]

icos = ico_points_from_penta()
###################################################

def show_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in points:
        ax.scatter(i[0], i[1], i[2])
    plt.show()


def show_named_points(named_points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for key, i in named_points.items():
        ax.scatter(i[0][0], i[0][1], i[0][2])
        ax.text(i[0][0], i[0][1], i[0][2], str(key))
    plt.show()

def show_named_points1(named_points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for key, i in named_points.items():
        ax.scatter(i[0], i[1], i[2])
        ax.text(i[0], i[1], i[2], str(key))
    plt.show()

#########################angles###########################
def show_points_by_angles(center, angles, labels=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if labels:
        for item, i in angles.items():
            ax.scatter(center[0] + cos(i[0]) * sin(i[1]), center[1] + sin(i[0]) * sin(i[1]), center[2] + cos(i[1]))
            ax.text(center[0] + cos(i[0]) * sin(i[1]), center[1] + sin(i[0]) * sin(i[1]), center[2] + cos(i[1]), item)
    else:
        for item, i in angles.items():
            ax.scatter(center[0] + cos(i[0]) * cos(i[1]), center[1] + sin(i[0]) * cos(i[1]), center[2] + sin(i[1]))
    plt.show()


def get_penta_angles():
    angles = [np.array([0, pi / 2])]
    init_angle = 0
    for i in range(5):
        init_angle += (2 * pi / 5) % (2 * pi)
        angles.append(np.array([init_angle, asin(sin(angle) / 2)]))
    # print(sin(angle)/2)
    angles.append(np.array([0, -pi / 2]))
    for i in range(5):
        init_angle += (2 * pi / 5) % (2 * pi)
        angles.append(np.array([init_angle, -asin(sin(angle) / 2)]))
    return angles
##############################Rotate dodecahedron#########################################
def Rx(x_angle):
    '''
    :param x_angle:
    :return: matrix for rotation around x-axis with x-angle
    '''
    R_x = np.eye(3)
    R_x[1][1], R_x[1][2], R_x[2][1], R_x[2][2] = [cos(x_angle), -sin(x_angle), sin(x_angle), cos(x_angle)]
    return R_x


def Ry(y_angle):
    '''
    :param y_angle:
    :return: matrix for rotation around y-axis with y-angle
    '''
    R_y = np.eye(3)
    R_y[0][0], R_y[0][2], R_y[2][0], R_y[2][2] = [cos(y_angle), sin(y_angle), -sin(y_angle), cos(y_angle)]
    return R_y


def Rz(z_angle):
    '''
    :param z_angle:
    :return: matrix for rotation around z-axis with z-angle
    '''
    R_z = np.eye(3)
    R_z[0][0], R_z[0][1], R_z[1][0], R_z[1][1] = [cos(z_angle), -sin(z_angle), sin(z_angle), cos(z_angle)]
    return R_z


def rotate_by_basis(point, y, z, **kwargs):
    '''
    :param point: array or one point to rotate
    :param y: first coordinate of rotation
    :param z: second coordinate of rotation
    :param kwargs: n_y and n_z, default 4
    :return: rotated point or points
    '''
    n_y_l = kwargs.get('n_y', n_y)
    n_z_l = kwargs.get('n_z', n_z)
    operator = Rz(d_hor_angle * z / n_z_l).dot(Ry(d_ver_angle * y / n_y_l))
    if isinstance(point, (np.ndarray)):
        return point.dot(operator)
    return [i.dot(operator) for i in point]


def rotate_ten_vars(point, i1, **kwargs):
    radius = kwargs.get('r', 0.7)
    b = list(product([-radius, 0, radius], repeat=2))
    operator = Rz(d_hor_angle*b[i1][0]).dot(Ry(d_ver_angle*b[i1][1]))
    pp_ = [i.dot(operator) for i in icos]
    for i in icos:
        pp_.append(i)
        # show_points(pp_)
    # operator = []
    if isinstance(point, (np.ndarray)):
        return point.dot(operator)
    return [i.dot(operator) for i in point]

def rotate_non_perpendicular(point, y, z, **kwargs): # fr=1..2; sr=6..10
    '''
    :param point: array or one point to rotate
    :param y: first coordinate of rotation
    :param z: second coordinate of rotation
    :param kwargs: n_y, n_z, fr, sr - parameters of exist basis
    :return:
    '''
    n_y_l = kwargs.get('n_y', n_y)
    n_z_l = kwargs.get('n_z', n_z)
    fr = kwargs.get('fr', 9)
    sr = kwargs.get('sr', 13)

    arounder = icos[fr]-icos[sr]
    W = np.zeros((3, 3))
    W[0, 1] = -arounder[2]
    W[0, 2] = arounder[1]
    W[1, 0] = arounder[2]
    W[1, 2] = - arounder[0]
    W[2, 0] = - arounder[1]
    W[2, 1] = arounder[0]
    phi = d_ver_angle*y/n_y_l
    rotated_matrix = np.eye(3)+np.sin(phi)*W+(2*(np.sin(phi/2)**2))*W**2
    operator = Rz(d_hor_angle * z / n_z_l).dot(rotated_matrix)
    if isinstance(point, (np.ndarray)):
        rot_p = point.dot(operator)
        return rot_p/np.linalg.norm(rot_p)
    rot_ps = [i.dot(operator) for i in point]
    return [i/np.linalg.norm(i) for i in rot_ps]

def rotate_radius_vars(point, i1, **kwargs):
    ''' Is it need to develope this method better?????
    :param point:
    :param i1:
    :param kwargs:
    :return:
    '''
    radius = kwargs.get('r', 0.7)
    ns = kwargs.get('ns', 3)
    b = list((np.linspace(-radius, radius, ns)))
    operator = Rz(ns*d_hor_angle*b[i1]).dot(Ry(ns*d_ver_angle*(r-b[i1])))
    pp_ = [i.dot(operator) for i in icos]
    for i in icos:
        pp_.append(i)
        # show_points(pp_)
    # operator = []
    if isinstance(point, (np.ndarray)):
        return point.dot(operator)
    return [i.dot(operator) for i in point]



def point_to_angles(point):
    '''
    :param point: coordinates in Decart basis
    :return: two angles in spherical coordinates (without knowing radius)
    '''
    if isinstance(point[0], (int, float)):
        if point[0] != 0:
            theta = atan(point[1] / point[0])
            if point[0] < 0:
                theta += pi
        else:
            if point[1] > 0:
                theta = pi / 2
            elif point[1] < 0:
                theta = -pi / 2
            else:
                theta = 0
        return np.array([theta, acos((point[2]))])
    else:
        angles = []
        for i in point:
            angles.append(point_to_angles(i))
        return angles


def check_diff_rotate(n_y, n_z):
    '''
    :param n_y: count of vertical possible rotation positions
    :param n_z: count of horizontal possible rotation positions
    :return: minimal distance between one point in two different rotations
    '''
    dhs = []
    for i in range(n_y):
        for j in range(n_z):
            dhs.append(np.array(rotate_by_basis(icos, i, j, **{'n_y': n_y, 'n_z': n_z})))
    diff = 1
    for i0, i in enumerate(dhs):
        for j0, j in enumerate(dhs):
            if i0 != j0 and i0 + j0 != 0:
                for k in range(12):
                    if np.linalg.norm(i[k] - j[k]) != 0:
                        diff = min(diff, np.linalg.norm(i[k] - j[k]))
    return diff


###############Rotate const for search#########################
step_rot = check_diff_rotate(n_y, n_z) * 0.5


# print(step_rot)
# step_rot = 0.0525
#################Checkers#############


def find_section(p0, p1, basis0=np.zeros(2), let_accurance=step_rot, all_posibility=False, method='first', **kwargs):
    '''
    :param p0: this point has already basis
    :param p1: not important basis of p1
    :param basis0: basis of p0-center atom
    :return: section of p0 atom in which there's p1
    '''
    vec = p1 - p0
    vec = np.array([i/np.linalg.norm(vec) for i in vec])
    if method == 'first':
        pp_ = rotate_by_basis(icos, basis0[0], basis0[1], **kwargs)
    elif method == 'incline':
        pp_ = rotate_non_perpendicular(icos, basis0[0], basis0[1], **kwargs)
    elif method == 'ten':
        pp_ = rotate_ten_vars(icos, basis0, **kwargs)
    if all_posibility:
        return min([[np.linalg.norm(ppx - vec), ix] for ix, ppx in enumerate(pp_)])[1]
    while True:
        for num, i in enumerate(pp_):
            if np.linalg.norm(i - vec) <= let_accurance:
                return num
        let_accurance *= 1.1


def get_reversed_section_and_basis(s, b0, method='first', **kwargs):
    if method == 'first':
        return find_basis(np.zeros(3), rotate_by_basis(icos[s], b0[0], b0[1], **kwargs), **kwargs)
    elif method == 'ten':
        return find_basis(np.zeros(3), rotate_ten_vars(icos[s], b0, **kwargs), **kwargs)
    return find_basis(np.zeros(3), rotate_non_perpendicular(icos[s], b0[0], b0[1], **kwargs), **kwargs)


def find_basis(point, connected, method='first', **kwargs):
    '''
    :param point: point for search basis
    :param connected: atoms which have bonds with point
    :return: basis for point (y, z) by min of max different between point and center of section
    '''
    diffs = []
    n_y_l = kwargs.get('n_y', n_y)
    n_z_l = kwargs.get('n_z', n_z)
    if method =='first':
        for j in sorted(product(range(n_y_l), range(n_z_l)), key=lambda x: sum(x)):
            diff = []
            for i in connected:
                v = i - point
                v /= np.linalg.norm(v)
                diff.append(min([np.linalg.norm(v - ppx) for ppx in rotate_by_basis(icos, j[0], j[1], **kwargs)]))
            diffs.append([max(diff), j])
        return min(diffs)[1]
    if method == 'incline':
        for j in sorted(product(range(n_y_l), range(n_z_l)), key=lambda x: sum(x)):
            diff = []
            for i in connected:
                v = i - point
                v /= np.linalg.norm(v)
                diff.append(min([np.linalg.norm(v - ppx) for ppx in
                                     rotate_non_perpendicular(icos, j[0], j[1], **kwargs)]))
            diffs.append([max(diff), j])
    if method == 'ten':
        for j in range(9):
            diff = []
            for i in connected:
                v = i - point
                v /= np.linalg.norm(v)
                diff.append(min([np.linalg.norm(v - ppx) for ppx in
                                 rotate_ten_vars(icos, j, **kwargs)]))
            diffs.append([max(diff), j])
        # return min(diff)[1]
    return min(diffs)[1]

def get_error_of_point(point, method='first', **kwargs):
    basis = find_basis(np.zeros(3), [point], method=method, **kwargs)
    section = find_section(np.zeros(3), point, basis0=basis, method=method, **kwargs)
    if method == 'first':
        pp_ = rotate_by_basis(icos, basis[0], basis[1], **kwargs)
    elif method == 'incline':
        pp_ = rotate_non_perpendicular(icos, basis[0], basis[1], **kwargs)
    elif method == 'ten':
        pp_ = rotate_ten_vars(icos, basis, **kwargs)
    # section = find_section(np.zeros(3), point,method=method)
    return np.linalg.norm(point-pp_[section])
    # if method == 'first':
    #     return np.linalg.norm(point-)
    # elif method == 'ten':
    #     return find_basis(np.zeros(3), rotate_ten_vars(icos[s], b0, **kwargs), **kwargs)
    # return find_basis(np.zeros(3), rotate_non_perpendicular(icos[s], b0[0], b0[1], **kwargs), **kwargs)



if __name__ == '__main__':
    # print(icos[0])
    # for i in rotate_by_basis(icos, 1, 1, n_y=4, n_z=4):
    #     print('\t'.join([str(a) for a in i]))
    # print(get_error_of_point(rotate_by_basis(icos[3], 1, 2, n_y=4, n_z=4), method='first', n_y=4, n_z=4))
    # n_y = 10
    # n_z = 11
    for j in sorted(product(range(1,15), repeat=2), key=lambda x: sum(x)):
        for inner in sorted(product(range(0, j[0]-1), range(0,j[1]-1)), key=lambda x: sum(x)):
            for ix in icos:
                if get_error_of_point(rotate_by_basis(ix, inner[0], inner[1], n_y=j[0], n_z=j[1]), method='first', n_y=j[0], n_z=j[1]) != 0:
                    print(j, inner)

    # print(rotate_non_perpendicular(icos[3], 7, 2, n_y=n_y, n_z=n_z))
    # print(get_error_of_point(rotate_non_perpendicular(icos[5], 7, 10, n_y=n_y, n_z=n_z), method='incline', n_y=n_y, n_z=n_z))
    print(get_error_of_point(rotate_ten_vars(icos[5], 8), method='ten'))
    # print(get_error_of_point(rotate_ten_vars(icos[3], 1), method='ten'))
####################find_section_and_rotate tests#####################################
    # for i in range(12): # p0 - rotated point, p1 - center
    #     for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
    #         if (i not in [0, 6]) or (i in [0, 6] and j[1] == 0):
    #             if find_section_and_rotate(rotate_by_basis(pp[i], j[0], j[1]), np.zeros(3))!=(i,j[0],j[1]):
    #                 print(i, j[0], j[1], find_section_and_rotate(rotate_by_basis(pp[i], j[0], j[1]), np.zeros(3)))
    # sections = [] # check unique
    # for i in range(12):
    #     for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
    #         if (i not in [0, 6]) or (i in [0, 6] and j[1] == 0):
    #             section = get_reversed_section_and_basis(i, j)
    #             if section in sections:
    #                 print(sections)
    #             sections.append(section)
##########################find_only_section_test######################################
    # for i in range(12):
    #     for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
    #         if i in [0, 6] and j[1] == 0:
    #             if not ((i) == find_section(np.zeros(3), rotate_by_basis(pp[i], j[0], j[1])), {'n_y':j[0], 'n_z': j[1]}):
    #                 print(i, j[0], j[1])
    # print(find_basis(np.zeros(3), rotate_by_basis(pp[i], j[0], j[1])), {'n_y':j[0], 'n_z': j[1]})
############################################################################
    # bs = find_basis(np.array([0, 0, 0]), rotate_by_basis(pp[4], 1, 2))
    # print(bs)
    # show_points(pp)
    # show_named_points1({i: icos[i] for i in range(20)})
    # for i in icos:
    #     print(np.linalg.norm(i/np.linalg.norm(i)))
    pass