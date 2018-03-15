from math import acos, sin, pi, cos, atan, asin
import matplotlib.pyplot as plt
import numpy as np
from itertools import product

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
d_ver_angle = pi / 3
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
    initial_angle = 0
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


def rotate_by_basis(point, y, z):
    '''
    :param point: array or one point to rotate
    :param y: first coordinate of rotation
    :param z: second coordinate of rotation
    :return: rotated point or points
    '''
    operator = Rz(d_hor_angle * z / n_z).dot(Ry(d_ver_angle * y / n_y))
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
            dhs.append(np.array(rotate_by_basis(pp, i, j)))
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
#################Checkers################################
def find_section_and_rotate(p0, p1):
    '''
    :param p0: another point, we'll find section in which this point relatively p1
    :param p1: point for this one we are looking for basis
    :return: section for p1-p0 bond (p0 is in section of p1) and rotated basis of p1
    '''
    wanted = p0 - p1
    sums = []
    available_sections_y = tuple(range(n_y))
    available_sections_z = tuple(range(n_z))
    for num, i in enumerate(pp):
        diff = np.linalg.norm(wanted - i)
        if diff < 1.2 * a / (np.sin(pi / 3)):#h*(2*(1+5**(-0.5)))**0.5/(2*sin(pi/3)):
            diffs = [rotate_by_basis(pp[num], n_y - 1, n_z - 1) - wanted,
                     rotate_by_basis(pp[num], n_y - 1, 0) - wanted,
                     rotate_by_basis(pp[num], 0, n_z - 1) - wanted,
                     rotate_by_basis(pp[num], 0, 0) - wanted]
            s = sum([np.linalg.norm(i) for i in diffs])
            sums.append([s, num])
    # possible problems: 0 and 6 sections
    sums.sort()
    while len(sums):
        ps = sums.pop(0)
        for j in sorted(product(available_sections_y, available_sections_z), key=lambda x: sum(x)):
            if np.linalg.norm(rotate_by_basis(pp[ps[1]], j[0], j[1]) - wanted) < step_rot:
                return ps[1], j[0], j[1]

    # return ps[1], None


def find_section(p0, p1, basis0=np.zeros(2), let_accurance=step_rot):
    '''
    :param p0: this point has already basis
    :param p1: not important basis of p1
    :param basis0: basis of p0-center atom
    :return: section of p0 atom in which there's p1
    '''
    vec = p1 - p0
    pp_ = rotate_by_basis(pp, basis0[0], basis0[1])
    for num, i in enumerate(pp_):
        if np.linalg.norm(i - vec) <= let_accurance:
            return num


def get_reversed_section_and_basis(s, b0):
    return find_section_and_rotate(np.zeros(3), rotate_by_basis(pp[s], b0[0], b0[1]))
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

# def show_named_points(named_points):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     for key, i in named_points.items():
#         ax.scatter(i[0], i[1], i[2])
#         ax.text(i[0], i[1], i[2], str(key))
#     plt.show()

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


########################################################

if __name__ == '__main__':
####################find_section_and_rotate tests#####################################
    for i in range(12): # p0 - rotated point, p1 - center
        for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
            if (i not in [0, 6]) or (i in [0, 6] and j[1] == 0):
                if find_section_and_rotate(rotate_by_basis(pp[i], j[0], j[1]), np.zeros(3))!=(i,j[0],j[1]):
                    print(i, j[0], j[1], find_section_and_rotate(rotate_by_basis(pp[i], j[0], j[1]), np.zeros(3)))
    sections = [] # check unique
    for i in range(12):
        for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
            if (i not in [0, 6]) or (i in [0, 6] and j[1] == 0):
                section = get_reversed_section_and_basis(i, j)
                if section in sections:
                    print(sections)
                sections.append(section)
##########################find_only_section_test######################################
    for i in range(12):
        for j in sorted(product((0, 1, 2, 3), repeat=2), key=lambda x: sum(x)):
            if i in [0, 6] and j[1] == 0:
                if not ((i)==find_section(np.zeros(3), rotate_by_basis(pp[i], j[0], j[1])), [j[0],j[1]]):
                    print(i, j[0], j[1])
############################################################################