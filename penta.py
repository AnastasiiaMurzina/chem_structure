from math import acos, sin, pi, cos, atan, asin
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations
import numpy as np

r = 1
a = 4 * r * (10 + 22 / 5 ** 0.5) ** (-0.5)  # length of the side
d = a * sin(3 * pi / 10) / sin(2 * pi / 5)  # distance between the center and the pentagon angle
h = d * sin(3 * pi / 10)  # hight to the middle of the side
angle = acos(- 5 ** (-0.5))  # between two surfaces
center_up = np.array([0, 0, 1])
center_down = np.array([0, 0, -1])
center = np.array([0, 0, 0])
min_distance = h / 2 # for decoding position

def get_side_centers_hotizontal(point):
    initial_angle = 0
    centers = []
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) + point[0], h * sin(initial_angle) + point[1], point[2]]))
        initial_angle += 2 * pi / 5
    return centers

def get_up_surface_centers(point_from):
    initial_angle = 0
    centers = [point_from]
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) * (1 - cos(angle)) + point_from[0], h * sin(initial_angle) * (1 - cos(angle)) + point_from[1], h * sin(angle) + point_from[2]]))
        initial_angle += 2 * pi / 5
    return centers

def get_down_surface_centers(point_from):
    initial_angle = 0
    centers = [point_from]
    for i in range(5):
        centers.append(np.array([h * cos(initial_angle) * (1 - cos(angle)) + point_from[0], h * sin(initial_angle) * (1 - cos(angle)) + point_from[1], -h * sin(angle) + point_from[2]]))
        initial_angle += 2 * pi / 5
    return centers

def show_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in points:
        ax.scatter(i[0], i[1], i[2])
    plt.show()

def get_penta_points(point_center):
    center_up = np.array([i for i in point_center])
    center_up[2] += 1
    center_down = np.array([i for i in point_center])
    center_down[2] -= 1
    return np.concatenate((get_down_surface_centers(center_up), get_up_surface_centers(center_down)), axis=0)

def get_dictionary_coordinates(penta_points):
    dictionary = {}
    for i, j in enumerate(penta_points):
        dictionary.update({i+1: j})
    return dictionary

def section_num_by_coords(center_point, penta_point):
    dictionary = get_dictionary_coordinates(get_penta_points(center_point))
    for key, item in dictionary.items():
        diff = np.array(item) - np.array(penta_point)
        diff = np.linalg.norm(diff)
        if diff < d:
            return key

def horizontal_rotate_points(points, alpha):
    for key, item in points.items():
        x, y, z = item
        x, y = [x * cos(alpha) - y * sin(alpha), x*sin(alpha)+y*cos(alpha)]
        points[key] = np.array([x,y,z])
    return points

#########################angles###########################
def show_points_by_angles(center, angles, labels=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if labels:
        for item, i in angles.items():
            ax.scatter(center[0]+cos(i[0])*sin(i[1]), center[1] + sin(i[0])*sin(i[1]), center[2] + cos(i[1]))
            ax.text(center[0]+cos(i[0])*sin(i[1]), center[1] + sin(i[0])*sin(i[1]), center[2] + cos(i[1]), item)
    else:
        for item, i in angles.items():
            ax.scatter(center[0] + cos(i[0]) * cos(i[1]), center[1] + sin(i[0]) * cos(i[1]), center[2] + sin(i[1]))

    plt.show()

def get_penta_angles():
    angles = [np.array([0,pi/2])]
    init_angle = 0
    for i in range(5):
        init_angle += (2 * pi / 5) % (2 * pi)
        angles.append(np.array([init_angle, asin(sin(angle)/2)]))
    # print(sin(angle)/2)
    angles.append(np.array([0, -pi/2]))
    for i in range(5):
        init_angle += (2 * pi / 5) % (2 * pi)
        angles.append(np.array([init_angle, -asin(sin(angle)/2)]))
    return angles

def get_dictionary_angles(center_point, penta_points):
    dictionary = {}
    dictionary.update({1: np.array([0, 0])})
    alpha = 0
    for i in range(5):
        alpha += 2 * pi / 5
        dictionary.update({i+2: np.array([alpha, pi / 4])})
    alpha = 0
    dictionary.update({7: np.array([0, pi])})
    for i in range(5):
        alpha += 2 * pi / 5
        dictionary.update({i + 8: np.array([alpha, 3 * pi / 4])})
    return dictionary

# BUG HERE! BE careful

def rotate(angles, alpha, betta):
    for key, item in angles.items():
        angles[key] = np.array([(angles[key][0]+alpha)%(2*pi), angles[key][1]+betta if abs(angles[key][1]+betta) < pi/2 else (pi - angles[key][1]-betta)])
    return angles


def horizontal_rotate_angles(angles, alpha):
    for key, item in angles.items():
        pass
        # angles[key] = np.array([angles[key][0]*cos(alpha)])

def vertical_rotate_angles(angles, betta):
    for key, item in angles.items():
        if item[0]>pi:
            v_angle = angles[key][1] + betta
        else:
            v_angle = angles[key][1] - betta
        if abs(v_angle) > pi/2:
            if v_angle>0:
                angles[key][1] = pi - v_angle
            else:
                angles[key][1] = - pi - v_angle
            angles[key][0] = (angles[key][0] + pi) % (2 * pi)
    return angles

########################################################

if __name__ == '__main__':
    my_arr = get_penta_points(center)
    darr = (get_dictionary_coordinates(my_arr))
    show_points([item for key,item in darr.items()])
    # darr = horizontal_rotate_points(darr, pi/2)

    # show_points([item for key,item in darr.items()])

    a2 = get_penta_angles()
    ann_arr = get_dictionary_coordinates(a2)
    print(ann_arr)
    ann_arr = vertical_rotate_angles(ann_arr, pi/6)
    show_points_by_angles(np.array([0,0,0]), ann_arr, labels=True)

    # ann_arr = rotate(ann_arr, 0, pi/12)
    # print(ann_arr)
    # show_points_by_angles(np.array([0,0,0]), ann_arr, labels=True)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for i in a2:
    #     p = np.array([cos(i[0])*cos(i[1]), sin(i[0])*cos(i[1]), sin(i[1])])
    #     ax.scatter(cos(i[0])*cos(i[1]), sin(i[0])*cos(i[1]), sin(i[1]))


    # plt.show()

#######################FALL TEST###############################################
    # zero = np.array([0, 0, 0])
    # for key, item in darr.items():
    #     print(section_num_by_coords(zero, item), section_num_by_coords(item, zero))
#######################################################################################
    # print(section_num_by_coords(np.array([0,0,0]), np.array([0.9,0,-0.5])))
    # print(get_dictionary_angles(center, my_arr))
    # check_zone = np.array([0.29, -0.9, 0.45])
    # key = section_num_by_coords(center, check_zone)
    # print(key)
    # key2 = section_num_by_coords(check_zone, center)
    # print(key2)

    # show_points_by_angles(center, get_dictionary_angles(center, my_arr))
    # distances = []
    # for i, j in combinations(my_arr, 2):
    #     distances.append((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)
    #
    # print(min(distances), max(distances))
