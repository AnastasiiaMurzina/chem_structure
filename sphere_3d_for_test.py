import matplotlib.pyplot as plt
import numpy as np
import cmath
from cmath import pi
from mpl_toolkits.mplot3d import Axes3D
def dec_point_to_sph(point):
    return [cmath.atan(point[1]/point[0]), cmath.acos(point[2]/cmath.sqrt(point[0]**2+point[1]**2+point[2]**2))]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
zero_point = [0,0,0]
ax.scatter(zero_point[0], zero_point[1], zero_point[2])
first_point = [1,0,0]
ax.scatter(first_point[0], first_point[1], first_point[2])
sfp = dec_point_to_sph(first_point)
n=2
if n==2:
    angle_elevation = 0+sfp[1]
    angle_azimuth = pi + sfp[0]
    endx = zero_point[0] + cmath.cos(angle_azimuth) * cmath.sin(angle_elevation)
    endy = zero_point[0]+cmath.sin(angle_azimuth)*cmath.sin(angle_elevation)
    endz = zero_point[0]+cmath.cos(angle_elevation)
    print(endx.real, endy.real, endz.real)
    print(angle_elevation*180/pi, angle_azimuth*180/pi)
    ax.scatter(endx.real, endy.real, endz.real)
# for i in range(n):
#     angle_elevation = (i)*pi/(n) - pi/2
#     angle_azimuth = (i)*2*pi/(n)
#     if n//2==0:
#         angle_azimuth*=2
#         angle_elevation*=2
#     # print(angle_elevation, angle_azimuth)
#     endx = zero_point[0]+cmath.cos(angle_azimuth)*cmath.sin(angle_elevation)
#     endy = zero_point[0]+cmath.sin(angle_azimuth)*cmath.sin(angle_elevation)
#     endz = zero_point[0]+cmath.cos(angle_elevation)
#     print(endx.real, endy.real, endz.real)
#     print(angle_elevation*180/pi, angle_azimuth*180/pi)
#     ax.scatter(endx.real, endy.real, endz.real)
for i in range(n):
    angle1 = i*360/n
    angle2 = i*180/n - 90
    print(angle1,angle2)
fig.show()
plt.show()
