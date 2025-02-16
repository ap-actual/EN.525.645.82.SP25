from math import cos
from math import sin
from math import acos
from math import radians
from math import pi

# coordinates from problem
x0 = [39, -77]
x1 = [39, -100]

# radius of earth, in nmi
re = 3443.8

# get small circle radius via re_small = re * cos(lat)
re_small_circle = re * cos(radians(x0[0]))

# get distance in longitude, in rads
theta = radians(x1[1] - x0[1])

# compute small circle distance
dist_small_circle = re_small_circle * theta
print("---------------\nsmall circle dist (nmi):")
print(dist_small_circle)

# calculate ECEF position of each coordinate
X0 = [0, 0, 0]
X0[0] = re * cos(radians(x0[0])) * cos(radians(x0[1]))
X0[1] = re * cos(radians(x0[0])) * sin(radians(x0[1]))
X0[2] = re * sin(radians(x0[0]))

X1 = [0, 0, 0]
X1[0] = re * cos(radians(x1[0])) * cos(radians(x1[1]))
X1[1] = re * cos(radians(x1[0])) * sin(radians(x1[1]))
X1[2] = re * sin(radians(x1[0]))

# compute great circle distance
x_dot_prod = (X0[0] * X1[0]) + (X0[1] * X1[1]) + (X0[2] * X0[2])
theta      = acos(x_dot_prod / (re * re))

dist_great_circle = theta * re

print("---------------\ngreat circle dist (nmi):")
print(dist_great_circle)
