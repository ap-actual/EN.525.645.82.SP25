import math


# coordinates from problem
x1 = [39, 77]
x2 = [39, 100]

# radius of earth, in mi
re = 3963.1

# get small circle radius via re_small = re * sin(lat)
re_small_circle = re * math.sin(math.radians(x1[0]))

# get distance in longitude, in rads
theta = math.radians(x2[1] - x1[1])

# compute small circle distance
dist_small_circle = re_small_circle * theta
print(dist_small_circle)

# calculate ECEF position of each coordinate
X1_X = re * math.cos(math.radians(x1[0])) * math.cos(math.radians(x1[1]))
X1_Y = re * math.cos(math.radians(x1[0])) * math.sin(math.radians(x1[1]))
X1_Z = re * math.sin(math.radians(x1[0]))

X2_X = re * math.cos(math.radians(x2[0])) * math.cos(math.radians(x2[1]))
X2_Y = re * math.cos(math.radians(x2[0])) * math.sin(math.radians(x2[1]))
X2_Z = re * math.sin(math.radians(x2[0]))

# compute great circle distance
x_dot_prod = (X1_X * X2_X) + (X1_Y * X2_Y) + (X1_Z * X2_Z)
theta      = math.acos(x_dot_prod / (re * re))

dist_great_circle = theta * re
print(dist_great_circle)
