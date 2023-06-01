import numpy
import pylab

r_earth=6378136.999954619 #radius_equator

def get_bB(C, c, h):
    # imagine:
    #   a circle centered on point pb with radius c
    #   a point pc that is at distance a=c+h away from the center of the circle
    #   a point pa that is on the circle
    #   with the angle pa-to-pc-to-pb equal to C
    # find:
    #   the distance pa-to-pc, i.e. b
    #   the angle pa-to-pb-to-pc, i.e. B
    # upper case vars are angles, lower case vars are lengths of the sides opposite the angle
    a = c+h
    # law of cosines:
    # cos A = (b^2 + c^2 - a^2)/(2*b*c)
    # cos C = (a^2 + b^2 - c^2)/(2*a*b)
    # law of sines:
    # sin(A)/a = sin(B)/b = sin(C)/c
    A = numpy.pi - numpy.arcsin(a * numpy.sin(C) / c) # law of sines for obtuse triangles
    b = c * numpy.cos(A) + a*numpy.cos(C) # proof via law of cosines
    B = numpy.arcsin(b * numpy.sin(C)/c) # law of sines
    return b, B


def htw_sr_me_from_alt_fov(altitude, fov_half_angle_deg):
    max_fov_half_angle_deg = numpy.rad2deg(numpy.arcsin(r_earth/(r_earth + altitude)))
    if fov_half_angle_deg >= max_fov_half_angle_deg:
        slant_range, theta = get_bB(C=numpy.deg2rad(max_fov_half_angle_deg), c=r_earth, h=altitude)
        min_elevation = 0.0
    else:
        slant_range, theta = get_bB(C=numpy.deg2rad(fov_half_angle_deg), c=r_earth, h=altitude)
        min_elevation = 90 - fov_half_angle_deg - numpy.rad2deg(theta)
    half_track_width = r_earth * numpy.sin(theta)
    return half_track_width, slant_range, min_elevation




nsamps = 50
fov_half_angles = numpy.linspace(45, 70, nsamps)
altitudes = numpy.linspace(400000, 1000000, nsamps)

X, Y = numpy.meshgrid(fov_half_angles, altitudes)
(rows, cols) = numpy.shape(X)

half_track_width_list=[]
slant_range_list=[]
min_elevation_list=[]
for row in range(rows):
    half_track_width_row=[]
    slant_range_row = []
    min_elevation_row = []
    for col in range(cols):
        fov_half_angle = X[row,col]
        altitude = Y[row,col]
        max_fov_half_angle = numpy.rad2deg(numpy.arcsin(r_earth/(r_earth + altitude)))
        if fov_half_angle >= max_fov_half_angle:
            slant_range, theta = get_bB(C=numpy.deg2rad(max_fov_half_angle), c=r_earth, h=altitude)
            min_elevation = 0.0
        else:
            slant_range, theta = get_bB(C=numpy.deg2rad(fov_half_angle), c=r_earth, h=altitude)
            min_elevation = 90 - fov_half_angle - numpy.rad2deg(theta)
        half_track_width = r_earth * numpy.sin(theta)
        half_track_width_row.append(half_track_width/1000) # km
        slant_range_row.append(slant_range/1000) # km
        min_elevation_row.append(min_elevation)
    half_track_width_list.append(half_track_width_row)
    slant_range_list.append(slant_range_row)
    min_elevation_list.append(min_elevation_row)
half_track_width=numpy.array(half_track_width_list)
slant_range=numpy.array(slant_range_list)
min_elevation=numpy.array(min_elevation_list)
Y = Y/1000 # km

f=pylab.figure()
ax1 = pylab.subplot(111)
#
c_htws = ax1.contour(X, Y, half_track_width, levels=[200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200], colors='r')
ax1.clabel(c_htws, fontsize=9, inline=1, fmt='%3.1f')
ax1.plot([X[0][0], X[0][0]], [Y[0][0], Y[0][0]], 'r', label="Half Track Width (km)")
#
c_srs = ax1.contour(X, Y, slant_range,
                       levels=[200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200],
                       colors='g')
ax1.clabel(c_srs, fontsize=9, inline=1, fmt='%3.1f')
ax1.plot([X[0][0], X[0][0]], [Y[0][0], Y[0][0]], 'g', label="Slant Range (km)")
#
c_mes = ax1.contour(X, Y,
                    min_elevation,
                    levels=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60],
                    colors='b')
ax1.clabel(c_mes, fontsize=9, inline=1, fmt='%2.1f')
ax1.plot([X[0][0], X[0][0]], [Y[0][0], Y[0][0]], 'b', label="Min Elevation (deg)")
#

#ax1.legend(loc='upper right')
ax1.legend()
#
ax1.grid(True)
ax1.set_xlabel("fov half angle (deg)")
ax1.set_ylabel("altitude (km)")
ax1.set_title("Contours of Half Track Width, Slant Range, and Min Elevation")
pylab.show()

    
