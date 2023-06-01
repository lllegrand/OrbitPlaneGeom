import numpy
import pylab
import copy

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


def get_radiation_pattern_poly(do_plot=True):
    rp = numpy.array([
    # data pulled from Allwave Corporation's test results at 2.48GHz
    # vertical plane
    [-90, -11.5], [-85, -10.5], [-80, -9], [-75,  -8], [-70, -7], [-65, -6], [-60, -5], [-55, -4], [-50, -3], [-45, -2.5], [-40, -2],
    [-35, -1.8], [-30, -1.2], [-25, -1], [-20, -1], [-15, -1], [-10, -1], [-5, -1], [0, -1], [5, -1], [10, -.5], [15, -.1], [20, 0],
    [25, -0.2], [30, -0.5], [35, -1.2], [40, -2.0], [45, -2.8], [50, -4], [55, -5.4], [60, -6.5], [70, -9], [75, -10], [80, -11], [85, -12], [90, -13],
    # horizontal plane
    [-90, -11], [-85, -10.5], [-80, -9.8], [-75, -8.5], [-70, -7.5], [-65, -6.2], [-60, -5], [-55, -3.6], [-50, -2.5], [-45,  -1.7],
    [-40,  -1.3], [-35,  -1.0], [-30,  -0.5], [-25,  -0.2], [-20,  -0.2], [-15,  -0.8], [-10,  -1.0], [-5,  -1.0], [0,  -1], [5,  -1],
    [10,  -0.8], [15,  -1.0], [20,  -1.0], [25,  -1.2], [30,  -1.5], [35,  -1.5], [40,  -1.5], [45,  -1.5], [50,  -2.0],
    [55,  -2.6], [60,  -3.5], [65,  -4.2], [70,  -5.5], [75,  -7.0], [80,  -8.5], [90,  -11.5]])
    #
    xs = abs(rp[:,0])
    ys = rp[:,1]
    if do_plot:
        pylab.figure()
        pylab.plot(xs, ys, '*', label="Measurements")
    #
    zeros = numpy.zeros(len(rp))
    ones = numpy.ones(len(rp))
    #
    x2s = xs*xs
    x3s = x2s*xs
    x4s = x3s*xs
    #
    A = numpy.array([x4s, x3s, x2s, zeros, ones])
    #
    c = numpy.dot(numpy.linalg.pinv(A).T, ys)
    #
    thetas = numpy.linspace(0,90,91)
    yfit = numpy.polyval(c, thetas)
    if do_plot:
        pylab.plot(thetas, yfit, label="Initial Fit")
    #
    # shift up so venter is zero db
    p = copy.copy(c)
    p[-1] = 0.0
    yfit2 = numpy.polyval(p, thetas)
    if do_plot:
        pylab.plot(thetas, yfit2, label="Shifted Fit")
        pylab.grid()
        pylab.xlabel("Off Axis Angle (deg)")
        pylab.ylabel("Gain (dB)")
        pylab.title("Patch Antenna Gain, Measured and Modeled")
        pylab.legend()
    return p


def get_half_track_slant_range_elevation(fov_half_angle, altitude):
    max_fov_half_angle = numpy.arcsin(r_earth/(r_earth + altitude))
    if fov_half_angle >= max_fov_half_angle:
        slant_range, theta = get_bB(C=max_fov_half_angle, c=r_earth, h=altitude)
        min_elevation = 0.0
    else:
        slant_range, theta = get_bB(C=fov_half_angle, c=r_earth, h=altitude)
        min_elevation = numpy.pi/2 - fov_half_angle - theta
    half_track_width = r_earth * numpy.sin(theta)
    return half_track_width, slant_range, min_elevation


def get_slant_range_and_fov(altitude, half_track_width):
    B = numpy.arcsin(half_track_width/r_earth)
    d1 = r_earth * numpy.cos(B)
    d2 = altitude + r_earth - d1
    slant_range = numpy.sqrt(d2**2 + half_track_width**2)
    fov_half_angle = numpy.arctan2(half_track_width, d2)
    return slant_range, fov_half_angle


def get_min_altitude(half_track_width):
    B = numpy.arcsin(half_track_width/r_earth)
    d1 = r_earth * numpy.cos(B)
    d2 = half_track_width * numpy.tan(B)
    min_altitude = d1 + d2 - r_earth
    return min_altitude


def get_altitudes_fovs_dbs(half_track_width, radpoly):
    baseline_slant_range = 600 * 1000
    max_altitude = 2000 * 1000
    #
    min_altitude = get_min_altitude(half_track_width)
    altitudes = numpy.linspace(min_altitude, max_altitude, 5000)
    fov_half_angles = []
    dbs = []
    for altitude in altitudes:
        slant_range, fov_half_angle = get_slant_range_and_fov(altitude, half_track_width)
        db = numpy.polyval(radpoly, numpy.rad2deg(fov_half_angle)) + 20*numpy.log10(baseline_slant_range/slant_range)
        dbs.append(db)
        fov_half_angles.append(fov_half_angle)
    return numpy.array(altitudes), numpy.array(fov_half_angles), numpy.array(dbs)


radpoly = get_radiation_pattern_poly(do_plot=True)


# plot the db tradeoff for a single half_track width
half_track_width = 800000
radpoly = get_radiation_pattern_poly(do_plot=False)
altitudes, fov_half_angles, dbs = get_altitudes_fovs_dbs(half_track_width, radpoly)
#
pylab.figure()
pylab.subplot(2,1,1)
pylab.title("Options for " + str(half_track_width/1000) + "km Half Track Width")
pylab.plot(numpy.rad2deg(fov_half_angles), altitudes/1000)
pylab.grid()
#pylab.xlabel("fov half angle")
pylab.ylabel("altitude (km)")
#
pylab.subplot(2,1,2)
pylab.plot(numpy.rad2deg(fov_half_angles), dbs)
pylab.grid()
pylab.xlabel("fov half angle (deg)")
pylab.ylabel("gain (dB)")




# plot the optimum altitude, fov for a range of half_track widths
half_track_widths = numpy.linspace(400000, 2000000, 100)
radpoly = get_radiation_pattern_poly(do_plot=False)
best_altitudes = []
best_dbs = []
best_fov_half_angles=[]
for half_track_width in half_track_widths:
    altitudes, fov_half_angles, dbs = get_altitudes_fovs_dbs(half_track_width, radpoly)
    idx = numpy.argmax(dbs)
    best_dbs.append(dbs[idx])
    best_altitudes.append(altitudes[idx])
    best_fov_half_angles.append(fov_half_angles[idx])
best_altitudes = numpy.array(best_altitudes)
best_dbs = numpy.array(best_dbs)
#
pylab.figure()
pylab.subplot(3,1,1)
pylab.plot(half_track_widths/1000, best_dbs)
pylab.grid()
pylab.ylabel("relative gain (db)")
pylab.title("Optimum Altitude and FOV for Max Gain")

pylab.subplot(3,1,2)
pylab.plot(half_track_widths/1000, best_altitudes/1000)
pylab.grid()
pylab.ylabel("best altitude (km)")


pylab.subplot(3,1,3)
pylab.plot(half_track_widths/1000, numpy.rad2deg(best_fov_half_angles))
pylab.grid()
pylab.xlabel("half track width (km)")
pylab.ylabel("best fov half angle (deg)")
pylab.ylim((46, 47))

pylab.show()
