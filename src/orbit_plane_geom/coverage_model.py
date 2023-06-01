import numpy
import geom3d

class CoverageModel(object):
    '''
    A standalone computational geometry model for orbital planes
    and coveage"
    '''


    r_earth=6378137.0 #radius_equator


    def __init__(self, inclination_deg, n_planes, dar_angle_deg, half_track_width):
        #
        self.inclination_deg = inclination_deg
        self.n_planes = n_planes
        self.dar_angle_deg = dar_angle_deg
        self.half_track_width = half_track_width
        #
        self.earth = geom3d.Sphere(center=[0,0,0], radius=CoverageModel.r_earth)
        self.orbital_planes = CoverageModel.get_orbit_planes(inclination_deg, n_planes, dar_angle_deg)


    @classmethod
    def from_altitude_angle(cls, inclination_deg, n_planes, dar_angle_deg, altitude, fov_pyramid_face_half_angle):
        half_track_width = CoverageModel.get_half_track_width(CoverageModel.r_earth, altitude, fov_pyramid_face_half_angle)
        result = cls(inclination_deg, n_planes, dar_angle_deg, half_track_width)
        result.altitude = altitude
        result.fov_pyramid_face_half_angle = fov_pyramid_face_half_angle
        return result


    def get_coverage_count_edges(self, latitude):
        ecef_xyz = CoverageModel.get_ecef_xyz(self.earth, latitude, longitude=-180)
        initial_count = CoverageModel.get_coverage_count(self.earth, ecef_xyz, self.orbital_planes, self.half_track_width)
        longitudes, events = CoverageModel.get_events(self.earth, latitude, self.orbital_planes, self.half_track_width)
        if len(longitudes) == 0:
            assert len(events) == 0
            longitudes = numpy.array([-180])
            counts = numpy.array([initial_count])
        else:
            counts = initial_count + numpy.cumsum(events)
        return longitudes, counts


    def get_coverage_gaps(self, latitude):
        longitudes, counts = self.get_coverage_count_edges(latitude)
        gaps = []
        prev_long = longitudes[-1]-360.
        prev_count = counts[-1]
        for long, count in zip(longitudes, counts):
            if prev_count == 0:
                gaps.append(long-prev_long)
            prev_long = long
            prev_count = count
        return gaps


    def is_full_coverage(self, latitude):
        gaps = self.get_coverage_gaps(latitude)
        if len(gaps) == 0:
            return True
        elif numpy.all(numpy.array(gaps) < 2*numpy.pi/(24*3600)):
            return True
        else:
            return False


    def get_max_full_coverage_latitude(self, lat_step=1.0):
        # early out
        latitude = -lat_step
        while latitude+lat_step <= 90 and self.is_full_coverage(latitude+lat_step):
            latitude = latitude+lat_step
        return latitude


    ####
    # utils
    ####


    def get_orbit_planes(inclination_deg, n_planes, dar_angle_deg):
        dar_angle = numpy.deg2rad(dar_angle_deg)
        inclination = numpy.deg2rad(inclination_deg)
        planes = []
        for i_plane in range(n_planes):
            theta = i_plane*dar_angle/n_planes
            planes.append(geom3d.Plane(point=[0,0,0], normal=[numpy.sin(inclination)*numpy.cos(theta), numpy.sin(inclination)*numpy.sin(theta), numpy.cos(inclination)]))
        return planes


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


    def get_half_track_width(r_earth, altitude, fov_pyramid_face_half_angle):
        slant_range, arc_angle = CoverageModel.get_bB(C=fov_pyramid_face_half_angle, c=r_earth, h=altitude)
        half_track_width = r_earth*numpy.sin(arc_angle)
        return half_track_width


    def minus_pi_to_pi(theta):
        while theta < -numpy.pi:
            theta = theta + 2*numpy.pi
        while theta >= numpy.pi:
            theta = theta - 2*numpy.pi
        assert -numpy.pi <= theta and theta < numpy.pi
        return theta


    def get_ecef_xyz(earth, latitude, longitude):
        lat_rad = numpy.deg2rad(latitude)
        long_rad = numpy.deg2rad(longitude)
        line = geom3d.Line(point=[0,0,0],
               direction=[numpy.cos(lat_rad)*numpy.cos(long_rad), numpy.cos(lat_rad)*numpy.sin(long_rad), numpy.sin(lat_rad)])
        xyz = earth.line_intersection(line)[0]
        return xyz


    def get_coverage_count(earth, xyz, orbital_planes, half_track_width):
        # count how many times xyz is between the edge planes
        count = 0
        for plane in orbital_planes:
            a_plane = geom3d.Plane.from_offset(plane, -half_track_width)
            a = numpy.dot(xyz - a_plane.point, a_plane.normal)
            b_plane = geom3d.Plane.from_offset(plane, half_track_width)
            b = numpy.dot(xyz - b_plane.point, b_plane.normal)
            if a*b <= 0:
                count += 1
        return count


    def get_events(earth, latitude, orbital_planes, half_track_width):
        # The coverage of a given orbital plane is bounded my two parallel planes
        # half_track_width on either side of the orbital plane.
        # this function takes a constant latitude slice of this coverage
        # Returns:
        # longitudes: a sorted increasing list of longitudes of coverage zone edges
        # events: for each edge_lon, the incremental coverage count
        #
        lat_rad = numpy.deg2rad(latitude)
        latitude_plane = geom3d.Plane.from_abcd(abcd=[0,0,1,-CoverageModel.r_earth*numpy.sin(lat_rad)])
        longitudes = numpy.array([])
        events = numpy.array([])
        #
        for plane in orbital_planes:
            a_plane = geom3d.Plane.from_offset(plane, -half_track_width)
            a_line = a_plane.plane_intersection(latitude_plane)
            a_points = earth.line_intersection(a_line)
            #
            b_plane = geom3d.Plane.from_offset(plane, half_track_width)
            b_line = b_plane.plane_intersection(latitude_plane)
            b_points = earth.line_intersection(b_line)
            #
            if not a_points is None:
                a_lons = numpy.arctan2(a_points[:,1], a_points[:,0])
                longitudes = numpy.append(longitudes, numpy.rad2deg(CoverageModel.minus_pi_to_pi(a_lons[0])))
                events = numpy.append(events, 1)
                longitudes = numpy.append(longitudes, numpy.rad2deg(CoverageModel.minus_pi_to_pi(a_lons[1])))
                events = numpy.append(events, -1)
            if not b_points is None:
                b_lons = numpy.arctan2(b_points[:,1], b_points[:,0])
                longitudes = numpy.append(longitudes, numpy.rad2deg(CoverageModel.minus_pi_to_pi(b_lons[0])))
                events = numpy.append(events, -1)
                longitudes = numpy.append(longitudes, numpy.rad2deg(CoverageModel.minus_pi_to_pi(b_lons[1])))
                events = numpy.append(events, 1)
        # consolidate repeats
        uniq_longs = list(set(longitudes))
        if len(uniq_longs) < len(longitudes):
            uniq_events = []
            for long in uniq_longs:
                idxs = numpy.where(longitudes == long)[0]
                uniq_events.append(sum(events[idxs]))
            longitudes = numpy.array(uniq_longs)
            events = numpy.array(uniq_events)
        # sort
        if len(longitudes) > 0:
            order =  numpy.argsort(longitudes)
            longitudes = numpy.array(longitudes)[order]
            events = numpy.array(events)[order]
        #
        assert len(longitudes) == len(events)
        assert sum(events) == 0
        return longitudes, events


    def count_edges_to_curve(longitudes, counts):
        assert len(longitudes) == len(counts)
        longitudes_curve = [-180]
        prev_count = counts[-1]
        counts_curve = [prev_count]
        for longitude, count in zip(longitudes, counts):
            longitudes_curve.append(longitude)
            counts_curve.append(prev_count)
            longitudes_curve.append(longitude)
            counts_curve.append(count)
            prev_count = count
        longitudes_curve.append(180)
        counts_curve.append(prev_count)
        return longitudes_curve, counts_curve


    def get_min_half_track_width_for_full_coverage(inclination_deg, n_planes, dar_angle_deg, max_latitude_of_interest):
        def get_lat(width):
            model = CoverageModel(inclination_deg, n_planes, dar_angle_deg, width)
            return model.get_max_full_coverage_latitude()
        lo_width = 0
        assert get_lat(lo_width) < max_latitude_of_interest
        hi_width = CoverageModel.r_earth
        assert get_lat(hi_width) > max_latitude_of_interest
        while hi_width - lo_width > 0.1:
            mid_width = (hi_width + lo_width)/2
            if get_lat(mid_width) < max_latitude_of_interest:
                lo_width = mid_width
            else:
                hi_width = mid_width
        return mid_width


    def get_best_dar(inclination_deg, half_track_width, n_planes, lat_step, dar_resolution):
        dar = 180.
        dar_step = 8 
        best_dar = None
        while dar_step >= dar_resolution:
            model = CoverageModel(inclination_deg, n_planes, dar, half_track_width)
            fcl = model.get_max_full_coverage_latitude(lat_step)
            print(dar, dar_step, fcl)
            while fcl >= 0 and dar <= 360-dar_step:
                best_dar = dar        
                dar = dar+dar_step
                model = CoverageModel(inclination_deg, n_planes, dar, half_track_width)
                fcl = model.get_max_full_coverage_latitude(lat_step)
                print(dar, dar_step, fcl)
            assert not best_dar is None, "No coverage at any dar"
            dar_step = dar_step/2
            dar = best_dar + dar_step
        model = CoverageModel(inclination_deg, n_planes, best_dar, half_track_width)
        best_fcl = model.get_max_full_coverage_latitude(lat_step)
        print(best_dar, best_fcl)
        return best_dar


    def get_full_coverage_latitude_contour_data(inclination_deg, n_planes, wrange=None, nsamps=10, lat_step=1.0):
        dars = numpy.linspace(180, 360, nsamps)
        #
        if wrange is None:
            inclination = numpy.deg2rad(inclination_deg)
            min_full_coverage_half_track_width = CoverageModel.r_earth*numpy.sin(numpy.pi/(2*n_planes))*numpy.sin(inclination) + 10
            widths = numpy.linspace(min_full_coverage_half_track_width, 2*min_full_coverage_half_track_width, nsamps)
        else:
            widths = numpy.linspace(wrange[0], wrange[1], nsamps)
        #
        #def get_full_coverage_latitude(self):
        X, Y = numpy.meshgrid(dars, widths)
        (rows, cols) = numpy.shape(X)
        fcl_list=[]
        for row in range(rows):
            fcl_row = []
            for col in range(cols):
                dar = X[row,col]
                width = Y[row,col]
                model = CoverageModel(inclination_deg=inclination_deg,
                                            n_planes=n_planes,
                                            dar_angle_deg=dar,
                                            half_track_width=width)
                fcl = model.get_max_full_coverage_latitude(lat_step)
                fcl_row.append(fcl)
            fcl_list.append(fcl_row)
            print(row)
        fcl = numpy.array(fcl_list)
        return X, Y, fcl


    def do_full_coverage_latitude_contour_plot(inclination_deg, n_planes, nsamps=50, lat_steps=1.0):
        import pylab

        X, Y, fcl = CoverageModel.get_full_coverage_latitude_contour_data(inclination_deg, n_planes, wrange=None, nsamps=nsamps, lat_step=1.0)

        f=pylab.figure()
        ax1 = pylab.subplot(111)
        #
        max_fcl = max(fcl.flatten())
        levels = list(numpy.arange(0, max_fcl, 10)) + [max_fcl]
        c_fcls = ax1.contour(X, Y/1000, fcl,
                             levels=levels)
         #                    colors='r')
        ax1.clabel(c_fcls, fontsize=9, inline=1, fmt='%2.1f')
        #
        lines = [c_fcls.collections[0]]
        labels = ['Full Coverage Latitude']
        #ax1.legend(lines,labels,loc='lower right')

        #
        ax1.grid(True)
        ax1.set_xlabel("Dar Angle (deg)")
        ax1.set_ylabel("Half Track Width (km)")
        title = "Full Coverage Latitudes, Inclination: " + str(inclination_deg) + " deg, Planes: " + str(n_planes)
        ax1.set_title(title)
        #
        filename = "inc"+str(inclination_deg)+"n"+str(n_planes)
        pylab.savefig(filename)
        #pylab.show()

#######
### Demos
#######

def do_coverage_model_visualize():
    inclination_deg = 45
    n_planes = 4
    dar_angle_deg = 180
    half_track_width=0

    my_model = CoverageModel(inclination_deg,
                                    n_planes,
                                    dar_angle_deg,
                                    half_track_width)
    earth = my_model.earth
    latitude = 30
    orbital_planes = my_model.orbital_planes

    vis = geom3d.Visualizer()

    vis.add_actor(earth)

    ecef_xyz = CoverageModel.get_ecef_xyz(earth, latitude, longitude=-180)
    vis.add_points([ecef_xyz])

    for plane in orbital_planes:
         circle = earth.plane_intersection(plane)
         vis.add_actor(circle)
         #
         plane1 = geom3d.Plane.from_offset(plane, half_track_width)
         circle1 = earth.plane_intersection(plane1)
         vis.add_actor(circle1)
         #
         plane2 = geom3d.Plane.from_offset(plane, -half_track_width)
         circle2 = earth.plane_intersection(plane2)
         vis.add_actor(circle2)
    vis.show()


def do_coverage_model_plot_coverage_counts():
    import pylab

    inclination_deg = 90
    n_planes = 3
    dar_angle_deg = 200
    #half_track_width=600000
    half_track_width = CoverageModel.r_earth*numpy.sin(numpy.pi/(2*n_planes))*numpy.sin(numpy.deg2rad(inclination_deg)) + 10

    my_model = CoverageModel(inclination_deg,
                                    n_planes,
                                    dar_angle_deg,
                                    half_track_width)

    latitude = 0
    gaps = my_model.get_coverage_gaps(latitude)
    print("gaps:", gaps)
    longitudes, counts = my_model.get_coverage_count_edges(latitude)

    xs,ys = CoverageModel.count_edges_to_curve(longitudes, counts)
    pylab.plot(xs, ys)
    pylab.show()


def do_coverage_model_fcl_contour_plots():
    #for inclination_deg in [30, 40, 50, 60, 70, 80, 90, 100]:
    #for inclination_deg in [30, 40]:
    #for inclination_deg in [50, 60]:
    #for inclination_deg in [70, 80]:
    for inclination_deg in [90, 100]:
        for n_planes in [6,7,8,9,10,11,12,13,14,15,16]:
            CoverageModel.do_full_coverage_latitude_contour_plot(inclination_deg, n_planes, nsamps=50, lat_steps=1.0)
            print(inclination_deg, n_planes)


if __name__ == "__main__":
    #do_coverage_model_visualize()
    do_coverage_model_plot_coverage_counts()
    #do_coverage_model_fcl_contour_plots()
