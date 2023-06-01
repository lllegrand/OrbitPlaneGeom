import orbit_plane_geom 
import pytest

test_args = ['-x']

test_args.append('src/orbit_plane_geom/coverage_model.py')
test_args.append('src/orbit_plane_geom/periodic_binary_signal.py')

pytest.main(test_args)
