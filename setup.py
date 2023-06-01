from setuptools import setup, find_packages

setup(
    name='OrbitPlaneGeom',
    version='1.0',
    author='louis',
    description='Satellite constellation coverage via computational geometry',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'numpy',
        'matplotlib',
        'geom3d @ git+ssh://git@github.com/HubbleNetwork/Geom3D.git'
    ]
)
