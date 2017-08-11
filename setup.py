#!/usr/bin/env python
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('wave_1d_fd_abc', parent_package, top_path)
    config.add_subpackage('wave_1d_fd_abc')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
