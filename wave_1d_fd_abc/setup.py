#!/usr/bin/env python
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('wave_1d_fd_abc', parent_package, top_path)
    config.add_extension(name='oneway', sources=['oneway.f90'], extra_f90_compile_args=['-Wall', '-Wextra', '-pedantic', '-Werror', '-march=native', '-O2', '-std=f95', '-fbounds-check'])
    config.add_extension(name='oneway_plain', sources=['oneway_plain.f90'], extra_f90_compile_args=['-Wall', '-Wextra', '-pedantic', '-Werror', '-march=native', '-O2', '-fbounds-check'])
    config.add_extension(name='twoway_plain', sources=['twoway_plain.f90'], extra_f90_compile_args=['-Wall', '-Wextra', '-pedantic', '-Werror', '-march=native', '-O2', '-fbounds-check'])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
