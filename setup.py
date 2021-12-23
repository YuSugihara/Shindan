#!/usr/bin/env python3

from distutils.core import setup
from vidkit.__init__ import __version__

setup(name='vidkit',
      version='{}'.format(__version__),
      description='Vid-kit: plant virus diagnosis tool-kit',
      author='Yu Sugihara',
      author_email='yu57th@gmail.com',
      url='https://github.com/YuSugihara/Vid-kit',
      license='GPL',
      packages=['vidkit'],
      entry_points={'console_scripts': [
            'vidkit = vidkit.vidkit:main',
            ]
        }
    )
