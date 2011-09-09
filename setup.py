#!/usr/bin/env python
from distutils.core import setup

setup(name = 'enrico',
      version = '0.1',
      author = 'Christoph Deil, David Sanchez',
      author_email = 'Christoph.Deil@mpi-hd.mpg.de, David.Sanchez@mpi-hd.mpg.de',
      url = 'TODO',
      license = 'BSD',
      description = 'Enrico helps you with your Fermi data analysis',
      long_description=open('README.txt').read(),
      packages = ['enrico', 
                  'enrico.sandbox', 
                  'enrico.user_contrib',
                  ],
      scripts = ['scripts/enrico',
                 ],
      data_files=[('enrico', ['data/default.conf']),
                  ],
      )

