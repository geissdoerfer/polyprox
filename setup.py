from setuptools import setup, Extension

mod_algs = Extension('algorithms', sources = ['src/algorithms.c'])

setup(name='polyprox',
      version='0.1',
      description='Polygonal curve approximation tools',
      author='Kai Geissdoerfer',
      author_email='kai.geissdoerfer@mailbox.tu-berlin.de',
      url = 'https://github.com/geissdoerfer/polyprox',
      download_url = 'https://github.com/geissdoerfer/polyprox/archive/0.1.tar.gz',
      ext_modules= [mod_algs],
      packages=['polyprox'],
      install_requires=['numpy']
     )
