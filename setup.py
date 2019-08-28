from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        __builtins__.__NUMPY_SETUP__ = False
        import numpy

        self.include_dirs.append(numpy.get_include())


mod_algs = Extension("algorithms", sources=["src/algorithms.c"])
setup(
    name="polyprox",
    version="0.4",
    description="Polygonal curve approximation tools",
    author="Kai Geissdoerfer",
    author_email="kai.geissdoerfer@mailbox.tu-berlin.de",
    url="https://github.com/geissdoerfer/polyprox",
    ext_modules=[mod_algs],
    packages=["polyprox"],
    setup_requires=["numpy", "pytest-runner"],
    cmdclass={"build_ext": build_ext},
    tests_require=["pytest>=3.9", "rdp"],
)
