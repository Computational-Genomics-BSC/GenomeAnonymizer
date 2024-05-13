from setuptools import setup, find_packages, Extension
# If you have .pyx things to cythonize
from Cython.Build import cythonize
from pysam import get_include as pysam_get_include
from pysam import get_defines as pysam_get_defines


"""setup(
    ext_modules=cythonize(["anonymizer_methods.py",
                           "short_read_tumor_normal_anonymizer.py"],
                          language_level="3"),
    # compiler_directives={'language_level': "3"}
)"""
extensions = [
    Extension(
        "pileup_io", ["src/GenomeAnonymizer/pileup_io.pyx"],
        include_dirs=pysam_get_include(),
        define_macros=pysam_get_defines()),
    ]
setup(
    ext_modules=cythonize(extensions,
                          language_level="3"),
    # compiler_directives={'language_level': "3"}
)