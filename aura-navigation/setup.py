
#!/usr/bin/python3

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.sysconfig import customize_compiler

import sys
import os
os.environ["CC"] = "g++-9"
os.environ["CXX"] = "g++-9"

# this class removes warning of 
# cc1plus: warning: command line option “-Wstrict-prototypes” is valid for Ada/C/ObjC but not for C++
# see: https://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
class my_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)

# class get_pybind_include(object):
#     """Helper class to determine the pybind11 include path
#     The purpose of this class is to postpone importing pybind11
#     until it is actually installed, so that the ``get_include()``
#     method can be invoked. """

#     def __str__(self):
#         import pybind11
#         return pybind11.get_include()

setup(name="aurauas_navigation",
      version="1.4",
      description="Navigation Tools",
      author="Curtis L. Olson",
      author_email="curtolson@flightgear.org",
      url="https://github.com/AuraUAS",
      #py_modules=["props", "props_json", "props_xml"],
      #package_dir = {"": "lib"},
      packages = find_packages(),
      cmdclass = {'build_ext': my_build_ext},
      setup_requires=['numpy'],
      ext_modules=[
          Extension("aurauas_navigation.structs",
                    include_dirs=["/usr/local/include/"], # eigen3
                    define_macros=[("HAVE_PYBIND11", "1")],
                    #extra_compile_args=["-std=c++11"],
                    sources=["src/nav_common/structs.cpp"],
                    depends=["src/nav_common/structs.h"]),
        #   Extension("aurauas_navigation.ekf15",
        #             extra_compile_args=["-std=c++11"],
        #             include_dirs=["/usr/local/include/eigen3/"],
        #             sources=["src/nav_ekf15/pybind11.cpp",
        #                      "src/nav_ekf15/EKF_15state.cpp",
        #                      "src/nav_common/nav_functions.cpp"],
        #             depends=["src/nav_ekf15/EKF_15state.h",
        #                      "src/nav_common/constants.h",
        #                      "src/nav_common/nav_functions.h"]),
        #   Extension("aurauas_navigation.ekf15_mag",
        #             #extra_compile_args=["-std=c++11"],
        #             include_dirs=["/usr/local/include/"],
        #             sources=["src/nav_ekf15_mag/pybind11.cpp",
        #                      "src/nav_ekf15_mag/EKF_15state.cpp",
        #                      "src/nav_common/nav_functions.cpp",
        #                      "src/nav_common/coremag.cc"],
        #             depends=["src/nav_ekf15_mag/EKF_15state.h",
        #                      "src/nav_common/constants.h",
        #                      "src/nav_common/nav_functions.h",
        #                      "src/nav_common/coremag.h"]),
          Extension("aurauas_navigation.ekf17",
                    # extra_compile_args=["-std=c++11"],
                    include_dirs=["/usr/local/include/"],
                    sources=["src/nav_tightly_ekf17/pybind11.cpp",
                             "src/nav_tightly_ekf17/nav-functions.cpp",
                             "src/nav_tightly_ekf17/uNavINS.cpp"],
                    depends=["src/nav_tightly_ekf17/nav-functions.h",
                             "src/nav_tightly_ekf17/uNavINS.h"]),                
          Extension("aurauas_navigation.openloop",
                    #language="c++",
                    #extra_compile_args=["-std=c++11"],
                    #link_flags = ['-fPIC', '-lstdc++'],
                    #extra_compile_args = ['-std=c++11', '-fPIC'],
                    include_dirs=["/usr/local/include/"],
                    sources=["src/nav_openloop/pybind11.cpp",
                             "src/nav_openloop/openloop.cpp",
                             "src/nav_openloop/glocal.cpp",
                             "src/nav_common/nav_functions.cpp",
                             "src/nav_common/coremag.cc"],
                    depends=["src/nav_openloop/openloop.h",
                             "src/nav_openloop/glocal.h",
                             "src/nav_common/constants.h",
                             "src/nav_common/nav_functions.h",
                             "src/nav_common/coremag.h"]),
          Extension("aurauas_navigation.uNavINS",
                    #extra_compile_args=["-std=c++11"],
                    include_dirs=["/usr/local/include/"],
                    sources=["src/UASLab_RAPTRS/pybind11.cpp",
                             "src/UASLab_RAPTRS/nav-functions.cpp",
                             "src/UASLab_RAPTRS/uNavINS.cpp"],
                    depends=["src/UASLab_RAPTRS/nav-functions.h",
                             "src/UASLab_RAPTRS/uNavINS.h"]),
        #   Extension("aurauas_navigation.uNavINS_BFS",
        #             extra_compile_args=["-std=c++11"],
        #             include_dirs=["/usr/local/include/eigen3/"],
        #             sources=["src/BFS_raptrs/pybind11.cpp",
        #                      "src/BFS_raptrs/uNavINS.cpp"],
        #             depends=["src/BFS_raptrs/uNavINS.h"])
      ],
     )

