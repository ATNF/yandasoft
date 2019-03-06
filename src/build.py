# @file
# build script for AutoBuild

from askapdev.rbuild.builders import Scons as Builder

b = Builder(".")
b.build()

