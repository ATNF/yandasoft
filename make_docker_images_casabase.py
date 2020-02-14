#!/usr/bin/env python3
# 
# Creates base image containing Casacore and Casarest,
# with specified MPI implementations.
# Note that -i (or --image) option must be used to actually create the image.
# Otherwise, only the Dockerfile is created.
# 
# Author: Paulus Lahur
#
#------------------------------------------------------------------------------
# SETTINGS
#
# Docker images will be created for all:
# - Specific machine targets
# - All MPI targets in generic machine target
#
# Note that, for specific machine, MPI target is already specified.
# Example:
# machine_targets = ["generic", "galaxy"]
machine_targets = ["galaxy"]

# MPI implementations for "generic" machine.
# A valid mpi target is either "mpich" or "openmpi-X.Y.Z", 
# where X, Y and Z are version numbers (major, minor and revision).
# Example:
# mpi_targets = ["mpich", "openmpi-2.1.6", "openmpi-3.1.4", "openmpi-4.0.2"]
mpi_targets = []

# Git repository of Yandasoft
git_repository = "https://github.com/ATNF/yandasoft.git"

# Docker image name is in this format: 
# target_prepend + (specific machine target OR mpi target for generic machine) + target_append
#image_prepend = "lahur/yandasoft-"
image_prepend = "lahur/casabase-"
image_append = ":latest"

# Set True if this is just a dry run. No Docker image will be created.
dry_run = True

#------------------------------------------------------------------------------
# CODE

import sys
import argparse
import subprocess

# Sanitizing parameters

machine_targets = list(map(str.lower, machine_targets))
mpi_targets = list(map(str.lower, mpi_targets))


class DockerClass:
    def set_file_name(self, file_name):
        self.file_name = file_name

    def set_content(self, content):
        self.content = content

    def set_image(self, image):
        self.image = image

    def write(self):
        '''Write dockerfile'''
        f = open(self.file_name, "w")
        f.write(self.content)
        f.close()


def get_mpi_type_and_version(mpi_name):
    '''
    Given the full name of MPI, return the MPI type: mpich or openmpi
    as well as the version.
    '''
    if (len(mpi_name) > 5):
        if (mpi_name[0:5] == "mpich"):
            return ("mpich", mpi_name[6:])
        elif (mpi_name[0:8] == "openmpi-"):
            return ("openmpi", mpi_name[8:])
        else:
            raise ValueError("Illegal MPI name", mpi_name)
    elif (len(mpi_name) == 5):
        if (mpi_name == "mpich"):
            return("mpich", "")
        else:
            raise ValueError("Illegal MPI name", mpi_name)
    else:
        raise ValueError("Illegal MPI name", mpi_name)


def get_mpi_type(mpi_name):
    '''
    Given the full name of MPI implementation, return the type: 
    MPICH, OpenMPI or None
    '''
    if (mpi_name == "mpich"):
        return "mpich"
    elif (mpi_name[0:8] == "openmpi-"):
        return "openmpi"
    else:
        print("ERROR: Illegal MPI name: ", mpi_name)
        return None


def get_openmpi_version(mpi_name):
    '''
    Given the full name of MPI implementation, return OpenMPI version.
    Return None if not OpenMPI.
    '''
    if (mpi_name[0:8] == "openmpi-"):
        return mpi_name[8:]
    else:
        print("ERROR: This is not OpenMPI: ", mpi_name)
        return None


def main():
    '''
    The main code does 3 things:
    1) Create Dockerfiles for all MPI implementations.
    2) If required, create the corresponding Docker images.
    3) Create sample SLURM batch files.
    '''

    parser = argparse.ArgumentParser(description="Make Docker images for various MPI implementations")
    parser.add_argument('-i', '--image', help='Create Docker images', action='store_true')
    parser.add_argument('-s', '--slurm', help='Create sample batch files for SLURM', action='store_true')
    args = parser.parse_args()

    print("Making Dockerfiles for all targets ...")

    docker_targets = []

    header = ("# This file is automatically created by " + __file__ + "\n")

    cmake_ver = "3.15.7"
    cmake_source = "cmake-" + cmake_ver + ".tar.gz"

    common_top_part = (
    "RUN apt-get update\n"
    "RUN apt-get upgrade -y\n"
    "RUN apt-get autoremove -y\n"
    "RUN apt-get install -y build-essential\n"
    "RUN apt-get install -y gfortran\n" 
    "RUN apt-get install -y g++\n"
    "RUN apt-get install -y libncurses5-dev\n"
    "RUN apt-get install -y libreadline-dev\n"
    "RUN apt-get install -y flex\n"
    "RUN apt-get install -y bison\n"
    "RUN apt-get install -y libopenblas-dev\n"        
    "RUN apt-get install -y liblapacke-dev\n"
    "RUN apt-get install -y libcfitsio-dev\n"
    "RUN apt-get install -y wcslib-dev\n"
    "RUN apt-get install -y libhdf5-serial-dev\n" 
    "RUN apt-get install -y libfftw3-dev\n" 
    "RUN apt-get install -y libpython3-dev\n" 
    "RUN apt-get install -y libpython2.7-dev\n"
    "RUN apt-get install -y python-pip\n"           
    "RUN apt-get install -y python-numpy\n"
    "RUN apt-get install -y python-scipy\n"
    "RUN apt-get install -y libboost-python-dev\n" 
    "RUN apt-get install -y libboost-dev\n"         
    "RUN apt-get install -y libboost-filesystem-dev\n" 
    "RUN apt-get install -y libboost-program-options-dev\n" 
    "RUN apt-get install -y libboost-signals-dev\n"
    "RUN apt-get install -y libboost-system-dev\n"  
    "RUN apt-get install -y libboost-thread-dev\n"   
    "RUN apt-get install -y libboost-regex-dev\n"  
    "RUN apt-get install -y libcppunit-dev\n"   
    "RUN apt-get install -y git\n"
    "RUN apt-get install -y libffi-dev\n"     
    "RUN apt-get install -y libgsl-dev\n"        
    "RUN apt-get install -y liblog4cxx-dev\n"           
    "RUN apt-get install -y make\n"
    "RUN apt-get install -y patch\n"           
    "RUN apt-get install -y subversion\n"          
    "RUN apt-get install -y wget\n"          
    "RUN apt-get install -y docker\n"       
    "RUN apt-get install -y libxerces-c-dev\n"
    "RUN apt-get install -y libcurl4-openssl-dev\n"
    "# Make cmake from source\n"
    "RUN mkdir /usr/local/share/cmake\n"
    "WORKDIR /usr/local/share/cmake\n"
    "RUN wget https://github.com/Kitware/CMake/releases/download/v" + cmake_ver + "/" + cmake_source + "\n"
    "RUN tar -zxf " + cmake_source + "\n"
    "WORKDIR /usr/local/share/cmake/cmake-" + cmake_ver + "\n"
    "RUN ./bootstrap --system-curl\n"
    "RUN make\n"
    "RUN make install\n")

    common_bottom_part = (
    "RUN mkdir /usr/local/share/casacore\n"
    "RUN mkdir /usr/local/share/casacore/data\n"
    "WORKDIR /usr/local/share/casacore/data\n"
    "RUN wget ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar\n"
    "RUN mv WSRT_Measures.ztar WSRT_Measures.tar.gz\n"
    "RUN gunzip WSRT_Measures.tar.gz\n"
    "RUN tar -xf WSRT_Measures.tar\n"
    "RUN rm WSRT_Measures.tar\n"
    "RUN mkdir /var/lib/jenkins\n"
    "RUN mkdir /var/lib/jenkins/workspace\n"
    "WORKDIR /home\n"
    "RUN git clone " + git_repository + "\n"
    "WORKDIR /home/yandasoft\n"
    "RUN ./build_all.sh -C \"-DDATA_DIR=/usr/local/share/casacore/data\"\n" 
    "RUN ./build_all.sh -r\n")

    for machine_target in machine_targets:
        if machine_target == "generic":
            base_system_part = ("FROM ubuntu:bionic\n")

            for mpi_target in mpi_targets:
                (mpi_type, mpi_ver) = get_mpi_type_and_version(mpi_target)

                if (mpi_type == "mpich"):
                    if (mpi_ver == ""):
                        # if MPICH version is not specified, get the precompiled generic version
                        mpi_part = "RUN apt-get install -y mpich\n"

                    else:
                        # else (if version is specified), download the source from website and build           
                        mpich_dir = "https://www.mpich.org/static/downloads/" + mpi_ver

                        # TODO: Check whether the version is correct and the file exists

                        mpi_part = (
                        "WORKDIR /home\n"
                        "RUN wget " + mpich_dir + "/" + mpi_target + ".tar.gz\n"
                        "RUN gunzip " + mpi_target + ".tar.gz\n"
                        "RUN tar -xf " + mpi_target + ".tar\n"
                        "WORKDIR /home/" + mpi_target + "\n"
                        "RUN ./configure --prefix=\"/home/$USER/mpich-install\n"
                        "RUN make\n"
                        "RUN make install\n"
                        "ENV PATH=$PATH:/home/$USER/mpich-install/bin\n"
                        "ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/$USER/mpich-install/lib/:/usr/local/lib\n")

                elif (mpi_type == "openmpi"):
                    # Download the source from OpenMPI website and build

                    openmpi_common_top_part = (
                    "WORKDIR /home\n")

                    openmpi_common_bottom_part = (
                    "RUN ./configure\n"
                    "RUN make all install\n"
                    "ENV LD_LIBRARY_PATH=/usr/local/lib\n")

                    openmpi_ver = mpi_target[8:]

                    # TODO: Check whether the version number is correct

                    # Directory name for OpenMPI download
                    openmpi_dir = "https://download.open-mpi.org/release/open-mpi/v" + openmpi_ver[0:3]

                    # TODO: Check whether this file exist

                    openmpi_version_part = (
                    "RUN wget " + openmpi_dir + "/" + mpi_target + ".tar.gz\n"
                    "RUN gunzip " + mpi_target + ".tar.gz\n"
                    "RUN tar -xf " + mpi_target + ".tar\n"
                    "WORKDIR /home/" + mpi_target + "\n")

                    mpi_part = openmpi_common_top_part + openmpi_version_part + openmpi_common_bottom_part

                else:
                    print("ERROR: unknown MPI target: ", mpi_target)
                    quit()

                docker_target = DockerClass()
                docker_target.set_file_name("Dockerfile-casabase-" + mpi_target)
                docker_target.set_content(header + base_system_part + common_top_part + mpi_part + common_bottom_part)
                docker_target.set_image(image_prepend + mpi_target + image_append)
                docker_targets.append(docker_target)
            # Next mpi target

        elif (machine_target == "galaxy"):
            # Galaxy (of Pawsey) has Docker image with its MPICH implementation already baked into 
            # an Ubuntu base.
            base_system_part = ("FROM pawsey/mpi-base:latest\n")

            docker_target = DockerClass()
            docker_target.set_file_name("Dockerfile-casabase-" + machine_target)
            docker_target.set_content(header + base_system_part + common_top_part + common_bottom_part)
            docker_target.set_image(image_prepend + machine_target + image_append)
            docker_targets.append(docker_target)

        else:
            print("ERROR: unknown machine target: ", machine_target)
            quit()

    print("Docker target count:", len(docker_targets))

    print()

    if args.image:
        print("Making Docker images for all targets ...")
    else:
        print("This is a dry run. No Docker image will be created")

    for docker_target in docker_targets:
        docker_target.write()
        docker_command = ("docker build -t " + docker_target.image + " -f " + docker_target.file_name + " .")
        if args.image:
            subprocess.run(docker_command, shell=True)
        else:
            subprocess.run("echo " + docker_command, shell=True)


if (__name__ == "__main__"):
    if sys.version_info[0] == 3:
        main()
    else:
        raise ValueError("Must use Python 3")
