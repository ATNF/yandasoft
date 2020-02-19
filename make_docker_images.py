#!/usr/bin/env python3
#
# Create Docker images for all cases as specified in the settings below.
# Possible targets are:
# - Specific machine (currently only Galaxy, with Cray MPICH)
# - Generic machine, with these MPI implementations
#   - MPICH
#   - OpenMPI of various versions
#
# Note that the images are split into base image (for components that are 
# seldom changed) and the final image (Yandasoft).
# This is done for deployment efficiency.
#
# Usage:
# 1) Make sure that the targets are correct (see SETTINGS section below).
# 2) Execute the script.
#    If you want only the final images (since base image is still the same):
#    ./make_docker_image.py -f
#
#    If you want to make both base and final images:
#    ./make_docker_image.py -bf
#
#    If you want no image, just Dockerfiles, for example for dry run: 
#    ./make_docker_image.py
#
# Author: Paulus Lahur
#
#------------------------------------------------------------------------------
# SETTINGS
#
# Set machine targets in the list below.
# Currently, we target generic HPCs and Galaxy.
machine_targets = ["generic", "galaxy"]

# Set MPI implementations for "generic" machine in the list below.
# Note that a specific machine already has its MPI specified.
# A valid mpi target is either "mpich" or "openmpi-X.Y.Z", 
# where X, Y and Z are version numbers (major, minor and revision).
# Our current MPI targets:
mpi_targets = ["mpich", "openmpi-4.0.2", "openmpi-3.1.4", "openmpi-2.1.6"]

#------------------------------------------------------------------------------
# CODE
# TODO: Add logging
# TODO: Add timing
# TODO: Add error handling, as this is going to be used within CI/CD

import sys
import argparse
import subprocess

# Git repository of Yandasoft
git_repository = "https://github.com/ATNF/yandasoft.git"

# Header for all automatically generated Dockerfiles
header = ("# This file is automatically created by " + __file__ + "\n")

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


def make_base_image(machine, mpi, prepend, append, actual):
    '''
    Make base image for components that are seldom changed:
    base OS, upgrades, standard libraries and apps, Casacore and Casarest.
    '''
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
    "RUN tar -zxf WSRT_Measures.tar.gz\n"
    "RUN rm WSRT_Measures.tar.gz\n"
    "RUN mkdir /var/lib/jenkins\n"
    "RUN mkdir /var/lib/jenkins/workspace\n"
    "WORKDIR /home\n"
    "RUN git clone " + git_repository + "\n"
    "WORKDIR /home/yandasoft\n"
    "RUN ./build_all.sh -C \"-DDATA_DIR=/usr/local/share/casacore/data\"\n" 
    "RUN ./build_all.sh -r\n")

    if machine == "generic":
        base_system_part = ("FROM ubuntu:bionic\n")

        (mpi_type, mpi_ver) = get_mpi_type_and_version(mpi)

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
                "RUN wget " + mpich_dir + "/" + mpi + ".tar.gz\n"
                "RUN tar -zxf " + mpi + ".tar.gz\n"
                "WORKDIR /home/" + mpi + "\n"
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

            openmpi_ver = mpi[8:]

            # TODO: Check whether the version number is correct

            # Directory name for OpenMPI download
            openmpi_dir = "https://download.open-mpi.org/release/open-mpi/v" + openmpi_ver[0:3]

            # TODO: Check whether this file exist

            openmpi_version_part = (
            "RUN wget " + openmpi_dir + "/" + mpi + ".tar.gz\n"
            "RUN tar -zxf " + mpi + ".tar.gz\n"
            "WORKDIR /home/" + mpi + "\n")

            mpi_part = openmpi_common_top_part + openmpi_version_part + openmpi_common_bottom_part

        else:
            print("ERROR: unknown MPI target: ", mpi)
            return None

        docker_target = DockerClass()
        docker_target.set_file_name("Dockerfile-casabase-" + mpi)
        docker_target.set_content(header + base_system_part + common_top_part + mpi_part + common_bottom_part)
        docker_target.set_image(prepend + mpi + append)

    elif (machine == "galaxy"):
        # Galaxy (of Pawsey) has Docker image with its MPICH implementation already baked into 
        # an Ubuntu base.
        base_system_part = ("FROM pawsey/mpi-base:latest\n")

        docker_target = DockerClass()
        docker_target.set_file_name("Dockerfile-casabase-" + machine)
        docker_target.set_content(header + base_system_part + common_top_part + common_bottom_part)
        docker_target.set_image(prepend + machine + append)

    else:
        print("ERROR: unknown machine target: ", machine)
        return None

    docker_target.write()
    docker_command = ("docker build -t " + docker_target.image + " -f " + docker_target.file_name + " .")
    if actual:
        subprocess.run(docker_command, shell=True)
    else:
        subprocess.run("echo " + docker_command, shell=True)

    return docker_target



def make_final_image(machine, mpi, prepend, append, base_image, actual):
    '''
    Make the final image on top of base image.
    '''

    common_bottom_part = (
    "WORKDIR /home/yandasoft\n"
    "RUN git pull " + git_repository + "\n"
    "RUN ./build_all.sh -a -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -y -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -e -O \"-DHAVE_MPI=1\"\n")

    base_part = ("FROM " + base_image + "\n")

    if machine == "generic":
        (mpi_type, mpi_ver) = get_mpi_type_and_version(mpi)

        if (mpi_type == "mpich"):
            if (mpi_ver == ""):
                # if MPICH version is not specified, get the precompiled generic version
                #base_part = ("FROM csirocass/casabase-mpich:latest\n")
                mpi_part = "RUN apt-get install -y mpich\n"

            else:
                # Otherwise, specific version of MPICH
                #base_part = ("FROM csirocass/casabase-mpich-" + mpi_ver + ":latest\n")
                
                # Download the source from MPICH website and build from source     
                mpich_dir = "https://www.mpich.org/static/downloads/" + mpi_ver

                # TODO: Check whether the version is correct and the file exists

                mpi_part = (
                "WORKDIR /home\n"
                "RUN wget " + mpich_dir + "/" + mpi + ".tar.gz\n"
                "RUN tar -zxf " + mpi + ".tar.gz\n"
                "WORKDIR /home/" + mpi + "\n"
                "RUN ./configure --prefix=\"/home/$USER/mpich-install\n"
                "RUN make\n"
                "RUN make install\n"
                "ENV PATH=$PATH:/home/$USER/mpich-install/bin\n"
                "ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/$USER/mpich-install/lib/:/usr/local/lib\n")

        elif (mpi_type == "openmpi"):
            #base_part = ("FROM csirocass/casabase-openmpi-" + mpi_ver + ":latest\n")

            # Download the source from OpenMPI website and build from source
            openmpi_common_top_part = (
            "WORKDIR /home\n")

            openmpi_common_bottom_part = (
            "RUN ./configure\n"
            "RUN make all install\n"
            "ENV LD_LIBRARY_PATH=/usr/local/lib\n")

            openmpi_ver = mpi[8:]

            # TODO: Check whether the version number is correct

            # Directory name for OpenMPI download
            openmpi_dir = "https://download.open-mpi.org/release/open-mpi/v" + openmpi_ver[0:3]

            # TODO: Check whether this file exist

            openmpi_version_part = (
            "RUN wget " + openmpi_dir + "/" + mpi + ".tar.gz\n"
            "RUN tar -zxf " + mpi + ".tar.gz\n"
            "WORKDIR /home/" + mpi + "\n")

            mpi_part = openmpi_common_top_part + openmpi_version_part + openmpi_common_bottom_part

        else:
            print("ERROR: unknown MPI target: ", mpi)
            return None

        docker_target = DockerClass()
        docker_target.set_file_name("Dockerfile-yandasoft-" + mpi)
        docker_target.set_content(header + base_part + mpi_part + common_bottom_part)
        docker_target.set_image(prepend + mpi + append)

    elif (machine == "galaxy"):
        #base_part = ("FROM csirocass/casabase-galaxy:latest\n")

        docker_target = DockerClass()
        docker_target.set_file_name("Dockerfile-yandasoft" + machine)
        docker_target.set_content(header + base_part + common_bottom_part)
        docker_target.set_image(prepend + machine + append)

    else:
        print("ERROR: unknown machine target: ", machine)
        return None

    docker_target.write()
    docker_command = ("docker build -t " + docker_target.image + " -f " + docker_target.file_name + " .")
    if actual:
        subprocess.run(docker_command, shell=True)
    else:
        subprocess.run("echo " + docker_command, shell=True)

    return docker_target



def make_batch_file(machine, mpi):
    '''
    Make sample batch files for SLURN
    '''
    # Currently not being used

    batch_common_part = (
    "#!/bin/bash -l\n"
    "## This file is automatically created by " + __file__ + "\n"
    "#SBATCH --ntasks=5\n"
    "##SBATCH --ntasks=305\n"
    "#SBATCH --time=02:00:00\n"
    "#SBATCH --job-name=cimager\n"
    "#SBATCH --export=NONE\n\n"
    "module load singularity/3.5.0\n")

    mpi_type = get_mpi_type(mpi)
    if (mpi_type == "mpich"):
        module = "mpich/3.3.0"
        image = "yandasoft-mpich_latest.sif"
        batch_mpi_part = (
        "module load " + module + "\n\n"
        "mpirun -n 5 singularity exec " + image +
        " cimager -c dirty.in > dirty_${SLURM_JOB_ID}.log\n")

    elif (mpi_type == "openmpi"):
        openmpi_ver = get_openmpi_version(mpi)
        if (openmpi_ver != None):
            module = "openmpi/" + openmpi_ver + "-ofed45-gcc"
            image = "yandasoft-" + openmpi_ver + "_latest.sif"
            batch_mpi_part = (
            "module load " + module + "\n\n"
            "mpirun -n 5 -oversubscribe singularity exec " + image +
            " cimager -c dirty.in > dirty_${SLURM_JOB_ID}.log\n")

    batch_file = "sample-" + machine + "-" + mpi + ".sbatch"
    print("Making batch file:", batch_file)
    f = open(batch_file, "w")
    f.write(batch_common_part + batch_mpi_part)
    f.close()


def main():
    '''
    The main code
    '''
    parser = argparse.ArgumentParser(
        description="Make Docker images for various MPI implementations",
        epilog="The targets can be changed from inside the script (the SETTINGS section)")
    parser.add_argument('-b', '--base_image', help='Create base image', action='store_true')
    parser.add_argument('-f', '--final_image', help='Create final image', action='store_true')
    #parser.add_argument('-s', '--slurm', help='Create sample batch files for SLURM', action='store_true')
    args = parser.parse_args()

    base_prepend = "csirocass/casabase-"
    base_append = ":latest"
    final_prepend = "csirocass/yandasoft-"
    final_append = ":latest"

    print("Making Dockerfiles for all targets ...")
    if args.base_image:
        print("Making base images ...")
    else:
        print("Base image will not be made")

    if args.final_image:
        print("Making final images ...")
    else:
        print("Final image will not be made")

    for machine in machine_targets:
        if machine == "generic":
            for mpi in mpi_targets:
                docker = make_base_image(machine, mpi, base_prepend, base_append, args.base_image)
                if docker != None:
                    docker = make_final_image(machine, mpi, final_prepend, final_append, docker.image, 
                        args.final_image)
        else:
            # Specific machine
            docker = make_base_image(machine, None, base_prepend, base_append, args.base_image)
            if docker != None:
                docker = make_final_image(machine, None, final_prepend, final_append, docker.image, 
                    args.final_image)


if (__name__ == "__main__"):
    if sys.version_info[0] == 3:
        main()
    else:
        raise ValueError("Must use Python 3")
