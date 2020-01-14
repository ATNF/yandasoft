#!/usr/bin/env python3
#
# This script automates 3 things:
# 1) Create Dockerfiles for all MPI implementations.
# 2) If required, create the corresponding Docker images.
# 3) Create sample SLURM batch files.
# 
# Author: Paulus Lahur
#
#------------------------------------------------------------------------------
# SETTINGS
#
# Set target MPI implementations in the list below.
# How to write the target:
# - MPICH must be written as: "mpich".
# - OpenMPI must be in this format: "openmpi-X.Y.Z", where X, Y and Z are
#   version's major, minor and revision number, respectively.
# Note that the targets are set internally rather than externally 
# (ie. as arguments) for these reasons:
# - Minimise error checking
# - The target will be passed on as argument in shell commands.
#   Bad target might be mistaken as another command.

mpi_targets = ["MPICH", "openmpi-2.1.6", "openmpi-3.1.4", "openmpi-4.0.2"]

# Docker image name is in this format: 
# target_prepend + mpi_target + target_append

target_prepend = "lahur/yandasoft-"
target_append = ":latest"

# Set True if this is just a dry run. No Docker image will be created.

dry_run = True

# HPC name that will be used as the prepend to the batch file name.

HPC = "pearcey"

#------------------------------------------------------------------------------
# CODE

# Sanitizing parameters

mpi_targets = list(map(str.lower, mpi_targets))


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
    import subprocess

    print("There are", len(mpi_targets), "MPI targets")

    print("Making Dockerfiles for all targets ...")

    dockerfiles = []

    common_top_part = (
    "# This file is automatically created by " + __file__ + "\n"
    "FROM ubuntu:bionic\n"
    "RUN apt-get update\n"
    "RUN apt-get upgrade -y\n"
    "RUN apt-get autoremove -y\n"
    "RUN apt-get install -y cmake\n"          
    "RUN apt-get install -y flex bison\n"
    "RUN apt-get install -y gfortran\n"           
    "RUN apt-get install -y git\n"            
    "RUN apt-get install -y g++\n"                
    "RUN apt-get install -y libboost-dev\n"         
    "RUN apt-get install -y libboost-python-dev\n" 
    "RUN apt-get install -y libboost-filesystem-dev\n" 
    "RUN apt-get install -y libboost-program-options-dev\n" 
    "RUN apt-get install -y libboost-signals-dev\n"
    "RUN apt-get install -y libboost-system-dev\n"  
    "RUN apt-get install -y libboost-thread-dev\n"   
    "RUN apt-get install -y libboost-regex-dev\n"  
    "RUN apt-get install -y libcppunit-dev\n"   
    "RUN apt-get install -y libcfitsio-dev\n"      
    "RUN apt-get install -y libffi-dev\n"     
    "RUN apt-get install -y libfftw3-dev\n"         
    "RUN apt-get install -y libgsl-dev\n"        
    "RUN apt-get install -y liblog4cxx-dev\n"           
    "RUN apt-get install -y libopenblas-dev\n"        
    "RUN apt-get install -y libpython-dev\n"
    "RUN apt-get install -y make\n"
    "RUN apt-get install -y patch\n"           
    "RUN apt-get install -y python-pip\n"           
    "RUN apt-get install -y subversion\n"          
    "RUN apt-get install -y wget\n"          
    "RUN apt-get install -y docker\n"       
    "RUN apt-get install -y python-numpy\n"   
    "RUN apt-get install -y python-scipy\n"
    "RUN apt-get install -y wcslib-dev\n"
    "RUN apt-get install -y libxerces-c-dev\n")            

    common_bottom_part = (
    "RUN mkdir /usr/local/share/casacore\n"
    "RUN mkdir /usr/local/share/casacore/data\n"
    "WORKDIR /usr/local/share/casacore/data\n"
    "RUN wget ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar\n"
    "RUN mv WSRT_Measures.ztar WSRT_Measures.tar.gz\n"
    "RUN gunzip WSRT_Measures.tar.gz\n"
    "RUN tar -xvf WSRT_Measures.tar\n"
    "RUN rm WSRT_Measures.tar\n"
    "RUN mkdir /var/lib/jenkins\n"
    "RUN mkdir /var/lib/jenkins/workspace\n"
    "WORKDIR /home\n"
    "RUN git clone https://ord006@bitbucket.csiro.au/scm/askapsdp/yandasoft.git\n"
    "WORKDIR /home/yandasoft\n"
    "RUN ./build_all.sh -C \"-DDATA_DIR=/usr/local/share/casacore/data\"\n" 
    "RUN ./build_all.sh -r\n"
    "RUN ./build_all.sh -a -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -y -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -e -O \"-DHAVE_MPI=1\"\n")

    for mpi_target in mpi_targets:
        dockerfile = "Dockerfile-" + mpi_target
        dockerfiles.append(dockerfile)
        print("Making Dockerfile:", dockerfile)

        if (mpi_target == "mpich"):
            mpi_part = "RUN apt-get install -y mpich\n"

        elif (mpi_target[0:8] == "openmpi-"):
            # Note that OpenMPI is more complicated than MPICH, because:
            # - There are multiple versions
            # - OpenMPI must be built from the source code
            # - The source code must be downloaded from OpenMPI website first
            # - The version dictates the directory it is downloaded from

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
            "RUN tar -xvf " + mpi_target + ".tar\n"
            "WORKDIR /home/" + mpi_target + "\n")

            mpi_part = openmpi_common_top_part + openmpi_version_part + openmpi_common_bottom_part

        else:
            print("ERROR: unknown MPI target: ", mpi_target)
            quit()


        f = open(dockerfile, "w")
        f.write(common_top_part + mpi_part + common_bottom_part)
        f.close()

    print()

    if dry_run:
        print("This is a dry run. No Docker image will be created")
    else:
        print("Making Docker images for all targets ...")

    for mpi_target in mpi_targets:
        dockerfile = "Dockerfile-" + mpi_target
        target = target_prepend + mpi_target + target_append
        docker_command = "docker build -t " + target + " -f " + dockerfile + " ."
        if dry_run:
            subprocess.run("echo " + docker_command, shell=True)
        else:
            subprocess.run(docker_command, shell=True)

            # TODO: Deal with possible error

    # Consider: Add automatic upload to DockerHub? This requires Docker login.

    print()
    print("Making sample batch files ...")

    for mpi_target in mpi_targets:
        batch_common_part = (
        "#!/bin/bash -l\n"
        "## This file is automatically created by " + __file__ + "\n"
        "#SBATCH --ntasks=5\n"
        "##SBATCH --ntasks=305\n"
        "#SBATCH --time=02:00:00\n"
        "#SBATCH --job-name=cimager\n"
        "#SBATCH --export=NONE\n\n"
        "module load singularity/3.5.0\n")

        mpi_type = get_mpi_type(mpi_target)
        if (mpi_type == "mpich"):
            module = "mpich/3.3.0"
            image = "yandasoft-mpich_latest.sif"
            batch_mpi_part = (
            "module load " + module + "\n\n"
            "mpirun -n 5 singularity exec " + image +
            " cimager -c dirty.in > dirty_${SLURM_JOB_ID}.log\n")

        elif (mpi_type == "openmpi"):
            openmpi_ver = get_openmpi_version(mpi_target)
            if (openmpi_ver != None):
                module = "openmpi/" + openmpi_ver + "-ofed45-gcc"
                image = "yandasoft-" + openmpi_ver + "_latest.sif"
                batch_mpi_part = (
                "module load " + module + "\n\n"
                "mpirun -n 5 -oversubscribe singularity exec " + image +
                " cimager -c dirty.in > dirty_${SLURM_JOB_ID}.log\n")
            else:
                break
        else:
            break

        batch_file = "sample-" + HPC + "-" + mpi_target + ".sbatch"
        print("Making batch file:", batch_file)
        f = open(batch_file, "w")
        f.write(batch_common_part + batch_mpi_part)
        f.close()


if (__name__ == "__main__"):
    main()