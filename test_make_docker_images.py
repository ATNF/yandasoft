# This is to test the script make_docker_images.py.
# Run this script using pytest.
# Author: Paulus Lahur <paulus.lahur@csiro.au>


import pytest
import os
import subprocess
import make_docker_images as m

def test_is_proper_name():
    assert m.is_proper_name("good_name") == True
    assert m.is_proper_name("bad name") == False

def test_set_recipe_name():
    d = m.DockerClass()
    d.set_recipe_name("good_name") 
    assert d.recipe_name == "good_name"
    with pytest.raises(ValueError):
        d.set_recipe_name("bad name") 
    
def test_set_recipe():
    d = m.DockerClass()
    d.set_recipe("recipe") 
    assert d.recipe == "recipe"
    with pytest.raises(ValueError):
        d.set_recipe("")
    with pytest.raises(TypeError):
        d.set_recipe(3) 

def test_set_image_name():
    d = m.DockerClass()
    d.set_image_name("good_name") 
    assert d.image_name == "good_name"
    with pytest.raises(ValueError):
        d.set_image_name("bad name") 

def test_write_recipe():
    d = m.DockerClass()
    with pytest.raises(ValueError):
        d.write_recipe()
    d.set_recipe_name("recipe_name")
    with pytest.raises(ValueError):
        d.write_recipe()
    d.set_recipe("recipe")
    d.write_recipe()
    with open("recipe_name", "r") as file:
        recipe = file.read()
        assert recipe == d.recipe
        # Clean up
        os.remove("recipe_name")

def test_get_build_command():
    d = m.DockerClass()
    with pytest.raises(ValueError):
        d.get_build_command()
    d.set_recipe_name("recipe_name")
    with pytest.raises(ValueError):
        d.get_build_command()
    d.set_image_name("image_name")
    assert d.get_build_command() == "docker build -t image_name -f recipe_name ."

def test_build_image():
    image_name = "docker_image_test"
    recipe_name = "Dockerfile_test"
    recipe = "FROM alpine:3.7\n"
    d = m.DockerClass()
    d.set_recipe_name(recipe_name)
    d.set_image_name(image_name)
    d.set_recipe(recipe)
    d.write_recipe()
    d.build_image()
    # Check whether the image is actually created
    assert subprocess.check_output(["docker", "images", "-q", image_name], text=True) != ""
    # Clean up
    subprocess.run(["docker", "rmi", image_name])
    # Note that the base image is not removed, just in case it's being used by others
    os.remove(recipe_name)

def test_get_mpi_type_and_version():
    with pytest.raises(TypeError):
        m.get_mpi_type_and_version(None)
    with pytest.raises(ValueError):
        m.get_mpi_type_and_version("s")
    with pytest.raises(ValueError):
        m.get_mpi_type_and_version("mpicz")
    with pytest.raises(ValueError):
        m.get_mpi_type_and_version("mpi_way_too_long_it_should_trigger_error")
    assert m.get_mpi_type_and_version("mpich") == ("mpich", "")
    assert m.get_mpi_type_and_version("mpich-3.3.2") == ("mpich", "3.3.2")
    assert m.get_mpi_type_and_version("openmpi-4.0.2") == ("openmpi", "4.0.2")


def test_make_base_image():
    d = m.make_base_image("generic", "mpich", "testbase", ":latest", False)

    ref = (
    "# This file is automatically created by ./make_docker_images.py\n"
    "FROM ubuntu:bionic\n"
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
    "RUN wget https://github.com/Kitware/CMake/releases/download/v3.15.7/cmake-3.15.7.tar.gz\n"
    "RUN tar -zxf cmake-3.15.7.tar.gz\n"
    "WORKDIR /usr/local/share/cmake/cmake-3.15.7\n"
    "RUN ./bootstrap --system-curl\n"
    "RUN make\n"
    "RUN make install\n"
    "RUN apt-get install -y mpich\n"
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
    "RUN git clone https://github.com/ATNF/yandasoft.git\n"
    "WORKDIR /home/yandasoft\n"
    "RUN ./build_all.sh -C \"-DDATA_DIR=/usr/local/share/casacore/data\"\n"
    "RUN ./build_all.sh -r\n")

    recipe_lines = d.recipe.splitlines()
    ref_lines = ref.splitlines()
    n = len(ref_lines)
    assert len(recipe_lines) == n
    # Check every single line except the first line
    for i in range(1,n):
        assert recipe_lines[i] == ref_lines[i]
    # Clean up
    os.remove(d.recipe_name)


def test_make_final_image():
    d = m.make_final_image("generic", "mpich", "testfinal", ":latest", "base_image", False)

    ref = (
    "# This file is automatically created by ./make_docker_images.py\n"
    "FROM base_image\n"
    "RUN apt-get install -y mpich\n"
    "WORKDIR /home/yandasoft\n"
    "RUN git pull https://github.com/ATNF/yandasoft.git\n"
    "RUN ./build_all.sh -a -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -y -O \"-DHAVE_MPI=1\"\n"
    "RUN ./build_all.sh -e -O \"-DHAVE_MPI=1\"\n")

    recipe_lines = d.recipe.splitlines()
    ref_lines = ref.splitlines()
    n = len(ref_lines)
    assert len(recipe_lines) == n
    # Check every single line except the first line
    for i in range(1,n):
        assert recipe_lines[i] == ref_lines[i]
    # Clean up
    os.remove(d.recipe_name)
