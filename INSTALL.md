# How to Install the package

Yandasoft is supplied with a bash script `build_all.sh` - the purpose of this is to pull and build all the non-system dependencies. It can also build the system dependencies for you - but this is not recommended. 

## System dependencies <BETA>

### Ubuntu

* cmake                         `# many` \
* flex bison                    `# casacore` \
* gfortran                      `# many` \
* git                           `# many` \
* g++                           `# many` \
* libboost-dev                  `# casacore` \
* libboost-filesystem-dev       `# yandasoft` \
* libboost-program-options-dev  `# base-askap` \
* libboost-signals-dev          `# base-askap` \
* libboost-system-dev           `# casarest` \
* libboost-thread-dev           `# casarest` \
* libcfitsio-dev                `# casacore` \
* libfftw3-dev                  `# casacore` \
* libgsl-dev                    `# many` \
* liblog4cxx-dev                `# yandasoft` \
* libopenblas-dev               `# casacore` \
* libopenmpi-dev                `# adios, casacore, oskar` \
* libpython-dev                 `# casacore, few python packages` \
* make                          `# many` \
* patch                         `# lofar-common` \
* python-pip                    `# so we can pip install virtualenv` \
* subversion                    `# lofar-blob, lofar-common` \
* wcslib-dev                    `# casacore`

### OSX

TBD

## Running `build_all.sh`

`build_all.sh -h`


```
Usage: ./build_all.sh [options]

Options:
 -s <system>   Target system, supported values are osx (default), centos, ubuntu
 -x <compiler> Compiler suite to use, supported values are gcc (default) and clang
 -m <cmake>    Cmake binary to use, must be >= 3.1, defaults to cmake
 -j <jobs>     Number of parallel compilation jobs, defaults to 1
 -p <prefix>   Prefix for installation, defaults to /usr/local
 -w <workdir>  Working directory, defaults to .
 -S 	     Install system dependencies. 
 -W            Remove the working directory at the end of the build
 -C <opts      Install Casacore + cmake options
 -A <opts>     Install ASKAP dependencies + cmake options
 -R <opts>     Install casarest + cmake options
 -Y <opts>     Install YandaSoft + cmake options
 -U 	     clean and uninstall yandasoft and dependencies (except casacore/casarest)
 -P            Use Python 3 
```

Typically you would install the system dependencies, then casacore/casarest, then the askap dependencies, then yandasoft. This can be done by simply:

`build_all.sh -c -a -y`

Assuming you want everything in ```/usr/local```. You can specify a prefix. If you find you need to provide cmake options the equivelant is:

`build_all.sh -C <opt> -A <opt> -Y <opt>`

If you want to clean up the mess:

`build_all.sh -U`

Will clean and uninstall everything. It will leave any downloaded repositories in place though.
