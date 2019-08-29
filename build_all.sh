#!/bin/bash
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
#
# (c) Copyright 2019 CSIRO
# Australia Telescope National Facility (ATNF)
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# PO Box 76, Epping NSW 1710, Australia
# atnf-enquiries@csiro.au
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#
# A build script for yandasoft heavily based - some would say copied from the 
# build script used for the full JACAL build on summit.
#
# I have simplified it for use as a continuous integration script for use in 
# organising the build of yandasoft.
#
# authors:
# Original author Rodrigo Tobar
# Changes to concentrate on the yandasoft build: Stephen Ord
#

print_usage() {
	echo "Usage: $0 [options]"
	echo
	echo "Options:"
	echo " -s <system>     Target system, supported values are osx (default), centos, ubuntu"
	echo " -x <compiler>   Compiler suite to use, supported values are gcc (default) and clang"
	echo " -m <cmake>      Cmake binary to use, must be >= 3.1, defaults to cmake"
	echo " -j <jobs>       Number of parallel compilation jobs, defaults to 1"
	echo " -p <prefix>     Prefix for installation, defaults to /usr/local"
	echo " -w <workdir>    Working directory, defaults to ."
	echo " -S 	       Install system dependencies. "
	echo " -W              Remove the working directory at the end of the build"
	echo " -C <opts> | -c  Install Casacore + cmake options"
	echo " -A <opts> | -a  Install ASKAP dependencies + cmake options"
	echo " -R <opts> | -r  Install casarest + cmake options"
	echo " -Y <opts> | -y  Install YandaSoft + cmake options"
	echo " -E <opts> | -e  Install Extra (analysis + pipelinetasks) + cmake options"
	echo " -U 	       clean and uninstall yandasoft and dependencies (except casacore/casarest)"
        echo " -O <opts>       Options to apply to all builds."
        echo " -P              Use Python 3 "
}

try() {
	"$@"
	status=$?
	if [ $status -ne 0 ]; then
		echo "Command exited with status $status, aborting build now: $@" 1>&2
		exit 1;
	fi
}

check_supported_values() {
	val_name=$1
	given_val=$2
	shift; shift
	for supported in "$@"; do
		if [ "$given_val" == "$supported" ]; then
			return
		fi
	done
	echo "Unsupported $val_name: $given_val" 1>&2
	echo "Supported $val_name values are: $@" 1>&2
	exit 1
}

system=osx
compiler=gcc
cmake=cmake
jobs=1
prefix=/usr/local
workdir="$PWD"
remove_workdir=no
build_oskar=yes
use_python3=no

install_system_dependencies=no
install_casacore=no
install_casarest=no
install_askap_dependencies=no
install_yandasoft=no
install_extra=no
clean_askap_dependencies=no
clean_yandasoft=no
build_adios=no
casacore_version=master
casacore_opts=
casarest_version=components-only
casarest_opts=
askap_opts=
yandasoft_opts=
extra_opts=
opts=

if [ $# -eq 0 ]; then 
	print_usage
	exit 0
fi

while getopts "A:ah?s:cm:j:p:w:WPoiC:cR:rY:yO:USx:eE:" opt
do
	case "$opt" in
		[h?])
			print_usage
			exit 0
			;;
		s)
			system="$OPTARG"
			;;
		x)
			compiler="$OPTARG"
			;;
		m)
			cmake="$OPTARG"
			;;
		j)
			jobs="$OPTARG"
			;;
		p)
			prefix="$OPTARG"
			;;
		w)
			workdir="$OPTARG"
			;;
		W)
			remove_workdir=yes
			;;
		o)
			build_oskar=no
			;;
		S)
			install_system_dependencies=yes
			;;
		A)
			install_askap_dependencies=yes
			askap_opts="$OPTARG"
			;;
		a)		
			install_askap_dependencies=yes
			;;
		C)
			casacore_opts="$OPTARG"
			install_casacore=yes
			;;
		c)	
			install_casacore=yes
			;;
		e)	install_extra=yes
			;;
		E)	extra_opts="$OPTARG"
			install_extra=yes
			;;
		R)
			casarest_opts="$OPTARG"
			install_casarest=yes
			;;
		r)	
			install_casarest=yes
			;;
		Y)
			yandasoft_opts="$OPTARG"
			install_yandasoft=yes
			;;
		y)
				
			install_yandasoft=yes
			;;
		O)
			opts="$OPTARG"
			;;
		P)
			use_python3=yes
			;;
		U)	
			install_system_dependencies=no
			install_askap_dependencies=no
			install_casacore=no
			install_yandasoft=no
			clean_askap_dependencies=yes
			clean_yandasoft=yes
			;;

		*)
			print_usage 1>&2
			exit 1
			;;
	esac
done
yandasoft_opts="${yandasoft_opts} ${opts} ${askap_opts} ${extra_opts} -DCASACORE_ROOT_DIR=${prefix}"
check_supported_values system $system centos ubuntu osx
check_supported_values compiler $compiler gcc clang cray
check_supported_values casacore_version $casacore_version master 2.4.0 2.0.3
if [ $casacore_version != master ]; then
	build_adios=no
fi

if [ $EUID == 0 ]; then
	SUDO=
else
	SUDO=sudo
fi

if [ $system == centos ]; then
	cmake=cmake3
fi

install_s_dependencies() {
	if [ $system == ubuntu ]; then
		$SUDO apt update
		$SUDO apt install -y \
		    cmake                         `# many` \
		    flex bison                    `# casacore` \
		    gfortran                      `# many` \
		    git                           `# many` \
		    g++                           `# many` \
		    libboost-dev                  `# casacore` \
		    libboost-filesystem-dev       `# yandasoft` \
		    libboost-program-options-dev  `# base-askap` \
		    libboost-signals-dev          `# base-askap` \
		    libboost-system-dev           `# casarest` \
		    libboost-thread-dev           `# casarest` \
		    libcfitsio-dev                `# casacore` \
		    libffi-dev                    `# cryptography (python) -> paramiko -> daliuge` \
		    libfftw3-dev                  `# casacore` \
		    libgsl-dev                    `# many` \
		    liblog4cxx-dev                `# yandasoft` \
		    libopenblas-dev               `# casacore` \
		    libopenmpi-dev                `# adios, casacore, oskar` \
		    libpython-dev                 `# casacore, few python packages` \
		    make                          `# many` \
		    patch                         `# lofar-common` \
		    python-pip                    `# so we can pip install virtualenv` \
		    subversion                    `# lofar-blob, lofar-common` \
		    wcslib-dev                    `# casacore`
		$SUDO apt clean
		if [ $compiler == clang ]; then
			$SUDO apt install -y clang
		fi
	elif [ $system == centos ]; then
		$SUDO yum --assumeyes install \
		    boost-devel    `# casacore` \
		    cfitsio-devel  `# casacore` \
		    cmake3         `# many` \
		    fftw3-devel    `# casacore` \
		    flex bison     `# casacore, lofar-common` \
		    gcc-c++        `# many, including clang itself` \
		    git            `# many` \
		    gsl-devel      `# casacore, yandasoft` \
		    libffi-devel   `# cryptography (python) -> paramiko -> daliuge` \
		    log4cxx-devel  `# yandasoft` \
		    make           `# many` \
		    openblas-devel `# casacore` \
		    openmpi-devel  `# adios, casacore, oskar, yandasoft` \
		    openssl-devel  `# cryptography (see above)` \
		    patch          `# lofar-common` \
		    python-devel   `# casacore, oskar, few python packages` \
		    python-pip     `# so we can pip install virtualenv` \
		    svn            `# lofar-blob, lofar-common` \
		    wcslib-devel   `# casacore`
		if [ $compiler == clang ]; then
			$SUDO yum --assumeyes install clang
		fi
		if [ $use_python3 == yes ]; then
			$SUDO yum --assumeyes install python34 python34-devel python34-pip boost-python34-devel
		fi
		$SUDO yum clean all
	elif [ $system == osx ]; then
		echo "OSX system dependencies not supported"
	fi
}

repo2dir() {
	d=`basename $1`
	echo ${d%%.git}
}

# Nice-to-use macro
clean_and_uninstall() {
	sourcedir="$PWD"
	if [ -d `repo2dir $1` ]; then
		cd `repo2dir $1`
		if [ -d build ]; then
		   cd build
		   if [ -f install_manifest.txt ]; then
		   		try make clean
				try make uninstall
                   fi
		fi
	fi
	cd "$sourcedir"
}

build_and_install() {
	sourcedir="$PWD"
	ref=$2
	is_branch=yes
	is_merge=no
	if [[ "$ref" =~ COMMIT-(.*) ]]; then
		ref=${BASH_REMATCH[1]}
		is_branch=no
	elif [[ "$ref" =~ MERGE-(.*) ]]; then
		ref=${BASH_REMATCH[1]}
		is_branch=no
		is_merge=yes
	fi
	if [ ! -d `repo2dir $1` ]; then
		try git clone $1
		cd `repo2dir $1`
		if [ $is_branch == yes ]; then
			git checkout -b $ref origin/$ref
		elif [ $is_merge == yes ]; then
			git config user.email "you@example.com"
			git config user.name "Your Name"
			git merge --no-edit remotes/origin/$ref
		else
			git checkout -b working_copy
			git reset --hard $ref
		fi
	else
		cd `repo2dir $1`
		git pull
	fi
	shift; shift
	test -d build || try mkdir build
	cd build
	if [ $compiler == clang ]; then
		comp_opts="-DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
	elif [ $compiler == cray ]; then
		comp_opts="-DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc"
	else
		comp_opts="-DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
	fi

        cmakeline=" ${cmake} -DCMAKE_INSTALL_PREFIX="$prefix" $comp_opts "$@" .."
        echo "$cmakeline"
	try ${cmake} -DCMAKE_INSTALL_PREFIX="$prefix" $comp_opts "$@" ..
	try make all -j${jobs}
	try make install -j${jobs}
	cd "$sourcedir"
}

source_venv() {
	if [ ! -f $prefix/bin/activate ]; then
		if [ $use_python3 == yes ]; then
			python3 -m venv $prefix
		else
			if [ -n "`command -v virtualenv 2> /dev/null`" ]; then
				try virtualenv $prefix
			else
				try pip install --user virtualenv
				try ~/.local/bin/virtualenv $prefix
			fi
		fi
	fi
	source $prefix/bin/activate
}

original_dir="$PWD"
if [ ! -d "$workdir" ]; then
	try mkdir -p "$workdir"
fi
cd "$workdir"

if [ $install_system_dependencies == yes ]; then
	install_s_dependencies
fi

# CentOS, you cheecky...
if [ $system == centos ]; then
	source /etc/profile.d/modules.sh
	module load mpi
fi

# Setup our environment
export LD_LIBRARY_PATH=$prefix/:$LD_LIBRARY_PATH



#casacore and casarest

if [ $casacore_version == master -a $build_adios == yes ]; then
	casacore_opts+=" -DUSE_ADIOS2=ON"
	casacore_version=summit_demo
fi
if [ $casacore_version != master ]; then
	casacore_version=COMMIT-v$casacore_version
fi
if [ $install_casacore == yes ]; then
	build_and_install https://github.com/casacore/casacore $casacore_version -DBUILD_TESTING=OFF $casacore_opts
	if [ $casacore_version == summit_demo ]; then
		# Lets reset this back to master to save handling lots os special cases below!
		casacore_version=master
	fi

	if [ $casacore_version == master ]; then
		casarest_version=components-only
	elif [ $casacore_version == COMMIT-v2.4.0 ]; then
		casarest_version=COMMIT-467ed6d
	else
		casarest_version=COMMIT-v1.4.1
	fi
fi
if [ $install_casarest == yes ]; then
	build_and_install https://github.com/steve-ord/casarest $casarest_version -DBUILD_TESTING=OFF $casarest_opts
fi

if [ $clean_askap_dependencies == yes ]; then

	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/lofar-common.git
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/lofar-blob.git
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-askap.git 
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-logfilters.git 
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-imagemath.git 
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-askapparallel.git 
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-scimath.git 
	clean_and_uninstall https://bitbucket.csiro.au/scm/askapsdp/base-accessors.git 
fi

if [ $install_askap_dependencies == yes ]; then
	
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/lofar-common.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/lofar-blob.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-askap.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-logfilters.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-imagemath.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-askapparallel.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-scimath.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-accessors.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/base-components.git master $yandasoft_opts
fi
if [ $install_extra == yes ]; then
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/askap-pipelinetasks.git master $yandasoft_opts
	build_and_install https://bitbucket.csiro.au/scm/askapsdp/askap-analysis.git master $yandasoft_opts
fi
	


if [ $clean_yandasoft == yes ]; then
   startdir="$PWD"
   if [ -d build ]; then
	echo "yandasoft build directory already exists"
	cd build
	if [ -f install_manifest.txt ]; then
		make clean
		make uninstall
		
	fi
   fi
   cd "$startdir"
fi

if [ $install_yandasoft == yes ]; then
# Go, go, go, yandasoft!
  if [ $casacore_version == master ]; then
	yandasoft_opts+=" -DCMAKE_CXX_FLAGS=-Dcasa=casacore"
  fi

  if [ -d build ]; then
	echo "yandasoft build directory already exists"
	cd build
	if [ -f install_manifest.txt ]; then
		make uninstall
	fi
	cd ..
	rm -rf build
	mkdir build
  else
	mkdir build
  fi
  cd build
  if [ $compiler == clang ]; then
		comp_opts="-DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  elif [ $compiler == cray ]; then
		comp_opts="-DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc"
  else
		comp_opts="-DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc"
  fi
  try ${cmake} -DCMAKE_INSTALL_PREFIX="$prefix" $comp_opts $yandasoft_opts ..
  try make -j${jobs} all
  try make -j${jobs} install
fi


if [ $remove_workdir == yes ]; then
	cd $original_dir
	rm -rf "$workdir"
fi
