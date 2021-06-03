#!/bin/bash -l 
# 
#===============================================================
# Code to extract msdata for a given beam. Particularly useful
# for data that has been split by spectral windows. 
#                                                --wr, 14May2021
#===============================================================


usage()
{
    echo "usage: sh $0 [[[-p path] [-b beam] [-w workDir] [-d dryRun] [-t tapeDir] [-o outDirName]]| [-h]]"
}

# Main script begins here: 
# Set the defaults: 
path=""
beamID=""
workDir="$(dirname $(pwd))/$(basename $(pwd))"
dryRun=1

tapeDir=""
outDirName="msdata"

while [ "$1" != "" ]; do
    case $1 in
        -p | --path )           shift
                                path="$1"
                                ;;
        -b | --beam )           shift 
		                beamID="$1"
                                ;;
        -w | --workDir )        shift 
		                workDir="$1"
                                ;;
        -t | --tapeDir )        shift 
		                tapeDir="$1"
                                ;;
        -o | --outDirName )     shift 
		                outDirName="$1"
                                ;;
        -d | --dryRun )         shift 
		                DRY_RUN="$1"
                                ;;
        -h | --help )           usage
                                exit 1
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

workDir="$(dirname $workDir)/$(basename $workDir)"
path="$(dirname $path)/$(basename $path)"
outDirName="$(basename $outDirName)"
# Sanity checks: 
#================================================================
if [ ! -e "${workDir}" ];
then
    echo "ERROR - WorkDir not found - check specified workDir: $workDir"
    usage
    exit 1
else
    cd ${workDir}
    # Append full path to workDir:
    workDir="$(dirname $(pwd))/$(basename $(pwd))"
fi

if [ "${path}" == "" ] || [ ! -e "${path}" ];
then
    echo "ERROR - Invalid path to raw ms data: $path"
    echo "pwd: $(pwd)"
    usage
    exit 1
fi
if [ "${beamID}" == "" ];
then
    echo "ERROR - Invalid Beam ID specified: $beamID"
    usage
    exit 1
fi

declare -a msList="($(ls -1a "$path" | grep ".ms"))"

nMS=${#msList[@]}
if [ "${msList}" == "" ];
then
    echo "No msdata found in: $path "
    exit 1
else
    echo "Filtering names of beam-$beamID measurement sets from a total of $nMS ms"
fi

#================================================================
# Some book-keeping and logs:
logDir="${workDir}/logs"
scriptDir="${workDir}/scripts"
mkdir -p ${logDir} ${scriptDir}

dateStr=$(date +%Y-%m-%d-%H-%M-%S)
sedstr="s/\.sh/_copyScript_${dateStr}\.sh/g"

copyScript="$(echo $0 |sed -e ${sedstr})"
copyScript="${scriptDir}/${copyScript}"


sedstr="s/\.sh/_copyLog_${dateStr}\.log/g"
copyLog="$(echo $0 |sed -e ${sedstr})"
copyLog="${logDir}/${copyLog}"
echo "Copy script: ${copyScript}"
echo "    logfile: ${copyLog}"

# Define destination directory (relative to workDir OR specify a tapeDir):
if [ "${tapeDir}" == "" ];
then
    copyCommand="module load askaputils
"
    copyMethod="mcp.sh"
    copyCommand="${copyCommand}

cd ${workDir}
# Create outDirName and change dir to outDirName
mkdir -p ${outDirName}
cd ${outDirName}
"
else
    backup_mode="tape"
    copyMethod="put"
    copyCommand="cd ${tapeDir}
mkdir ${outDirName}
cd ${outDirName}
"
fi


echo "  workDir: ${workDir}"
echo "  destDir: ${outDirName}"
echo "      log: ${copyLog}"
echo "   script: ${copyScript}"

#================================================================

# For each ms, find out the matching beamID that user specified: 
module load mstool
msOutList=""

# Initialise empty script and log
touch "${copyScript}"
touch "${copyLog}"

for (( iMS=0; iMS<$nMS; iMS++ ))
do
    beamNow=$(msInfo.py -m "${path}/${msList[$iMS]}" -q beam)
    msNow="${msList[$iMS]}"
    echo "${iMS}/${nMS}. - BeamID for ${msNow}: ${beamNow}" >>"${copyLog}"
    if [ ${beamNow} == ${beamID} ];
    then
        msOutList="${msOutList} 
${path}/${msNow}"
        if [ "${backup_mode}" != "tape" ];then 
	    if [ "${outDirName}" != "" ];
	    then
                copyCommand="${copyCommand}
${copyMethod} ${path}/${msNow} ."
            fi
	else 
            copyCommand="${copyCommand}
${copyMethod} ${path}/${msNow}"
        fi
    fi
done
echo "The following ms datsets were found for beam-$beamID :" >>"${copyLog}"
echo "${msOutList}" >>"${copyLog}"

cat > "${copyScript}" <<EOFCOPY
$copyCommand
EOFCOPY

if [ "${DRY_RUN}" == "0" ];
then
    echo "Copying data using script: ${copyScript}"
    if [ "${backup_mode}" == "tape" ];
    then
	    echo "Will archive data to tape using pshell..."
            module load askaputils 
	    pshell -i ${copyScript}
    else
	    echo "Will copy data to local disk using mcp.sh..."
	    sh ${copyScript}
    fi
else
    echo "This was a dry run. To copy, set dry run option to false. "
fi
