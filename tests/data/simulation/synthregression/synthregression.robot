*** Settings ***
Library             Process
Library             DateTime
Documentation       Synthesis Regresssion Test Suite 
Default Tags        synthesis  

*** Test Cases ***
MSMFS with new imager
    [Tags]          testmsmfs      newimager
    [Template]      Run PythonTest ${thetest}
    testmsmfs  
MSMFS with channels gridded together
    [Tags]          combinedchannels    newimager
    [Template]      Run PythonTest ${thetest}
    testcombined
MSMFS with reversed channels gridded together
    [Tags]          combinedchannels    newimager
    [Template]      Run PythonTest ${thetest}
    testcombined-reversed
SPECTRALLINE with new imager
    [Tags]          spectralline    newimager
    [Template]      Run PythonTest ${thetest}
    testspectral
REVERSED 
    [Tags]          spectralline    newimager
    [Template]      Run PythonTest ${thetest}
    test-reversed-freq
wtermtest: WProject and WStack gridders
    [Tags]          wtermtest
    [Template]      Run PythonTest ${thetest}
    wtermtest
imageroutputtest: test selection of imager outputs
    [Tags]          testimageroutputs
    [Template]      Run PythonTest ${thetest}
    testimageroutputs
pbcorrtest: single field with mosaicing gridders (primary beam correction)
    [Tags]          pbcorrtest
    [Template]      Run PythonTest ${thetest}
    pbcorrtest
noisetest: testing the noise in the image is as expected
    [Tags]          noisetest
    [Template]      Run PythonTest ${thetest}  
    noisetest
facetingtest: test of faceted imaging with a spherical function gridder
    [Tags]          facetingtest
    [Template]      Run PythonTest ${thetest}  
    facetingtest
calibratortest: test of ccalibrator
    [Tags]          calibratortest    non-critical
    [Template]      Run PythonTest ${thetest}  
    calibratortest
leakagecalibtest: test of polarisation leakage calibration
    [Tags]          leakagecalibtest   non-critical
    [Template]      Run PythonTest ${thetest}  
    leakagecalibtest
1934-638: test source position and flux on real ATCA data"
    [Tags]          test1934  
    [Template]      Run PythonTest ${thetest}
    test1934  
ddcal: test DDCal calibration application in imager
    [Tags]          ddcal   newimager
    [Template]      Run PythonTest ${thetest}
    ddcal
delaysolver: test delaysolver app on real ATCA data"
    [Tags]          testdelaysolver
    [Template]      Run PythonTest ${thetest}
    testdelaysolver
mfsmodel: start imager in spectral mode with mfs model
    [Tags]          spectralline   newimager
    [Template]      Run PythonTest ${thetest}
    testmfsmodel
restart: restart imager in spectral mode
    [Tags]          spectralline   newimager
    [Template]      Run PythonTest ${thetest}
    testrestart
extended: compare standard and pixellist clean on extended source
    [Tags]          pixellist   newimager
    [Template]      Run PythonTest ${thetest}
    testcleanextended
*** Keywords ***
Run PythonTest ${thetest}
    ${start_time} =    Get Current Date
    Log To Console    ${thetest}
    [Tags]          ${thetest}
    ${result} =    Run Process   python3 ${thetest}.py    stdout=${thetest}.stdout.txt    stderr=${thetest}.stderr.txt    shell=True
    Should Be Equal As Integers      ${result.rc}     0
    ${end_time} =    Get Current Date
    ${elapsed_time} =    Subtract Date From Date  ${end_time}  ${start_time}
    Log To Console    Elapsed time: ${elapsed_time} seconds
