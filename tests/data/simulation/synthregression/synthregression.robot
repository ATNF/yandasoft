*** Settings ***
Library             Process
Documentation       Synthesis Regresssion Test Suite 
Default Tags        synthesis  
*** Test Cases ***
wtermtest: WProject and WStack gridders
    [Template]      Run PythonTest ${thetest}
    wtermtest
pbcorrtest: single field with mosaicing gridders (primary beam correction)
    [Template]      Run PythonTest ${thetest}
    pbcorrtest
noisetest: testing the noise in the image is as expected
    [Template]      Run PythonTest ${thetest}  
    noisetest
facetingtest: test of faceted imaging with a spherical function gridder
    [Template]      Run PythonTest ${thetest}  
    facetingtest
calibratortest: test of ccalibrator
    [Template]      Run PythonTest ${thetest}  
    calibratoretest
leakagecalibtest: test of polarisation leakage calibration
    [Template]      Run PythonTest ${thetest}  
    leakagecalibtest
1934-638: test source position and flux on real ATCA data"
    [Template]      Run PythonTest ${thetest}
    test1934  
MSMFS with new imager
    [Template]      Run PythonTest ${thetest}
    testmsmfs  
*** Keywords ***
Run PythonTest ${thetest}
    Log            ${thetest}
    [Tags]         ${thetest}
    ${result} =    Run Process   python   ${thetest}.py    shell=True   stdout=./stdout
    Should Be Equal As Integers      ${result.rc}     0
