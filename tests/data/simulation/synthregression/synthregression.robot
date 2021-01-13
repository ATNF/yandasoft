*** Settings ***
Library             Process
Documentation       Synthesis Regresssion Test Suite 
Default Tags        synthesis  
*** Test Cases ***
wtermtest: WProject and WStack gridders
    [Tags]          wtermtest
    [Template]      Run PythonTest ${thetest}
    wtermtest
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
MSMFS with new imager
    [Tags]          testmsmfs    non-critical
    [Template]      Run PythonTest ${thetest}
    testmsmfs  
*** Keywords ***
Run PythonTest ${thetest}
    Log            ${thetest}
    [Tags]         ${thetest}
    ${result} =    Run Process   python   ${thetest}.py    shell=True   stdout=${thetest}.stdout    stderr=${thetest}.stderr
    Should Be Equal As Integers      ${result.rc}     0
