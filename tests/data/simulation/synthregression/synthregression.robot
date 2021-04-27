*** Settings ***
Library             Process
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
SPECTRALLINE with new imager
    [Tags]          spectralline    newimager
    [Template]      Run PythonTest ${thetest}
    testspectral
REVERSED 
    [Tags]          spectralline    newimager
    [Template]      Run PythonTest ${thetest}
    test-reversed-freq

*** Keywords ***
Run PythonTest ${thetest}
    Log            ${thetest}
    [Tags]         ${thetest}
    ${result} =    Run Process   python3   ${thetest}.py    stdout=${thetest}.stdout.txt    stderr=${thetest}.stderr.txt    shell=True
    Should Be Equal As Integers      ${result.rc}     0
