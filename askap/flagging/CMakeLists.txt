#
# askap/flagging/
#
add_sources_to_yandasoft(
AmplitudeFlagger.cc
ElevationFlagger.cc
FlaggerFactory.cc
IFlagger.cc
MSFlaggingSummary.cc
SelectionFlagger.cc
StokesVFlagger.cc
)

install (FILES
AmplitudeFlagger.h
ElevationFlagger.h
FlaggerFactory.h
FlaggingStats.h
IFlagger.h
MSFlaggingSummary.h
SelectionFlagger.h
StokesVFlagger.h

DESTINATION include/askap/flagging
)
