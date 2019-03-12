add_library( messages OBJECT 
ContinuumWorkRequest.cc
ContinuumWorkUnit.cc
IMessage.cc
SpectralLineWorkRequest.cc
SpectralLineWorkUnit.cc
)

set_property(TARGET messages PROPERTY POSITION_INDEPENDENT_CODE ON)

install (FILES
ContinuumWorkRequest.h
ContinuumWorkUnit.h
IMessage.h
SpectralLineWorkRequest.h
SpectralLineWorkUnit.h
DESTINATION include/askap/messages
)
