add_library( messages OBJECT 
ContinuumWorkRequest.cc
ContinuumWorkUnit.cc
IMessage.cc
SpectralLineWorkRequest.cc
SpectralLineWorkUnit.cc
)
install (FILES
ContinuumWorkRequest.h
ContinuumWorkUnit.h
IMessage.h
SpectralLineWorkRequest.h
SpectralLineWorkUnit.h
DESTINATION include/askap/messages
)
