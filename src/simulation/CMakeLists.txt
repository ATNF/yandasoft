add_library( simulation OBJECT 
Simulator.cc
)

set_property(TARGET simulation PROPERTY POSITION_INDEPENDENT_CODE ON)

install (FILES
Simulator.h
DESTINATION include/askap/simulation
)
