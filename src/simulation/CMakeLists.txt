add_library( simulation OBJECT 
Simulator.cc
)
install (FILES
Simulator.h
DESTINATION include/askap/simulation
)
