add_executable(tmeasurementequation tmeasurementequation.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tmeasurementequation 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tmeasurementequation
	COMMAND tmeasurementequation
	)
