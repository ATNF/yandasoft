add_executable(topcal topcal.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(topcal 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME topcal
	COMMAND topcal
	)
