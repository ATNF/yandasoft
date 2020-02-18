add_executable(tdeconvolution tdeconvolution.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tdeconvolution 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tdeconvolution
	COMMAND tdeconvolution
	)
