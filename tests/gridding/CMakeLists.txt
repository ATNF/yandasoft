add_executable(tgridding tgridding.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tgridding 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tgridding
	COMMAND tgridding
	)
