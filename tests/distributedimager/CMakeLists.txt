add_executable(tdistributedimager tdistributedimager.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tdistributedimager
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tdistributedimager
	COMMAND tdistributedimager
	)
