add_executable(tsynthesisutils tsynthesisutils.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tsynthesisutils
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tsynthesisutils
	COMMAND tsynthesisutils
	)
