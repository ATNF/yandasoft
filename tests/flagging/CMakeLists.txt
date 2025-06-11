add_executable(tcflagger tcflagger.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tcflagger
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tcflagger
	COMMAND tcflagger
	)
