add_executable(tgridding tgridding.cc)
add_executable(tgriddingadapters tgriddingadapters.cc)
include_directories(${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(tgridding 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tgridding
	COMMAND tgridding
	)

target_link_libraries(tgriddingadapters 
	askap::yandasoft
	${CPPUNIT_LIBRARY}
)
add_test(
    NAME tgriddingadapters
    COMMAND tgriddingadapters
    )
