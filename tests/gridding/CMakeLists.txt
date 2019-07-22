add_executable(tgridding tgridding.cc)
target_link_libraries(tgridding 
	yanda_synthesis
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tgridding
	COMMAND tgridding
	)
