add_executable(topcal topcal.cc)
target_link_libraries(topcal 
	yanda_synthesis
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME topcal
	COMMAND topcal
	)
