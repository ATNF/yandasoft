add_executable(tdeconvolution tdeconvolution.cc)
target_link_libraries(tdeconvolution 
	yanda_synthesis
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tdeconvolution
	COMMAND tdeconvolution
	)
