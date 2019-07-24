add_executable(tmeasurementequation tmeasurementequation.cc)
target_link_libraries(tmeasurementequation 
	yanda_synthesis
	${CPPUNIT_LIBRARY}
)
add_test(
	NAME tmeasurementequation
	COMMAND tmeasurementequation
	)
