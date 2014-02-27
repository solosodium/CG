include_directories(
	${PROJECT_SOURCE_DIR}
	${SDL_INCLUDE_DIR}
)
add_executable(p1 main.cpp project.cpp)
target_link_libraries(p1 math application scene ${SDL_LIBRARY} ${PNG_LIBRARIES}
	${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})

if(APPLE)
	target_link_libraries(p1 SDLmain)
endif()

install(TARGETS p1 DESTINATION ${PROJECT_SOURCE_DIR}/..)
