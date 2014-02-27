add_executable(p2 main.cpp project.cpp shader.cpp framebuffer.cpp cubemap.cpp)
target_link_libraries(p2 application math scene tinyxml ${SDL_LIBRARY}
	${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})

if(APPLE)
	target_link_libraries(p2 SDLmain)
endif()

install(TARGETS p2 DESTINATION ${PROJECT_SOURCE_DIR}/..)
