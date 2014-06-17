add_executable(p5 collisions.cpp main.cpp physics.cpp planebody.cpp
				spherebody.cpp spring.cpp trianglebody.cpp)
target_link_libraries(p5 application math scene tinyxml ${SDL_LIBRARY}
                      ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}
                      ${GLEW_LIBRARIES})

if(APPLE)
    target_link_libraries(p5 SDLmain)
endif()

install(TARGETS p5 DESTINATION ${PROJECT_SOURCE_DIR}/..)
