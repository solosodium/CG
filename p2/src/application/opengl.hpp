/**
 * @file opengl.hpp
 * @brief Includes the correct opengl header files depending on the
 *  platform. Use this file to include any gl header files.
 */

#ifndef _462_OPENGL_HPP_
#define _462_OPENGL_HPP_

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <GL/glew.h>
#define NO_SDL_GLEXT
#include <SDL_opengl.h>

#endif /* _462_OPENGL_HPP_ */

