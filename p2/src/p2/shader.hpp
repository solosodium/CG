#ifndef _462_GLSL_SHADER_HPP_
#define _462_GLSL_SHADER_HPP_

#include "application/opengl.hpp"
#include <string>

namespace _462 {

class Shader
{
public:

    Shader();
	Shader(const char* vertex_shader, const char* fragment_shader);

    ~Shader();

	void compile_from_file(const char* vertex_shader, const char* fragment_shader);

    GLhandleARB get_program()
    {
        return program_;
    }

    void bind();

    void unbind();

protected:
	
    GLhandleARB program_;
};

} /* _462 */

#endif /* _462_GLSL_SHADER_HPP_ */

