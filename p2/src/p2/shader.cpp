#include "p2/shader.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace _462 {

static char* load_file(const char* file)
{
    std::ifstream infile;
    infile.open(file);

    if(infile.fail()){
        std::cout << "ERROR: cannot open file " << file << std::endl;
        infile.close();
    }

    infile.seekg(0, std::ios::end );
    int length = (int)infile.tellg();
    infile.seekg(0, std::ios::beg );

    char* buffer = new char[length];
    infile.getline(buffer, length, '\0');

    infile.close();

    return buffer;
}

/**
 * Load a file as either a vertex shader or a fragment shader, and attach
 * it to a program
 * @param file  The file to load
 * @param type  Either GL_VERTEX_SHADER_ARB, or GL_FRAGMENT_SHADER_ARB
 * @param program  The shading program to which the shaders are attached
 */
static void load_shader(const char* file, GLint type, GLhandleARB program)
{
    int result;
    char error_msg[2048];

    const char* src = load_file(file);
    // Create shader object
    GLhandleARB shader = glCreateShaderObjectARB(type);

    // Load Shader Sources
    glShaderSourceARB(shader, 1, &src, NULL);
    // Compile The Shaders
    glCompileShaderARB(shader);
    // Get compile result
    glGetObjectParameterivARB(shader, GL_OBJECT_COMPILE_STATUS_ARB, &result);
    if (!result)
	{
        glGetInfoLogARB(shader, sizeof(error_msg), NULL, error_msg);
        std::cout << "GLSL COMPILE ERROR(" << file << "): " << error_msg << std::endl;
        exit(2);
    }
	else std::cout << "Compiled shaders successfully" << std::endl;

    // Attach The Shader Objects To The Program Object
    glAttachObjectARB(program, shader);

    delete[] src;
}

static void load_shaders(GLhandleARB program, const char* vert_file,
						 const char* frag_file)
{
    // Load vertex shader
    std::cout << "loading vertex shader " << vert_file << std::endl;
    load_shader( vert_file, GL_VERTEX_SHADER_ARB, program );

    // Load fragment shader
    std::cout << "loading fragment shader " << frag_file << std::endl;
    load_shader( frag_file, GL_FRAGMENT_SHADER_ARB, program );

    glLinkProgramARB(program);

    GLint rv;
    glGetProgramiv(program, GL_LINK_STATUS, &rv);
    if (rv == GL_TRUE) std::cout << "Successfully linked shader" << std::endl;
    else std::cout << "FAILED to link shader" << std::endl;
}

Shader::Shader()
{
    // Create shader program
    program_  = glCreateProgramObjectARB();
}

Shader::~Shader()
{
    glDeleteObjectARB(program_);
}

void Shader::compile_from_file(const char* vertex_shader,
							   const char* fragment_shader)
{
    load_shaders(program_, vertex_shader, fragment_shader);
}

void Shader::bind()
{
    glUseProgramObjectARB(program_);
}

void Shader::unbind()
{
    glUseProgramObjectARB(0);
}

} /* _462 */
