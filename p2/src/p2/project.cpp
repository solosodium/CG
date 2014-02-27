#include "p2/project.hpp"
#include "scene/scene.hpp"

#include <fstream>
#define CUBEMAP_TEXTURE_NUM 2
#define BRDF_TEXTURE_NUM 4

// A namespace declaration. All proejct files use this namespace.
// Add this declration (and its closing) to all source/headers you create.
// Note that all #includes should be BEFORE the namespace declaration.
namespace _462
{

// definitions of functions for the GlslProject class

// constructor, invoked when object is created
GlslProject::GlslProject()
{
}

// destructor, invoked when object is destroyed
GlslProject::~GlslProject()
{
    // Warning: Do not throw exceptions or call virtual functions from deconstructors!
    // They will cause undefined behavior (probably a crash, but perhaps worse).
}

bool GlslProject::initialize( const SceneRenderer* renderer, int width, int height )
{
    bool rv = true;

    // copy renderer for later use
    this->renderer = renderer;
    GLint maxbuffers;
    glGetIntegerv( GL_MAX_DRAW_BUFFERS, &maxbuffers );
    std::cout << "max draw buffers: " << maxbuffers << std::endl;

    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);

    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );

    material.ambient = Color3( 1.0, 0.21, 0.2 );
    material.diffuse = Color3( 1.0, 0.43, 0.4 );
    material.specular = Color3( 1.0, 0.8, 0.4 );
    material.shininess = 30.0;
    int err;
    GLenum tex;
    GLhandleARB prog;

    material_shader = new Shader();
	// TODO: this line changes the loading of the shader
    
    // NOTE: if use a shader other than material shader, use one of the following lines
    material_shader->compile_from_file("shaders/material_vert.glsl", "shaders/material_frag.glsl");
    //material_shader->compile_from_file("shaders/outline_vert.glsl", "shaders/outline_frag.glsl");
    
    // setup the cubemap
    cubemap = Cubemap(cubemap_faces, CUBEMAP_TEXTURE_NUM);

    GLhandleARB program = material_shader->get_program();
    material_shader->bind();

    // The commented code shows how to get a cubemap into a shader. You need to make sure
    // your shader also has this variable declared, otherwise the program will not start.
    
    glUniform1iARB(err = glGetUniformLocationARB(program, "cubemap"), CUBEMAP_TEXTURE_NUM);
    //assert(err != -1);

    return rv; // return false to indicate initialization failure.
}

/**
 * Clean up the project. Free any memory, etc.
 */
void GlslProject::destroy()
{
    delete material_shader;
}

/**
 * Render the scene using shaders.
 */
void GlslProject::render()
{
    const Camera& cam = *renderer->get_camera();

    //glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glClearColor(0.7f, 0.7f, 0.7f, 1.0f);

    // render the environment map
    cubemap.render(cam);

    // render the surface of the object through material shader.
    material_shader->bind();
    GLhandleARB prog = material_shader->get_program();

    glClear(GL_DEPTH_BUFFER_BIT);
    renderer->render_scene();

    material_shader->unbind();
    glFlush();
}

} /* _462 */
