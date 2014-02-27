
#ifndef _462_GLSL_FRAMEBUFFER_HPP_
#define _462_GLSL_FRAMEBUFFER_HPP_

#include "application/opengl.hpp"

namespace _462 {

class Framebuffer
{
public:

    enum BufferType
    {
        BUFFER_COLOR_0 = 0,
        BUFFER_COLOR_1,
        BUFFER_DEPTH
    };

    enum { BUFFER_MAX = 3 };

    Framebuffer();

    ~Framebuffer();

    /**
     * Sets the size of the buffer. A recompile will be needed if size changes.
     */
    void set_size( GLuint width, GLuint height );

    /**
     * Sets the active texture number for a given buffer. If NULL, the given
     * buffer is not used. If non-null, uses that texture for the number.
     * A recompile will be needed.
     */
    void set_buffer( BufferType type, GLuint* texture_number);

    /**
     * Compiles the framebuffer if needed. May be invoked multiple times.
     */
    void compile();

    /**
     * Begin rendering to the framebuffer.
     * @remark Does nothing if a compile is needed or already in a render.
     */
    void begin_render();

    /**
     * Finish rendering to the framebuffer.
     * @remark Does nothing if a compile is needed or not in a render.
     */
    void end_render();

    /**
     * Binds all textures and renders the scene as a single quad taking up the
     * entire viewport.
     * @remark Does nothing if a compile is needed or already in a render.
     */
    void render();

private:

    /**
     * Stores handle for the texture object for each buffer.
     * Created lazily on compile as needed.
     */
    GLuint texture_id_[BUFFER_MAX];

    /**
     * The active texture to use for the given buffer, plus 1. Zero if not
     * to be used.
     */
    GLuint texture_number_[BUFFER_MAX];

    /**
     * Handle to the frame buffer
     */
    GLuint fbo_;

    /**
     * Current dimensions of the buffer.
     */
    GLuint width_, height_;

    /**
     * True if the buffer needs to be recompiled before use. Can be caused
     * by things such as resizing the buffer or changing the attached buffers.
     */
    bool compile_needed_;

    /**
     * True if currently in a render; begin_render has been invoked but not end_redner.
     */
    bool is_rendering_;
};

} /* _462 */

#endif /* _462_GLSL_FRAMEBUFFER_HPP_ */

