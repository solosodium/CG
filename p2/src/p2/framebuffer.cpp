
#include "p2/framebuffer.hpp"
#include <cstring>
#include <iostream>

namespace _462 {

static const GLenum InternalFormat[Framebuffer::BUFFER_MAX] = {
    GL_RGBA8,
    GL_RGBA8,
    GL_DEPTH_COMPONENT32,
};

static const GLenum ExternalFormat[Framebuffer::BUFFER_MAX] = {
    GL_RGBA,
    GL_RGBA,
    GL_DEPTH_COMPONENT,
};

static const GLenum AttachmentType[Framebuffer::BUFFER_MAX] = {
    GL_COLOR_ATTACHMENT0_EXT,
    GL_COLOR_ATTACHMENT1_EXT,
    GL_DEPTH_ATTACHMENT_EXT,
};

// HACK just assume its 4
#define MAX_TEXTURE_UNITS 4
static const GLenum ActiveTextures[MAX_TEXTURE_UNITS] = {
    GL_TEXTURE0,
    GL_TEXTURE1,
    GL_TEXTURE2,
    GL_TEXTURE3,
};

Framebuffer::Framebuffer()
    : width_(0), height_(0), compile_needed_(true), is_rendering_(false)
{
    // clear values
    memset( texture_id_, 0, sizeof texture_id_ );
    memset( texture_number_, 0, sizeof texture_number_ );

    // create a framebuffer object
    glGenFramebuffersEXT( 1, &fbo_ );
}

Framebuffer::~Framebuffer()
{
    // will just ignore non-initialized textures
    glDeleteTextures( BUFFER_MAX, texture_id_ );
    glDeleteFramebuffersEXT( 1, &fbo_ );
}

void Framebuffer::set_size( GLuint width, GLuint height )
{
    if ( width != width_ || height != height_ ) {
        width_ = width;
        height_ = height;
        compile_needed_ = true;
    }
}

void Framebuffer::set_buffer(BufferType type, GLuint* texture_number)
{
    // if null, clear and return
    if ( !texture_number_ || *texture_number >= MAX_TEXTURE_UNITS )
        texture_number_[type] = 0;

    size_t value = *texture_number + 1;

    // make sure another type is not mapped to the same texture number
    for ( int i = 0; i < BUFFER_MAX; i++ ) {
        if ( i != type && texture_number_[type] == value )
            return;
    }

    // map type to value
    texture_number_[type] = value;

    // XXX should try to set this only if needed
    compile_needed_ = true;
}

void Framebuffer::compile()
{
    GLenum rv;

    // don't compile if we don't need to, or if we can't
    if (!compile_needed_ || width_ == 0 || height_ == 0 || is_rendering_) {
        std::cout << "Cannot compile framebuffer at this time.\n";
        return;
    }

    glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, fbo_ );

    for ( size_t i = 0; i < BUFFER_MAX; i++ ) {
        if ( texture_number_[i] ) {
            // NOTE: fglrx driver crashes if we use TEXTURE_2D, so have to use TEXTURE_RECTANGLE
            // generate a texture, only if we need a new one
            if ( !texture_id_[i] )
                glGenTextures( 1, &texture_id_[i] );

            std::cout << "creating texture " << i << "\n";

            glBindTexture( GL_TEXTURE_RECTANGLE_ARB, texture_id_[i] );
            // allocate space
            glTexImage2D(
                GL_TEXTURE_RECTANGLE_ARB, 0, InternalFormat[i],
                width_, height_, 0, ExternalFormat[i], GL_UNSIGNED_BYTE, NULL );

            glTexParameteri( GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
            glTexParameteri( GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

            if (i == BUFFER_DEPTH) {
                glTexParameteri( GL_TEXTURE_RECTANGLE_ARB, GL_DEPTH_TEXTURE_MODE, GL_LUMINANCE );
                glTexParameteri( GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_COMPARE_MODE, GL_NONE );
            }
        }

        // bind to (or detach from) the fbo
        glFramebufferTexture2DEXT( GL_FRAMEBUFFER_EXT, AttachmentType[i], GL_TEXTURE_RECTANGLE_ARB, texture_id_[i], 0 );
    }

    // check if fbo complete
    rv = glCheckFramebufferStatusEXT( GL_FRAMEBUFFER_EXT );
    switch (rv)
    {
    case GL_FRAMEBUFFER_COMPLETE_EXT:
        std::cout << "FRAMEBUFFER INITIALIZED SUCCESSFULLY\n";
        // only clear on success
        compile_needed_ = false;
        break;
    default:
        std::cout << "ERROR INITIALIZING FRAMEBUFFER: " << rv << "\n";
        break;
    }

    // clear bound texture and fbo
    glBindTexture( GL_TEXTURE_RECTANGLE_ARB, 0 );
    glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, 0 );
}

void
Framebuffer::begin_render()
{
    if ( !compile_needed_ && !is_rendering_ ) {
        glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, fbo_ );
        glPushAttrib( GL_VIEWPORT_BIT );
        glViewport( 0, 0, width_, height_ );
        is_rendering_ = true;
    }
}

void
Framebuffer::end_render()
{
    if ( is_rendering_ ) {
        glPopAttrib();
        glBindFramebufferEXT( GL_FRAMEBUFFER_EXT, 0 );
        is_rendering_ = false;
    }
}

void
Framebuffer::render()
{
    if ( is_rendering_ || compile_needed_ )
        return;

    glPushAttrib( GL_ALL_ATTRIB_BITS );
    glPushClientAttrib( GL_CLIENT_ALL_ATTRIB_BITS );

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
    glOrtho( -1.0, 1.0, -1.0, 1.0, .1, 10.0 );
    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

    glEnable( GL_TEXTURE_RECTANGLE_ARB );
    glDisable( GL_LIGHTING );
    glDisable( GL_DEPTH_TEST );

    for ( size_t i = 0; i < BUFFER_MAX; i++ ) {
        if ( texture_number_[i] ) {
            glActiveTexture( ActiveTextures[texture_number_[i] - 1] );
            glBindTexture( GL_TEXTURE_RECTANGLE_ARB, texture_id_[i] );
        }
    }
    glActiveTexture( GL_TEXTURE0 );

    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glTexCoord2f(GLfloat(width_), 0.0f);
    glVertex3f(1.0f, -1.0f, -1.0f);
    glTexCoord2f(GLfloat(width_), GLfloat(height_));
    glVertex3f(1.0f,  1.0f, -1.0f);
    glTexCoord2f(0.0f, GLfloat(height_));
    glVertex3f(-1.0f,  1.0f, -1.0f);
    glEnd();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();

    for ( size_t i = 0; i < BUFFER_MAX; i++ ) {
        if ( texture_number_[i] ) {
            glActiveTexture( ActiveTextures[texture_number_[i] - 1] );
            glBindTexture( GL_TEXTURE_RECTANGLE_ARB, 0 );
        }
    }

    glPopAttrib();
}

} /* _462 */
