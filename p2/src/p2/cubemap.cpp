#include "p2/cubemap.hpp"
#include "math/vector.hpp"
#include "math/camera.hpp"
#include <cstring>
#include <iostream>
#include <float.h>

namespace _462 {

static void rev(unsigned char* data, int tex_width, int tex_height)
{
    int sz = tex_width*4;
    unsigned char* tmp = new unsigned char[sz];
    for (int i = 0; i < tex_height/2; i++)
    {
		memmove(tmp, &data[i*sz], sz);
		memmove(&data[i*sz], &data[(tex_height-i-1)*sz], sz);
		memmove(&data[(tex_height-i-1)*sz], tmp, sz);
    }
    delete[] tmp;
}
static unsigned int
store_texture(const char* filename, unsigned int tex_num,
			  unsigned int TEX_T = GL_TEXTURE_2D)
{
    GLuint texture;
    glActiveTexture(GL_TEXTURE0 + tex_num);
    glGenTextures(1, &texture);
    glBindTexture(TEX_T, texture);

    int tex_width, tex_height;
    unsigned char* data = imageio_load_image(filename, &tex_width, &tex_height);

    if (tex_width <= 0 || tex_height <= 0)
    {
	std::cout << "Error loading texture from file \"" << filename
		  << "\"." << std::endl;
	exit(0);
    }

    rev(data, tex_width, tex_height);

    glTexImage2D(TEX_T, 0, GL_RGBA8, tex_width, tex_height,
				 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

    if (data) free(data);

    glActiveTexture(GL_TEXTURE0);
    return texture;
}

Cubemap::Cubemap(std::string* filenames, unsigned int active_tex = 5)
{
    glEnable(GL_TEXTURE_CUBE_MAP);
	glActiveTexture(GL_TEXTURE0 + active_tex);
    for (int i = 0; i < 6; i++)
    {
		int face_id = GL_TEXTURE_CUBE_MAP_POSITIVE_X + i;
		texture_ids[i] = store_texture(filenames[i].c_str(), active_tex, face_id);
    }

    // set texture wrapping behaviour
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // min/max filter setting.
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	active_texture = active_tex;

	glDisable(GL_TEXTURE_CUBE_MAP);
}

void Cubemap::render(const Camera& cam)
{
	unsigned int TEX_ID = active_texture;
	glPushAttrib( GL_ALL_ATTRIB_BITS );
    glPushClientAttrib( GL_CLIENT_ALL_ATTRIB_BITS );

    const Vector3 dir = cam.get_direction();
    const Vector3 up = cam.get_up();

    // Reset and transform the matrix.
	glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(cam.get_fov_degrees(),
				   cam.get_aspect_ratio(),
				   cam.get_near_clip(),
				   DBL_MAX);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
    gluLookAt(0,0,0, dir.x, dir.y, dir.z, up.x, up.y, up.z);

	glActiveTexture(GL_TEXTURE0 + TEX_ID);

    glEnable(GL_TEXTURE_CUBE_MAP);
	glEnable(GL_TEXTURE_2D);
    //glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    //glDisable(GL_BLEND);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    // Just in case we set all vertices to white.
    glColor4f(1,1,1,1);

    float fe = 150.0;

    glBegin(GL_QUADS);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, -1.0);
    glVertex3f(fe, -fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, -1.0);
    glVertex3f(-fe, -fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, -1.0);
    glVertex3f(-fe, fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, -1.0);
    glVertex3f(fe, fe, -fe);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, 1.0);
    glVertex3f(fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, -1.0);
    glVertex3f(fe, -fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, -1.0);
    glVertex3f(fe, fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, 1.0);
    glVertex3f(fe, fe, fe);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, 1.0);
    glVertex3f(-fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, 1.0);
    glVertex3f(fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, 1.0);
    glVertex3f(fe, fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, 1.0);
    glVertex3f(-fe, fe, fe);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, -1.0);
    glVertex3f(-fe, -fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, 1.0);
    glVertex3f(-fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, 1.0);
    glVertex3f(-fe, fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, -1.0);
    glVertex3f(-fe, fe, -fe);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, -1.0);
    glVertex3f(-fe, fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, 1.0, 1.0);
    glVertex3f(-fe, fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, 1.0);
    glVertex3f(fe, fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, 1.0, -1.0);
    glVertex3f(fe, fe, -fe);

    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, -1.0);
    glVertex3f(-fe, -fe, -fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, -1.0, -1.0, 1.0);
    glVertex3f(-fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, 1.0);
    glVertex3f(fe, -fe, fe);
    glMultiTexCoord3f(GL_TEXTURE0 + TEX_ID, 1.0, -1.0, -1.0);
    glVertex3f(fe, -fe, -fe);

    glEnd();

    // Restore enable bits and matrix

	glMatrixMode(GL_PROJECTION);
    glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glDisable(GL_TEXTURE_CUBE_MAP);

	glActiveTexture(GL_TEXTURE0);

    glPopClientAttrib();
	glPopAttrib();

}

} /* _462 */
