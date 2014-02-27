
#ifndef _462_GLSL_CUBEMAP_HPP_
#define _462_GLSL_CUBEMAP_HPP_

#include "application/opengl.hpp"
#include "application/imageio.hpp"
#include "math/camera.hpp"
namespace _462 {

class Cubemap
{
public:
	Cubemap() {}
	Cubemap(std::string* filenames, unsigned int active_texture);

	void render(const Camera& cam);

	GLuint texture_ids[6];
	GLuint active_texture;
};

} /* _462 */

#endif /* _462_GLSL_FRAMEBUFFER_HPP_ */

