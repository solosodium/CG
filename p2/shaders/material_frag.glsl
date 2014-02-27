varying vec3 norm;
varying vec3 cam_dir;
varying vec3 color;

// Declare any additional variables here. You may need some uniform variables.

// TODO: declare the cubemap variable [samplerCube]
uniform samplerCube cubemap;

void main(void)
{
	// Your shader code here.
	// Note that this shader won't compile since gl_FragColor is never set.
	
	// compute the reflected ray
	vec3 ref_ray = reflect(cam_dir, normalize(norm));
	
	// set the color from the cubemap
	gl_FragColor = textureCube(cubemap, ref_ray);
}
