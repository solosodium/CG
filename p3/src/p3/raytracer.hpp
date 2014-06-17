/**
 * @file raytacer.hpp
 * @brief Raytracer class
 *
 * Implement these functions for project 3.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#ifndef _462_RAYTRACER_HPP_
#define _462_RAYTRACER_HPP_

#define MAX_DEPTH 5

#include "math/color.hpp"
#include "math/random462.hpp"

namespace _462 {

class Scene;
class Ray;
// change the name
struct Intersect;
class Raytracer
{
public:

    Raytracer();

    ~Raytracer();

    bool initialize(Scene* scene, size_t num_samples,
                    size_t width, size_t height);

    bool raytrace(unsigned char* buffer, real_t* max_time);

	/******************************************/
	/***** Added new functions interfaces *****/
	/******************************************/
	
	///////////////////
	// Key Functions //
	///////////////////
	
	// Function to cast ray to intercept
	// @param: Ray ray: the casted ray
	// @return: the color returned to this particular ray
	Color3 trace (Ray ray, size_t depth);
	
	// Function to shade a point
	// @param: Ray ray: the casted ray
	// @param: Intersect intersect: the intersect data
	// @return: the color of that point
	Color3 shade (Ray ray, Intersect intersect);
	
	/******************************************/

private:

    Color3 trace_pixel(const Scene* scene,
		       size_t x,
		       size_t y,
		       size_t width,
		       size_t height);

    // the scene to trace
    Scene* scene;

    // the dimensions of the image to trace
    size_t width, height;

    // the next row to raytrace
    size_t current_row;

    unsigned int num_samples;
};

} /* _462 */

#endif /* _462_RAYTRACER_HPP_ */
