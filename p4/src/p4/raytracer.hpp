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

namespace _462 {

/**
 * TODO: define global constants
 */
const size_t MAX_BOUNCE = 5;			// maximum bouncing number (random)

const size_t NUM_SHOOTING_RAY = 5;     // number of shooting rays
const size_t NUM_GATHERING_RAY = 5;    // number of gathering rays


class Scene;
class Ray;
struct Intersection;
class Raytracer
{
public:

    Raytracer();

    ~Raytracer();

    bool initialize( Scene* scene, size_t num_samples,
		     size_t width, size_t height );

    bool raytrace( unsigned char* buffer, real_t* max_time );

private:

    Color3 trace_pixel(const Scene* scene,
		       size_t x,
		       size_t y,
		       size_t width,
		       size_t height);

    Color3 trace(const Ray& r, unsigned int depth);

    Color3 direct_illumination(const Intersection& info);

    Color3 refraction(const Ray& r, const Intersection* info,
		      const Color3& c, unsigned int depth);

    /**
     * TODO: define new functions
     */

    // function to recursivly compute shooting ray, and get color
    Color3 shooting_ray_trace (const Ray& r, unsigned int depth, Intersection& info);

    // function to recursivly compute gathering ray, and get color
    Color3 gathering_ray_trace (const Ray& r, unsigned int depth, Intersection& info);

    // bidirectional version of refraction
    Color3 gathering_refraction(const Ray& r, Intersection& info,
              const Color3& c, unsigned int depth);

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
