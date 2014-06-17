/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * Implement these functions for project 4.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "raytracer.hpp"
#include "scene/scene.hpp"

#include <SDL_timer.h>
#include <iostream>
#include <random>

#ifdef OPENMP // just a defense in case OpenMP is not installed.

#include <omp.h>

#endif
namespace _462 {

// max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 8

static const unsigned STEP_SIZE = 8;

Raytracer::Raytracer()
    : scene(0), width(0), height(0) { }

// random real_t in [0, 1)
static inline real_t random()
{
    return real_t(rand())/RAND_MAX;
}

Raytracer::~Raytracer() { }

/**
 * Initializes the raytracer for the given scene. Overrides any previous
 * initializations. May be invoked before a previous raytrace completes.
 * @param scene The scene to raytrace.
 * @param width The width of the image being raytraced.
 * @param height The height of the image being raytraced.
 * @return true on success, false on error. The raytrace will abort if
 *  false is returned.
 */
bool Raytracer::initialize(Scene* scene, size_t num_samples,
			   size_t width, size_t height)
{
    /*
     * omp_set_num_threads sets the maximum number of threads OpenMP will
     * use at once.
     */
#ifdef OPENMP
    omp_set_num_threads(MAX_THREADS);
#endif
    this->scene = scene;
    this->num_samples = num_samples;
    this->width = width;
    this->height = height;

    current_row = 0;

    Ray::init(scene->camera);
    scene->initialize();

    const SphereLight* lights = scene->get_lights();

    // TODO any initialization or precompuation before the trace

    return true;
}

/************************************************************************/
/***************** Added new functions implementations ******************/
/************************************************************************/

///////////////////
// Key Functions //
///////////////////

// Function to cast ray to intercept
// @param: Ray ray: the casted ray
// @return: the color returned to this particular ray
Color3 Raytracer::trace(Ray ray, size_t depth)
{
	// base case for trace
	if (depth > MAX_DEPTH) {
		return scene->background_color;
	}
	
	// get all the geometries in the scene
	Geometry* const* geometries = scene->get_geometries();
	const size_t num_geometries = scene->num_geometries();

	// need to memorize the closest distance
	real_t min_dis = 999999.0;			// minimum distance
	Intersect intersect;				// min dis intersection
	intersect.intersect = false;		// defalut state is no intersection found
    
    /**
     * Potntially optimizable
     */
	// iterate through the geometries
	for (size_t k=0; k<num_geometries; k++) {
        // test geometry with intersection test
        Intersect inter_temp = geometries[k]->intersect(ray);
        // if intersect
        if (inter_temp.intersect) {
			real_t dis_temp = distance(ray.e, inter_temp.position);
			// if closer intersect point found
			if (dis_temp < min_dis) {
				min_dis = dis_temp;
				intersect = inter_temp;
			}
        }
	}

	// check if intersection found
	if (intersect.intersect) 
    {
        // from here it needs to cast two more rays depends on the material property
        if (intersect.refractive_index < EPSILON)
        {
            // Case: refractive index = 0 means it is opague, need reflective ray
            // first check if reflective is possible
            if (dot(normalize(-ray.d), intersect.normal) < 0.0) {
                // no reflection possible
                return shade(ray, intersect);
            }
            else {
                // get reflection ray
                Vector3 ref_d = ray.d - intersect.normal * dot(ray.d, intersect.normal) * 2.0;
                Ray ref = Ray(intersect.position, normalize(ref_d));
                // cast ray
                return shade(ray, intersect) + intersect.specular * trace(ref, depth+1);
            }
        }
        else
        {
            // Case: refractive index is something else, consider refraction

            // Decide whether it is shooting in or shooting out
            if (dot(normalize(-ray.d), intersect.normal) > 0.0)
            {
                // this ray comes from outside of the object

                // Step 0: calculate Fresnel coefficient (Fundie book pp.305)
                real_t R0 = (intersect.refractive_index - 1.0) * (intersect.refractive_index - 1.0) / 
                            (intersect.refractive_index + 1.0) / (intersect.refractive_index + 1.0);
                real_t cos_theta = dot(normalize(-ray.d), intersect.normal);
                real_t R = R0 + (1.0-R0)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta);

                /*
                // Step 1: calcualte reflect ray: incoming = ray; normal = intersect.normal
                Vector3 reflect_d = ray.d - intersect.normal * dot(ray.d, intersect.normal) * 2.0;
                Ray reflect = Ray(intersect.position, normalize(reflect_d));

                // Step 2: calculate refract ray: incoming = ray; normal = intersect.normal
                real_t n1_n2 = scene->refractive_index / intersect.refractive_index;
                Vector3 refract_d = n1_n2*ray.d + ( n1_n2*cos_theta - sqrt(1.0 - n1_n2*n1_n2*(1.0-cos_theta*cos_theta)) )*intersect.normal;
                Ray refract = Ray(intersect.position, normalize(refract_d));

                // Step 3: cast two rays
                return intersect.specular * R * trace(reflect, depth+1) + (1.0 - R) * trace(refract, depth+1);
                */

                
                // Fix: choose a random ray
                real_t random_number = random();

                if (random_number < R) {
                    // reflect ray only
                    
                    // Step 1: calcualte reflect ray: incoming = ray; normal = intersect.normal
                    Vector3 reflect_d = ray.d - intersect.normal * dot(ray.d, intersect.normal) * 2.0;
                    Ray reflect = Ray(intersect.position, normalize(reflect_d));

                    // Step 3: cast ray
                    return intersect.specular * R * trace(reflect, depth+1);
                }
                else {
                    // refract ray only

                    // Step 2: calculate refract ray: incoming = ray; normal = intersect.normal
                    real_t n1_n2 = scene->refractive_index / intersect.refractive_index;
                    Vector3 refract_d = n1_n2*ray.d + ( n1_n2*cos_theta - sqrt(1.0 - n1_n2*n1_n2*(1.0-cos_theta*cos_theta)) )*intersect.normal;
                    Ray refract = Ray(intersect.position, normalize(refract_d));

                    // Step 3: cast ray
                    return (1.0 - R) * trace(refract, depth+1);
                }
                
            }
            else
            {
                // this ray comes from inside of this object

                // check if it is total internal reflection
                real_t cos_theta = dot(normalize(-ray.d), -intersect.normal);    // cos theta

                if (intersect.refractive_index / scene->refractive_index * sqrt(1.0-cos_theta*cos_theta) >= 1.0) 
                {
                    // total internal reflection detects

                    // Step 1: calcualte reflect ray: incoming = ray; normal = -intersect.normal
                    Vector3 reflect_d = ray.d + intersect.normal * dot(ray.d, -intersect.normal) * 2.0;
                    Ray reflect = Ray(intersect.position, normalize(reflect_d));

                    // Step 2: cast ray
                    return trace(reflect, depth+1);
                }
                else 
                {
                    // Step 0: calculate Fresnel coefficient (Fundie book pp.305)
                    real_t R0 = (scene->refractive_index - 1.0) * (scene->refractive_index - 1.0) / 
                                (scene->refractive_index + 1.0) / (scene->refractive_index + 1.0);
                    real_t R = R0 + (1.0-R0)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta)*(1.0-cos_theta);

                    /*
                    // Step 1: calcualte reflect ray: incoming = ray; normal = intersect.normal
                    Vector3 reflect_d = ray.d + intersect.normal * dot(ray.d, -intersect.normal) * 2.0;
                    Ray reflect = Ray(intersect.position, normalize(reflect_d));

                    // Step 2: calculate refract ray: incoming = ray; normal = intersect.normal
                    real_t n1_n2 =  intersect.refractive_index / scene->refractive_index;
                    Vector3 refract_d = n1_n2*ray.d + ( n1_n2*cos_theta - sqrt(1.0 - n1_n2*n1_n2*(1.0-cos_theta*cos_theta)) )*(-intersect.normal);
                    Ray refract = Ray(intersect.position, normalize(refract_d));

                    // Step 3: cast ray
                    return R * trace(reflect, depth+1) + (1.0 - R) * trace(refract, depth+1);
                    */

                    
                    // Fix: choose a random ray
                    real_t random_number = random();
                    if (random_number < R) {
                        // reflect ray only
                        
                        // Step 1: calcualte reflect ray: incoming = ray; normal = intersect.normal
                        Vector3 reflect_d = ray.d + intersect.normal * dot(ray.d, -intersect.normal) * 2.0;
                        Ray reflect = Ray(intersect.position, normalize(reflect_d));

                        // Step 3: cast ray
                        return R * trace(reflect, depth+1);
                    }
                    else {
                        // refract ray only

                        // Step 2: calculate refract ray: incoming = ray; normal = intersect.normal
                        real_t n1_n2 =  intersect.refractive_index / scene->refractive_index;
                        Vector3 refract_d = n1_n2*ray.d + ( n1_n2*cos_theta - sqrt(1.0 - n1_n2*n1_n2*(1.0-cos_theta*cos_theta)) )*(-intersect.normal);
                        Ray refract = Ray(intersect.position, normalize(refract_d));

                        // Step 3: cast ray
                        return (1.0 - R) * trace(refract, depth+1);
                    }
                    
                }
            }
        }
	}
	else {
		return scene->background_color;
	}
}

// Function to shade a point
// @param: Ray ray: the casted ray
// @param: Intersect intersect: the intersect data
// @return: the color of that point
Color3 Raytracer::shade (Ray ray, Intersect intersect)
{
    //-----------------------------------------------
	// compute the ambient color
	Color3 ambient = intersect.ambient * scene->ambient_light;
	
    //-----------------------------------------------
	// compute the diffuse color
	Color3 diffuse = Color3::Black();
	// get all the light sources in the scene
	const SphereLight* lights = scene->get_lights();
	size_t num_lights = scene->num_lights();
	// iterate through all the lights
	for (size_t k=0; k<num_lights; k++) 
	{
		// Check if an object is in the way (shadow ray cast)
        size_t inter_acc = NUM_SHADOW_RAY;  // the accumulator of number of blocks
		// this part is Monte Carlo test for the sphere light
        for (size_t n=0; n<NUM_SHADOW_RAY; n++)
        {
            // find a random point in the light sphere
            real_t x, y, z;
            do {
            x = random_gaussian();
            y = random_gaussian();
            z = random_gaussian();
            } while (x*x + y*y + z*z > 1.0);
            Vector3 rand_pos = normalize(Vector3(x-0.5, y-0.5, z-0.5));
            rand_pos = rand_pos * lights[k].radius + lights[k].position;
            // cast a ray to check block
            Ray light_ray = Ray(intersect.position, normalize(rand_pos-intersect.position));
            /**
             * Potntially optimizable
             */
            // get hold on to all the geometries
            Geometry* const* geometries = scene->get_geometries();
            const size_t num_geometries = scene->num_geometries();
            // iterate through the geometries
            for (size_t m=0; m<num_geometries; m++) {
                Intersect it = geometries[m]->intersect(light_ray);
                if ( it.intersect && (distance(intersect.position, lights[k].position) > distance(intersect.position, it.position)) ) {
                    inter_acc -= 1;
                    break;
                }
            }
        }
        // calculate the blockage factor
        real_t block_factor = ((real_t)(inter_acc))/((real_t)(NUM_SHADOW_RAY));
		// calculate the light attenuation
        real_t dis = length(lights[k].position - intersect.position);    // distance from light to the intersection point
        // the light attenuation factor
        real_t att_factor = 1.0 / (lights[k].attenuation.constant + lights[k].attenuation.linear*dis + lights[k].attenuation.quadratic*dis*dis);
        
        // ABUG ?

        // calculate the diffuse factor using Bling Phong shading model
        real_t diff_factor = dot(normalize(lights[k].position-intersect.position), normalize(intersect.normal));
        // adjust diff_factor for smaller than 0.0 case
        diff_factor = (diff_factor < 0.0) ? 0.0 : diff_factor;

        // cumulate the color to the diffuse component
        diffuse += intersect.diffuse * block_factor * att_factor * diff_factor * lights[k].color;
	}

    //-----------------------------------------------
	// return the value
	return ( ambient + diffuse ) * intersect.texture;
}

/************************************************************************/

/**
 * Performs a raytrace on the given pixel on the current scene.
 * The pixel is relative to the bottom-left corner of the image.
 * @param scene The scene to trace.
 * @param x The x-coordinate of the pixel to trace.
 * @param y The y-coordinate of the pixel to trace.
 * @param width The width of the screen in pixels.
 * @param height The height of the screen in pixels.
 * @return The color of that pixel in the final image.
 */
Color3 Raytracer::trace_pixel(const Scene* scene,
			      size_t x,
			      size_t y,
			      size_t width,
			      size_t height)
{
    assert(x < width);
    assert(y < height);
    real_t dx = real_t(1)/width;
    real_t dy = real_t(1)/height;

    Color3 res = Color3::Black();
    unsigned int iter;
    
    for (iter = 0; iter < num_samples; iter++)
    {
        // pick a point within the pixel boundaries to fire our
        // ray through.
        real_t i = real_t(2)*(real_t(x)+random())*dx - real_t(1);
        real_t j = real_t(2)*(real_t(y)+random())*dy - real_t(1);

        Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));

        // TODO return the color of the given pixel
        // you don't have to use this stub function if you prefer to
        // write your own version of Raytracer::raytrace.

		// TESTING
		res += trace(r, 0);

    }
    return res*(real_t(1)/num_samples);
}

/**
 * Raytraces some portion of the scene. Should raytrace for about
 * max_time duration and then return, even if the raytrace is not copmlete.
 * The results should be placed in the given buffer.
 * @param buffer The buffer into which to place the color data. It is
 *  32-bit RGBA (4 bytes per pixel), in row-major order.
 * @param max_time, If non-null, the maximum suggested time this
 *  function raytrace before returning, in seconds. If null, the raytrace
 *  should run to completion.
 * @return true if the raytrace is complete, false if there is more
 *  work to be done.
 */
bool Raytracer::raytrace(unsigned char* buffer, real_t* max_time)
{
    // TODO Add any modifications to this algorithm, if needed.

    static const size_t PRINT_INTERVAL = 64;

    // the time in milliseconds that we should stop
    unsigned int end_time = 0;
    bool is_done;

    if (max_time)
    {
        // convert duration to milliseconds
        unsigned int duration = (unsigned int) (*max_time * 1000);
        end_time = SDL_GetTicks() + duration;
    }

    // until time is up, run the raytrace. we render an entire group of
    // rows at once for simplicity and efficiency.
    for (; !max_time || end_time > SDL_GetTicks(); current_row += STEP_SIZE)
    {
        // we're done if we finish the last row
        is_done = current_row >= height;
        // break if we finish
        if (is_done) break;

        int loop_upper = std::min(current_row + STEP_SIZE, height);

        // This tells OpenMP that this loop can be parallelized.
#pragma omp parallel for
        for (int c_row = current_row; c_row < loop_upper; c_row++)
        {
            /*
             * This defines a critical region of code that should be
             * executed sequentially.
             */
#pragma omp critical
            {
                if (c_row % PRINT_INTERVAL == 0)
                    printf("Raytracing (Row %d)\n", c_row);
            }

            for (size_t x = 0; x < width; x++)
            {
                // trace a pixel
                Color3 color = trace_pixel(scene, x, c_row, width, height);

                // write the result to the buffer, always use 1.0 as the alpha
                color.to_array(&buffer[4 * (c_row * width + x)]);

                
            }
        }
    }

    if (is_done) printf("Done raytracing!\n");

    return is_done;
}

} /* _462 */
