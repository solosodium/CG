/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * Project 3 solution code. DO NOT DISTRIBUTE.
 *
 * @author Harry. Q. Bovik (hqbovik)
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

static inline real_t normal()
{
    static std::default_random_engine generator;
    static std::normal_distribution<real_t> dist =
      std::normal_distribution<real_t>();
    return dist(generator);
}

static inline Vector3 point_on_sphere()
{
    Vector3 r;
    r.x = normal();
    r.y = normal();
    r.z = normal();
    return normalize(r);
}

// returns a vector on the uniform hemisphere.
static inline Vector3 cos_weighted_hemi(const Vector3& n)
{
    real_t Xi1 = random();
    real_t Xi2 = random();

    real_t theta = std::acos(std::sqrt(1.0-Xi1));
    real_t phi = 2.0 * M_PI * Xi2;

    real_t xs = sin(theta) * cos(phi);
    real_t ys = cos(theta);
    real_t zs = sin(theta) * sin(phi);

    Vector3 y(n.x, n.y, n.z);
    Vector3 h = y;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) {
        h.x= 1.0;
    } else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) {
        h.y= 1.0;
    } else {
        h.z= 1.0;
    }

    Vector3 x = normalize(cross(h, y));
    Vector3 z = normalize(cross(x, y));

    Vector3 dir = xs * x + ys * y + zs * z;
    return normalize(dir);
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

	/**
	 * TODO: BVH initialization should be here
	 */
	

    return true;
}

static inline Vector3 reflect(const Vector3& I, const Vector3& N)
{
    return normalize(I - real_t(2.0)*dot(N, I)*N);
}

inline Color3 Raytracer::direct_illumination(const Intersection& info)
{
  const SphereLight* lights = scene->get_lights();
  unsigned int N = scene->num_lights();

  Color3 res = Color3::Black();
  Intersection tinfo;
  for (unsigned int i = 0; i < N; i++)
  {
    tinfo.t = DBL_MAX;
    Vector3 lightDir = (lights[i].position +
        point_on_sphere()*lights[i].radius) - info.pos;
    real_t dist = length(lightDir);

    Ray shadow(info.pos + lightDir*0.00001, lightDir);
    if (!scene->intersect(shadow, &tinfo) || tinfo.t > real_t(1))
    {
      lightDir /= dist;

      real_t weight = std::max(dot(info.normal, lightDir), real_t(0.0));

      res += lights[i].color*info.diffuse*weight;
    }
  }
  return res;
}
static inline bool refract(const Vector3& d, const Vector3& n, real_t d_n,
                           real_t refI, Vector3* t)
{
    // create refraction ray
    real_t disc = 1.0 - refI*refI*(1 - d_n*d_n);
    if (disc < 0.0)
        return false; // total internal reflection

    real_t rad = sqrt(disc);

    if (d_n > 0.0)
        *t = refI*(d - n*d_n) + n*rad;
    else
        *t = refI*(d - n*d_n) - n*rad;

    return true;
}

Color3 Raytracer::refraction(const Ray& r, const Intersection* info,
			     const Color3& c, unsigned int depth)
{
    if (depth <= 0)
      return Color3::Black();

    real_t d_n = dot(r.d, info->normal);

    real_t objRI = info->refractive_index,
           refIndex = objRI/scene->refractive_index;

    if (d_n < 0.0)
    {
        refIndex = 1.0/refIndex; // inside object.. kinda hacky.
    }

    Vector3 t;

    Color3 rest = Color3::Black();
    real_t R = 1.0;
    Color3 cc = c;

    if (refract(r.d, info->normal, d_n, refIndex, &t)) {
      // do refraction calculations and fresnel effect.

      Ray nr = Ray(r.e + 0.00001*t, t);
      real_t R0 = (objRI-1.0)*(objRI-1.0)/((objRI + 1.0)*(objRI + 1.0));
      real_t f = 1 - fabs(d_n);
      real_t g = f*f;
      g*=g;

      R = R0 + (1.0 - R0)*g*f;

      if (random() < R) {
        Ray reflected;
        reflected.d = reflect(r.d, info->normal);
        reflected.e = info->pos + 0.0001*reflected.d;
        Color3 reflColor = trace(reflected, depth-1);
        cc = info->specular*reflColor;
        return cc;
      } else {
        return trace(nr, depth-1);
      }
    } else {
      Vector3 newDir = reflect(r.d, info->normal);
      Ray newr = Ray(r.e + 0.00001*newDir, newDir);
      cc = trace(newr, depth-1);
    }
    return (R*cc + rest);
}

Color3 Raytracer::gathering_refraction(const Ray& r, Intersection& info,
              const Color3& c, unsigned int depth)
{
    if (depth <= 0)
      return Color3::Black();

    real_t d_n = dot(r.d, info.normal);

    real_t objRI = info.refractive_index,
           refIndex = objRI/scene->refractive_index;

    if (d_n < 0.0)
    {
        refIndex = 1.0/refIndex; // inside object.. kinda hacky.
    }

    Vector3 t;

    Color3 rest = Color3::Black();
    real_t R = 1.0;
    Color3 cc = c;

    if (refract(r.d, info.normal, d_n, refIndex, &t)) {
      // do refraction calculations and fresnel effect.

      Ray nr = Ray(r.e + 0.00001*t, t);
      real_t R0 = (objRI-1.0)*(objRI-1.0)/((objRI + 1.0)*(objRI + 1.0));
      real_t f = 1 - fabs(d_n);
      real_t g = f*f;
      g*=g;

      R = R0 + (1.0 - R0)*g*f;

      if (random() < R) {
        Ray reflected;
        reflected.d = reflect(r.d, info.normal);
        reflected.e = info.pos + 0.0001*reflected.d;
        Color3 reflColor = gathering_ray_trace(reflected, depth-1, info);
        cc = info.specular*reflColor;
        return cc;
      } else {
        return gathering_ray_trace(nr, depth-1, info);
      }
    } else {
      Vector3 newDir = reflect(r.d, info.normal);
      Ray newr = Ray(r.e + 0.00001*newDir, newDir);
      cc = gathering_ray_trace(newr, depth-1, info);
    }
    return (R*cc + rest);
}

Color3 Raytracer::trace(const Ray& r, unsigned int depth)
{
    if (depth <= 0) return Color3::Black();

    Intersection info;

    if (!scene->intersect(r, &info))
      return scene->background_color;

    /**
     * TODO: render lights
     */
    unsigned int index = -1;
    real_t tt =scene->intersect_lights(r, index);
    if (index != -1) {
      if (tt > 0 && tt < info.t) {
        return Color3::White();
      }
    }

    //Color3 color = scene->ambient_light*info.ambient;

    /**
     * TODO: compute ambient term here
     */
    Color3 color = Color3::Black();

    // send shooting rays for each light
    for (size_t i=0; i<scene->num_lights(); i++) {
      // the rays should be sent repeated
      for (size_t k=0; k<NUM_SHOOTING_RAY; k++) {
        // generate random direction ray and depth
        Vector3 dir_l = point_on_sphere();
        Ray ray_l = Ray(scene->get_lights()[i].position, dir_l);
        unsigned int dep_l = (unsigned int) (random()*MAX_BOUNCE + 1);
        // create intersection informaton
        Intersection info_l;
        // get color of the return ray
        Color3 temp = shooting_ray_trace(ray_l, dep_l, info_l);
        // check if there is blockage
        Ray ray_s = Ray(info_l.pos, info.pos-info_l.pos);
        Intersection info_s;
        if (scene->intersect(ray_s, &info_s)) {
          // check intersection validity
          if (info_s.t>0.0 && info_s.idx!=info.idx && info_s.refractive_index>0.0) {
            color += (1.0/(real_t)(NUM_SHOOTING_RAY)) * temp;
          }
        }
      }
    }
    // send gathering rays from the contact point
    for (size_t k=0; k<NUM_GATHERING_RAY; k++) {
      // form random direction ray and depth
      Vector3 dir_o = cos_weighted_hemi(info.normal);
      Ray ray_o = Ray(info.pos, dir_o);
      unsigned int dep_o = (unsigned int) (random()*MAX_BOUNCE + 1);
      // get the color
      Intersection info_o;
      color += (1.0/(real_t)(NUM_GATHERING_RAY)) * gathering_ray_trace(ray_o, dep_o, info_o);
    }
    // adjust color for ambient
    color = color * info.ambient;

    Color3 reflColor = Color3::Black();
    
    color += direct_illumination(info);

    if (info.refractive_index > 0) {
      Ray reflected;
      reflected.d = reflect(r.d, info.normal);
      reflected.e = info.pos + 0.0001*reflected.d;
      reflColor = trace(reflected, depth-1);
      color = info.specular*reflColor;
      Ray R = Ray(info.pos + r.d*0.00001, r.d);
      color = refraction(R, &info, color, depth);
    } else if (info.specular != Color3::Black()) {
      Ray reflected;
      reflected.d = reflect(r.d, info.normal);
      reflected.e = info.pos + 0.00001*reflected.d;
      reflColor = trace(reflected, depth-1);
      color += info.specular*reflColor;
    }

    return info.texColor*color;
}

Color3 Raytracer::shooting_ray_trace (const Ray& r,  unsigned int depth, Intersection& info)
{
  // if pass limit, return nothing
  if (depth <=0) {
    return Color3::Black();
  }
  // check hitting
  if (!scene->intersect(r, &info)) {
    return scene->background_color;
  }
  // get the diffuse color on the current pixel that hit
  Color3 color = direct_illumination(info);
  // calculate reflected ray
  Ray reflected = Ray(info.pos, reflect(r.d, info.normal));
  // return the color with recursive call
  return color + shooting_ray_trace(reflected, depth-1, info);
}

Color3 Raytracer::gathering_ray_trace (const Ray& r, unsigned int depth, Intersection& info)
{
  // if pass limit, return nothing
  if (depth <=0) {
    return Color3::Black();
  }
  // check hitting
  if (!scene->intersect(r, &info)) {
    return scene->background_color;
  }
  // render light
  unsigned int index = -1;
  real_t tt =scene->intersect_lights(r, index);
  if (index != -1) {
    if (tt > 0 && tt < info.t) {
      return Color3::White();
    }
  }
  // ready the color
  Color3 color = Color3::Black();
  Color3 reflColor = Color3::Black();
  color += direct_illumination(info);
  if (info.refractive_index > 0) {
    Ray reflected;
    reflected.d = reflect(r.d, info.normal);
    reflected.e = info.pos + 0.0001*reflected.d;
    reflColor = gathering_ray_trace(reflected, depth-1, info);
    color = info.specular*reflColor;
    Ray R = Ray(info.pos + r.d*0.00001, r.d);
    color = gathering_refraction(R, info, color, depth);
  } else if (info.specular != Color3::Black()) {
    Ray reflected;
    reflected.d = reflect(r.d, info.normal);
    reflected.e = info.pos + 0.00001*reflected.d;
    reflColor = gathering_ray_trace(reflected, depth-1, info);
    color += info.specular*reflColor;
  }
  return info.texColor*color;
}

real_t squared_length(const Color3& color)
{
    return color.r*color.r + color.g*color.g + color.b*color.b;
}

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
      real_t i = real_t(2)*(real_t(x)+random())*dx - real_t(1);
      real_t j = real_t(2)*(real_t(y)+random())*dy - real_t(1);

      Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));
      //res += trace(r, MAX_DEPTH);
      res += trace(r, MAX_DEPTH);
    }

    //std::cout << "Pixel " << x << " " << y << std::endl;

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
