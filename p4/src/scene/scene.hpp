/**
 * @file scene.hpp
 * @brief Class definitions for scenes.
 *
 */

#ifndef _462_SCENE_SCENE_HPP_
#define _462_SCENE_SCENE_HPP_

#include "math/vector.hpp"
#include "math/quaternion.hpp"
#include "math/matrix.hpp"
#include "math/camera.hpp"
#include "scene/material.hpp"
#include "scene/mesh.hpp"
#include "ray.hpp"
#include <string>
#include <vector>
#include <cfloat>

namespace _462 {

/**
 * TODO: world coordinates bounds data structure
 */
struct Bounds
{
	// member values
    real_t x_min;
    real_t x_max;

    real_t y_min;
    real_t y_max;
    
    real_t z_min;
    real_t z_max;

    real_t x_center;
    real_t y_center;
    real_t z_center;

	// constructor with default values
	Bounds () : x_min(0.0), x_max(0.0),
                y_min(0.0), y_max(0.0),
                z_min(0.0), z_max(0.0),
                x_center(0.0), y_center(0.0), z_center(0.0) {}
};

struct Intersection
{
	Vector3 pos;
	Vector3 normal;
	Vector3 bary; // triangle barycentric coords

	Color3 ambient, specular, diffuse, texColor;

	real_t t, refractive_index;

	unsigned int idx;
};

class Geometry
{
public:
    Geometry();
    virtual ~Geometry();

    /*
       World transformation are applied in the following order:
       1. Scale
       2. Orientation
       3. Position
    */

    // The world position of the object.
    Vector3 position;

    // The world orientation of the object.
    // Use Quaternion::to_matrix to get the rotation matrix.
    Quaternion orientation;

    // The world scale of the object.
    Vector3 scale;

    Matrix4 invMat;

	/**
	 * TODO: keep invMat_inv matrix here for fast computation
	 */
	Matrix4 invMat_inv;

    Matrix3 normMat;

    /**
     * Renders this geometry using OpenGL in the local coordinate space.
     */
    virtual void render() const = 0;

    virtual bool intersect(const Ray& r, Intersection* info)
    {
      return false;
    }

    bool initialize();

    virtual void get_extended_info(Intersection* info);

    /**
     * TODO: define get bounds function interface for Geometry
     */
    virtual void get_bounds (Bounds& bounds) const = 0;

};

struct SphereLight
{
    struct Attenuation
    {
        real_t constant;
        real_t linear;
        real_t quadratic;
    };

    SphereLight();

    bool initialize();

    bool intersect(const Ray& r, real_t& t) const;

    // The position of the light, relative to world origin.
    Vector3 position;
    Quaternion orientation;
    Vector3 scale;

    Matrix4 invMat, mat;
    Matrix3 normMat;

    // The color of the light (both diffuse and specular)
    Color3 color;
    // attenuation
    Attenuation attenuation;
    real_t radius;
};

/**
 * TODO: bounding box data structure (in world coordinates)
 */
class BoundingBox
{
public:
	std::vector<size_t> indices;	// list of geometries indices
	
    Bounds bounds;                  // bounds of the bounding box

	BoundingBox* left;				// left child pointer
	BoundingBox* right;				// right child pointer
	
	// constructor
	BoundingBox () {
		left = NULL;
		right = NULL;
	}
};

/**
 * The container class for information used to render a scene composed of
 * Geometries.
 */
class Scene
{
public:

    /// the camera
    Camera camera;
    /// the background color
    Color3 background_color;
    /// the amibient light of the scene
    Color3 ambient_light;
    /// the refraction index of air
    real_t refractive_index;

    /// Creates a new empty scene.
    Scene();

    /// Destroys this scene. Invokes delete on everything in geometries.
    ~Scene();

	bool initialize();

	bool intersect(const Ray& r, Intersection* info);

    // accessor functions
    Geometry* const* get_geometries() const;
    size_t num_geometries() const;
    const SphereLight* get_lights() const;
    size_t num_lights() const;
    Material* const* get_materials() const;
    size_t num_materials() const;
    Mesh* const* get_meshes() const;
    size_t num_meshes() const;

    real_t intersect_lights(const Ray& r, unsigned int& idx);

    /// Clears the scene, and invokes delete on everything in geometries.
    void reset();

    // functions to add things to the scene
    // all pointers are deleted by the scene upon scene deconstruction.
    void add_geometry( Geometry* g );
    void add_material( Material* m );
    void add_mesh( Mesh* m );
    void add_light( const SphereLight& l );

	/**
	 * TODO: declare the parameters and functions for creating spatial-partitioning data structure
	 */
    std::vector<Bounds> geometry_bounds;    // bounds of all geometries

    void initialize_bounds ();              // function to initialize bounds for all geometries

	BoundingBox* root;                      // root bounding box

    BoundingBox* build_bounding_boxes ( std::vector<size_t> list );    // function to build bounding boxes

    real_t compute_surface_area ( std::vector<size_t> list );   // compute the surface area of the list (half + half)
	
    void intersect_bounding_boxes ( BoundingBox* root, const Ray& r, std::vector<size_t>& target );

private:

    typedef std::vector< SphereLight > SphereLightList;
    typedef std::vector< Material* > MaterialList;
    typedef std::vector< Mesh* > MeshList;
    typedef std::vector< Geometry* > GeometryList;

    // list of all lights in the scene
    SphereLightList point_lights;
    // all materials used by geometries
    MaterialList materials;
    // all meshes used by models
    MeshList meshes;
    // list of all geometries. deleted in dctor, so should be allocated on heap.
    GeometryList geometries;

private:

    // no meaningful assignment or copy
    Scene(const Scene&);
    Scene& operator=(const Scene&);

};

} /* _462 */

#endif /* _462_SCENE_SCENE_HPP_ */
