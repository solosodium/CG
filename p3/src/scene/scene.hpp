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
//#include "scene/bvh.hpp"
#include "ray.hpp"
#include <string>
#include <vector>
#include <cfloat>

namespace _462 {

/************************************************************************/
/***************** Added new functions implementations ******************/
/************************************************************************/

// number of shadow rays sent to light source
static const size_t NUM_SHADOW_RAY = 20;

// minimum t threshold
static const real_t MIN_T_THRESHOLD = 0.000001;

// a very small number
static const real_t EPSILON = 0.000001;


/************************************************************************/

/*******************************************/
/***** Section for intersection struct *****/
/*******************************************/

struct Intersect
{
	// mameber parameters
	
    bool intersect;     		// boolean indicator of intersection
    Vector3 position;   		// position of the intersection
    Vector3 normal;     		// normal of the intersection
	//const Material* material;	// material of the intersection
	
    // added members

    Color3 ambient;             // ambient material color
    Color3 diffuse;             // diffuse material color
    Color3 specular;            // specular material color
    real_t refractive_index;    // refractive index
    Color3 texture;             // texture color of the intersection
    
	// constructor
    
	Intersect() 
	{
        // TODO: something to add?
    }
};

/******************************************/

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

    // Inverse transformation matrix
	Matrix4 invMat;
    // Normal transformation matrix
	Matrix3 normMat;

    /**
     * Renders this geometry using OpenGL in the local coordinate space.
     */
    virtual void render() const = 0;

	/**
	 * Added virtual function for intersection
	 */
	virtual Intersect intersect (Ray ray) const = 0;

	bool initialize();
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

	bool intersect(const Ray& r, real_t& t);

    // The position of the light, relative to world origin.
    Vector3 position;
    // The color of the light (both diffuse and specular)
    Color3 color;
    // attenuation
    Attenuation attenuation;
    // radius of the light
	real_t radius;
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

    // accessor functions
    Geometry* const* get_geometries() const;
    size_t num_geometries() const;
    const SphereLight* get_lights() const;
    size_t num_lights() const;
    Material* const* get_materials() const;
    size_t num_materials() const;
    Mesh* const* get_meshes() const;
    size_t num_meshes() const;

    /// Clears the scene, and invokes delete on everything in geometries.
    void reset();

    // functions to add things to the scene
    // all pointers are deleted by the scene upon scene deconstruction.
    void add_geometry( Geometry* g );
    void add_material( Material* m );
    void add_mesh( Mesh* m );
    void add_light( const SphereLight& l );

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
