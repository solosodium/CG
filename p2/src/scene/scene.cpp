/**
 * @file scene.cpp
 * @brief Function definitions for scenes.
 *
 */

#include "scene/scene.hpp"
#include <ctime>
namespace _462 {


Geometry::Geometry():
    position( Vector3::Zero ),
    orientation( Quaternion::Identity ),
    scale( Vector3::Ones )
{

}

Geometry::~Geometry() { }



PointLight::PointLight():
    position( Vector3::Zero ),
    color( Color3::White )
{
    attenuation.constant = 1;
    attenuation.linear = 0;
    attenuation.quadratic = 0;
}


Scene::Scene()
{
    reset();
}

Scene::~Scene()
{
    reset();
}

bool Scene::subdivide_geometries()
{
    bool success = true;
    int n = geometries.size();
    std::clock_t start;
    start = std::clock();
    for (int i = 0; i < n; i++)
    {
	success &= geometries[i]->subdivide();
    }
    double duration = (std::clock() - start)/double(CLOCKS_PER_SEC);
    std::cout << "Subdivision took " << duration << " seconds." << std::endl;
    return success;
}

const Scene::GeometryList& Scene::get_geometries() const
{
    return geometries;
}

const Scene::PointLightList& Scene::get_lights() const
{
    return point_lights;
}

const Scene::MaterialList& Scene::get_materials() const
{
    return materials;
}

const Scene::MeshList& Scene::get_meshes() const
{
    return meshes;
}

void Scene::reset()
{
    for ( GeometryList::iterator i = geometries.begin(); i != geometries.end(); ++i ) {
        delete *i;
    }
    for ( MaterialList::iterator i = materials.begin(); i != materials.end(); ++i ) {
        delete *i;
    }
    for ( MeshList::iterator i = meshes.begin(); i != meshes.end(); ++i ) {
        delete *i;
    }

    geometries.clear();
    materials.clear();
    meshes.clear();
    point_lights.clear();

    camera = Camera();

    background_color = Color3::Black;
    ambient_light = Color3::Black;
    refractive_index = 1.0;
}

void Scene::add_geometry( Geometry* g )
{
    geometries.push_back( g );
}

void Scene::add_material( Material* m )
{
    materials.push_back( m );
}

void Scene::add_mesh( Mesh* m )
{
    meshes.push_back( m );
}

void Scene::add_light( const PointLight& l )
{
    point_lights.push_back( l );
}


} /* _462 */
