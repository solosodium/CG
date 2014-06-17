/**
 * @file model.cpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "scene/model.hpp"
#include "scene/material.hpp"
#include "application/opengl.hpp"
#include "scene/triangle.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>


namespace _462 {

Model::Model() : mesh( 0 ), material( 0 ) { }
Model::~Model() { }

void Model::render() const
{
    if ( !mesh )
        return;
    if ( material )
        material->set_gl_state();
    mesh->render();
    if ( material )
        material->reset_gl_state();
}

////////////////////////////////////////////////////////
// Added Function: intersection test for ray and model//
////////////////////////////////////////////////////////

Intersect Model::intersect (Ray ray) const
{
	// construct return value
    real_t min_dis = 9999999.9;         // minimum distance to begin with
	Intersect result = Intersect();     // final intersect object
	result.intersect = false;           // the default value is false

    // after transformaton, the center of the object should be (0,0,0)
    Vector3 E = invMat.transform_point(ray.e);        // ray eye coordinate
    Vector3 D = invMat.transform_vector(ray.d);       // ray direction vector

    // iterate through all the mesh triangles
    for (size_t kk=0; kk<mesh->num_triangles(); kk++)
    {
        // get the three vertices of the triangle
        MeshVertex first = mesh->vertices[ mesh->triangles[kk].vertices[0] ];
        MeshVertex second = mesh->vertices[ mesh->triangles[kk].vertices[1] ];
        MeshVertex third = mesh->vertices[ mesh->triangles[kk].vertices[2] ];

        // calculate parameters for the intersection test
        real_t a = first.position.x - second.position.x;
        real_t b = first.position.y - second.position.y;
        real_t c = first.position.z - second.position.z;

        real_t d = first.position.x - third.position.x;
        real_t e = first.position.y - third.position.y;
        real_t f = first.position.z - third.position.z;

        real_t g = D.x;
        real_t h = D.y;
        real_t i = D.z;

        real_t j = first.position.x - E.x;
        real_t k = first.position.y - E.y;
        real_t l = first.position.z - E.z;

        // value cache
        real_t ei_hf = e*i - h*f;
        real_t gf_di = g*f - d*i;
        real_t dh_eg = d*h - e*g;

        real_t M = a*ei_hf + b*gf_di + c*dh_eg;

        // check t intersection
        real_t ak_jb = a*k - j*b;
        real_t jc_al = j*c - a*l;
        real_t bl_kc = b*l - k*c;
        real_t t = -( f*ak_jb + e*jc_al +d*bl_kc ) / M;

        if (t < MIN_T_THRESHOLD) {
            continue;
        }

        // check gamma intersection
        real_t gamma = ( i*ak_jb + h*jc_al + g*bl_kc ) / M;
        if (gamma < 0.0 || gamma > 1.0) {
            continue;
        }

        // check beta intersection
        real_t beta = ( j*ei_hf + k*gf_di + l*dh_eg ) / M;
        if (beta < 0.0 || beta > 1.0-gamma) {
            continue;
        }

        // calculate current position
        // get normal transformation matrix
        Matrix4 invMat_inv;
        make_transformation_matrix(&invMat_inv, position, orientation, scale);
        // calcualte world position
        Vector3 local_pos = E + t*D;
        Vector3 world_pos = invMat_inv.transform_point(local_pos);

        // need to check if this is smaller than the current min_dis
        if (distance(ray.e, world_pos) >= min_dis) {
            continue;
        }

        // reach here means intersection found, you are allowed to modify returned result

        // update distance first
        min_dis = distance(ray.e, world_pos);

        // set intersection found flag
        result.intersect = true;

        // set intersection position
        result.position = world_pos;

        // compute and set normal
        Vector3 local_norm = normalize(first.normal + second.normal + third.normal);
        result.normal = normalize(normMat * local_norm);

        // calculate texture color at this point
        if (material->texture_filename.empty()) {
            // if none of the material has texture, return white
            result.texture = Color3::White();
        }
        else 
        {
            // get texture coordinates first
            Vector2 temp_tex = first.tex_coord * (1.0-beta-gamma) + 
                               second.tex_coord * beta + 
                               third.tex_coord * gamma;
            // Bug fix: here need to solve the texture wrap problem
            temp_tex.x = temp_tex.x - floor(temp_tex.x);
            temp_tex.y = temp_tex.y - floor(temp_tex.y);
            // get the width and height
            int width, height;
            material->get_texture_size(&width, &height);
            // set the texture color
            result.texture = material->get_texture_pixel( (int)(temp_tex.x*(real_t)(width)), (int)(temp_tex.y*(real_t)(height)) ); 
        }
    }

    // this can be done outside the for loop (save time)

    // calculate ambient color
    result.ambient = material->ambient;
    // calculate diffuse color
    result.diffuse = material->diffuse;
    // calculate specular color
    result.specular = material->specular;
    // calulate refractive index
    result.refractive_index = material->refractive_index;

	// return the value
	return result;
}

} /* _462 */
