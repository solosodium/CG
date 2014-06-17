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

bool Model::intersect(const Ray& r, Intersection* info)
{
	unsigned int N = mesh->num_triangles();
	const MeshVertex* verts = mesh->get_vertices();
	Intersection tmpinfo;
	info->t = DBL_MAX;
	info->idx = -1;
	for (unsigned int i = 0; i < N; i++)
	{
		const unsigned int* vIdxs = mesh->triangles[i].vertices;
		if (triangle_intersect(r, verts[vIdxs[0]].position,
          verts[vIdxs[1]].position, verts[vIdxs[2]].position, &tmpinfo))
		{
			if (tmpinfo.t < info->t)
			{
				*info = tmpinfo;
				info->idx = i;
			}
		}
	}
	return info->idx != -1;
}

void Model::get_extended_info(Intersection* info)
{
	const unsigned int* vIdxs = mesh->triangles[info->idx].vertices;

	info->ambient = material->ambient;
	info->diffuse = material->diffuse;
	info->specular = material->specular;
	info->refractive_index = material->refractive_index;

	const Vector2& tx0 = mesh->vertices[vIdxs[0]].tex_coord;
	const Vector2& tx1 = mesh->vertices[vIdxs[1]].tex_coord;
	const Vector2& tx2 = mesh->vertices[vIdxs[2]].tex_coord;

	info->texColor = get_tex_color(info->bary, *material, tx0, tx1, tx2);

	const Vector3& n0 = mesh->vertices[vIdxs[0]].normal;
	const Vector3& n1 = mesh->vertices[vIdxs[1]].normal;
	const Vector3& n2 = mesh->vertices[vIdxs[2]].normal;

	info->normal = lerp(info->bary, n0, n1, n2);
}

/**
 * TODO: implement get bounds function
 */
void Model::get_bounds (Bounds& bounds) const
{
    // loop through all the vertices to get the bounds
	for (size_t k; k<mesh->num_vertices(); k++) 
	{
		// transform to world coordinates
		Vector3 vert_world = invMat_inv.transform_point((mesh->get_vertices())[k].position);
		// fit into either one of the bounds
		if (vert_world.x < bounds.x_min) {
			bounds.x_min = vert_world.x;
		}
		else if (vert_world.x > bounds.x_max) {
			bounds.x_max = vert_world.x;
		}
		if (vert_world.y < bounds.y_min) {
			bounds.y_min = vert_world.y;
		}
		else if (vert_world.y > bounds.y_max) {
			bounds.y_max = vert_world.y;
		}
		if (vert_world.z < bounds.z_min) {
			bounds.z_min = vert_world.z;
		}
		else if (vert_world.z > bounds.z_max) {
			bounds.z_max = vert_world.z;
		}
	}
	// settle down the center values as well
    bounds.x_center = (bounds.x_min + bounds.x_max) / 2.0;
    bounds.y_center = (bounds.y_min + bounds.y_max) / 2.0;
    bounds.z_center = (bounds.z_min + bounds.z_max) / 2.0;
}

} /* _462 */
