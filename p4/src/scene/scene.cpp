/**
 * @file scene.cpp
 * @brief Function definitions for scenes.
 *
 * @author Eric Butler (edbutler)
 * @author Kristin Siu (kasiu)
 */

#include "scene/scene.hpp"

namespace _462 {


Geometry::Geometry():
    position(Vector3::Zero()),
    orientation(Quaternion::Identity()),
    scale(Vector3::Ones())
{

}

Geometry::~Geometry() { }

bool Geometry::initialize()
{
    make_inverse_transformation_matrix(&invMat, position, orientation, scale);
    //Matrix4 mat;
    //make_transformation_matrix(&mat, position, orientation, scale);
    //make_normal_matrix(&normMat, mat);
	/**
	 * TODO: initialize invMat_inv
	 */
	make_transformation_matrix(&invMat_inv, position, orientation, scale);
	make_normal_matrix(&normMat, invMat_inv);
	
    return true;
}

bool SphereLight::initialize()
{

    //std::cout << position.x << std::endl;
    make_inverse_transformation_matrix(&invMat, position, orientation, scale);
    //Matrix4 mat;
    make_transformation_matrix(&mat, position, orientation, scale);
    make_normal_matrix(&normMat, mat);
    return true;
}

SphereLight::SphereLight():
    position(Vector3::Zero()),
    color(Color3::White()),
    radius(real_t(0)),
    orientation(Quaternion::Identity()),
    scale(Vector3::Ones())
{
    attenuation.constant = 1;
    attenuation.linear = 0;
    attenuation.quadratic = 0;

}

bool solveQuadratic(const real_t &a, const real_t &b, const real_t &c, real_t &x0, real_t &x1) 
{
	real_t discr = b * b - 4 * a * c;
	if (discr < 0)
		return false;
	else if (discr == 0)
		x0 = x1 = - 0.5 * b / a;
	else {
		real_t q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
		x0 = q / a; x1 = c / q;
	}
	if (x0 > x1)
		std::swap(x0, x1);
	return true;
}

bool SphereLight::intersect(const Ray& tr, real_t& t) const
{
    Ray r;
    r.d = invMat.transform_vector(tr.d);
    r.e = invMat.transform_point(tr.e);

    //std::cout << r.e.x << " " << tr.e.x << std::endl;

    real_t a = squared_length(r.d);
    real_t b = real_t(2)*dot(r.d, r.e);
    real_t c = squared_length(r.e) - radius * radius;
    real_t disc = b*b - real_t(4)*a*c;

    if (disc < real_t(0)) return false;

    int df = b < 0 ? -1 : 1;
    real_t q = (-b + df*sqrt(disc))/real_t(2);

    real_t t0 = q/a;
    real_t t1 = c/q;

    if (t0 > t1) std::swap(t0, t1);

    if (t1 < real_t(0)) return false;

    t = t0 > real_t(0) ? t0 : t1;
    //info->normal = r.d*(info->t+1e-9) + r.e;
    //info->normal = normalize(normMat*info->normal);

    return true;
}

void Geometry::get_extended_info(Intersection* info)
{}

Scene::Scene()
{
    reset();
}

Scene::~Scene()
{
    reset();
}

bool Scene::initialize()
{
	bool res = true;
	for (unsigned int i = 0; i < num_geometries(); i++)
	    res &= geometries[i]->initialize();
	for (unsigned int i = 0; i < num_lights(); i++)
	    res &= point_lights[i].initialize();

    /**
     * TODO: additional initialization for bounds for geometries
     */
    initialize_bounds();        // initialize bounds for all geometries

    std::vector<size_t> list;
    for (size_t k=0; k<num_geometries(); k++) {
        list.push_back(k);
    }
    root = build_bounding_boxes( list );

	return res;
}

bool Scene::intersect(const Ray& r, Intersection* info)
{
	// TODO: this should be changed for BVH

    std::vector<size_t> target;
    intersect_bounding_boxes ( root, r, target );
	
    info->t = DBL_MAX;
    Intersection tmpinfo;
    Ray tr(r);
    int objIdx = -1;

    for (unsigned int i = 0; i < target.size(); i++)
    {
        tr.d = (geometries[target[i]]->invMat.transform_vector(r.d));
        tr.e = geometries[target[i]]->invMat.transform_point(r.e);

        if (geometries[target[i]]->intersect(tr, &tmpinfo) && tmpinfo.t < info->t)
        {
            *info = tmpinfo;
            objIdx = target[i];
        }
    }
    
    /*
    for (unsigned int i = 0; i < num_geometries(); i++)
    {
		tr.d = (geometries[i]->invMat.transform_vector(r.d));
		tr.e = geometries[i]->invMat.transform_point(r.e);

		if (geometries[i]->intersect(tr, &tmpinfo) && tmpinfo.t < info->t)
		{
	    	*info = tmpinfo;
	    	objIdx = i;
		}
    }
    */

    if (objIdx != -1)
    {
		info->pos = r.d*(info->t+1e-9) + r.e;
		geometries[objIdx]->get_extended_info(info);
		info->normal = normalize(geometries[objIdx]->normMat*info->normal);
		return true;
    }

    return false;
}

real_t Scene::intersect_lights(const Ray& r, unsigned int& idx)
{
    unsigned int N = num_lights();
    const SphereLight* lights = get_lights();

    real_t mint = DBL_MAX;

    for (unsigned int i = 0; i < N; i++)
    {
	real_t tmpt = DBL_MAX;
	if (lights[i].intersect(r, tmpt) && tmpt < mint)
	{
	    mint = tmpt;
	    idx = i;
	}
    }

    return mint;
}

Geometry* const* Scene::get_geometries() const
{
    return geometries.empty() ? NULL : &geometries[0];
}

size_t Scene::num_geometries() const
{
    return geometries.size();
}

const SphereLight* Scene::get_lights() const
{
    return point_lights.empty() ? NULL : &point_lights[0];
}

size_t Scene::num_lights() const
{
    return point_lights.size();
}

Material* const* Scene::get_materials() const
{
    return materials.empty() ? NULL : &materials[0];
}

size_t Scene::num_materials() const
{
    return materials.size();
}

Mesh* const* Scene::get_meshes() const
{
    return meshes.empty() ? NULL : &meshes[0];
}

size_t Scene::num_meshes() const
{
    return meshes.size();
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

    background_color = Color3::Black();
    ambient_light = Color3::Black();
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

void Scene::add_light( const SphereLight& l )
{
    point_lights.push_back( l );
}

/**
 * TODO: implement functions
 */
void Scene::initialize_bounds ()
{
    // loop through all geometries to create bounds
	for (size_t k=0; k<num_geometries(); k++) {
        Bounds bounds_temp;                     // placeholder
        geometries[k]->get_bounds(bounds_temp); // get bounds
        geometry_bounds.push_back(bounds_temp); // add to list

        //std::cout << bounds_temp.x_min << " " << bounds_temp.x_max << " : " <<
        //             bounds_temp.y_min << " " << bounds_temp.y_max << " : " <<
        //             bounds_temp.z_min << " " << bounds_temp.z_max << std::endl;
    }

    //std::cout << "Geos: " << num_geometries() << std::endl;
}

BoundingBox* Scene::build_bounding_boxes ( std::vector<size_t> list )
{
    // base case 0, list empty
    if (list.empty()) {
        BoundingBox* node = NULL;
        return node;
    }
    // base case 1, list has one
    else if (list.size() == 1) {
        BoundingBox* node = new BoundingBox();
        node->indices = list;
        node->bounds = geometry_bounds[list[0]];
        node->left = NULL;
        node->right = NULL;

        //std::cout << "Node" << std::endl;
        
        return node;
    }
    // general case
    else {
        BoundingBox* node = new BoundingBox();
        node->indices = list;
        // compute bounds
        for (size_t k=0; k<list.size(); k++) {
            // fit into either one of the bounds
            if (geometry_bounds[list[k]].x_min < node->bounds.x_min) {
                node->bounds.x_min = geometry_bounds[list[k]].x_min;
            }
            if (geometry_bounds[list[k]].x_max > node->bounds.x_max) {
                node->bounds.x_max = geometry_bounds[list[k]].x_max;
            }
            if (geometry_bounds[list[k]].y_min < node->bounds.y_min) {
                node->bounds.y_min = geometry_bounds[list[k]].y_min;
            }
            if (geometry_bounds[list[k]].y_max > node->bounds.y_max) {
                node->bounds.y_max = geometry_bounds[list[k]].y_max;
            }
            if (geometry_bounds[list[k]].z_min < node->bounds.z_min) {
                node->bounds.z_min = geometry_bounds[list[k]].z_min;
            }
            if (geometry_bounds[list[k]].z_max > node->bounds.z_max) {
                node->bounds.z_max = geometry_bounds[list[k]].z_max;
            }
        }
        // decide how to split the list
        // it could be by x, y or z
        std::vector<size_t> by_x = list;
        std::vector<size_t> by_y = list;
        std::vector<size_t> by_z = list;
        // bubble sort by x_center
        for (size_t i=0; i<list.size(); i++) {
            for (size_t k=i; k<list.size(); k++) {
                if (geometry_bounds[by_x[i]].x_center > geometry_bounds[by_x[k]].x_center) {
                    std::swap(by_x[i], by_x[k]);
                }
            }
        }
        // bubble sort by y_center
        for (size_t i=0; i<list.size(); i++) {
            for (size_t k=i; k<list.size(); k++) {
                if (geometry_bounds[by_y[i]].y_center > geometry_bounds[by_y[k]].y_center) {
                    std::swap(by_y[i], by_y[k]);
                }
            }
        }
        // bubble sort by z_center
        for (size_t i=0; i<list.size(); i++) {
            for (size_t k=i; k<list.size(); k++) {
                if (geometry_bounds[by_z[i]].z_center > geometry_bounds[by_z[k]].z_center) {
                    std::swap(by_z[i], by_z[k]);
                }
            }
        }
        // check which has the smallest
        real_t by_x_area = compute_surface_area(by_x);
        real_t by_y_area = compute_surface_area(by_y);
        real_t by_z_area = compute_surface_area(by_z);

        std::vector<size_t> left, right;

        if (by_x_area < by_y_area) {
            if (by_x_area < by_z_area) {
                // x smallest
                for (size_t k=0; k<by_x.size(); k++) {
                    if (k < by_x.size()/2)
                        left.push_back(by_x[k]);
                    else
                        right.push_back(by_x[k]);
                }
            }
            else {
                // z smallest
                for (size_t k=0; k<by_z.size(); k++) {
                    if (k < by_z.size()/2)
                        left.push_back(by_z[k]);
                    else
                        right.push_back(by_z[k]);
                }
            }
        }
        else {
            if (by_y_area < by_z_area) {
                // y smallest
                for (size_t k=0; k<by_y.size(); k++) {
                    if (k < by_y.size()/2)
                        left.push_back(by_y[k]);
                    else
                        right.push_back(by_y[k]);
                }
            }
            else {
                // z smallest
                for (size_t k=0; k<by_z.size(); k++) {
                    if (k < by_z.size()/2)
                        left.push_back(by_z[k]);
                    else
                        right.push_back(by_z[k]);
                }
            }
        }
        // set left and right and return
        node->left = build_bounding_boxes(left);
        node->right = build_bounding_boxes(right);

        //std::cout << "Node" << std::endl;

        return node;
    }
}

real_t Scene::compute_surface_area ( std::vector<size_t> list )
{
    Bounds lower_bounds, upper_bounds;      // the lower list and upper list
    // loop through all the elements
    for (size_t k=0; k<list.size(); k++) {
        // for lower half
        if (k < list.size()/2) {
            if (geometry_bounds[list[k]].x_min < lower_bounds.x_min) {
                lower_bounds.x_min = geometry_bounds[list[k]].x_min;
            }
            if (geometry_bounds[list[k]].x_max > lower_bounds.x_max) {
                lower_bounds.x_max = geometry_bounds[list[k]].x_max;
            }
            if (geometry_bounds[list[k]].y_min < lower_bounds.y_min) {
                lower_bounds.y_min = geometry_bounds[list[k]].y_min;
            }
            if (geometry_bounds[list[k]].y_max > lower_bounds.y_max) {
                lower_bounds.y_max = geometry_bounds[list[k]].y_max;
            }
            if (geometry_bounds[list[k]].z_min < lower_bounds.z_min) {
                lower_bounds.z_min = geometry_bounds[list[k]].z_min;
            }
            if (geometry_bounds[list[k]].z_max > lower_bounds.z_max) {
                lower_bounds.z_max = geometry_bounds[list[k]].z_max;
            }
        }
        // for upper half
        else {
            if (geometry_bounds[list[k]].x_min < upper_bounds.x_min) {
                upper_bounds.x_min = geometry_bounds[list[k]].x_min;
            }
            if (geometry_bounds[list[k]].x_max > upper_bounds.x_max) {
                upper_bounds.x_max = geometry_bounds[list[k]].x_max;
            }
            if (geometry_bounds[list[k]].y_min < upper_bounds.y_min) {
                upper_bounds.y_min = geometry_bounds[list[k]].y_min;
            }
            if (geometry_bounds[list[k]].y_max > upper_bounds.y_max) {
                upper_bounds.y_max = geometry_bounds[list[k]].y_max;
            }
            if (geometry_bounds[list[k]].z_min < upper_bounds.z_min) {
                upper_bounds.z_min = geometry_bounds[list[k]].z_min;
            }
            if (geometry_bounds[list[k]].z_max > upper_bounds.z_max) {
                upper_bounds.z_max = geometry_bounds[list[k]].z_max;
            }
        }
    }
    // calculate the surface area
    real_t lower_area = (lower_bounds.x_max-lower_bounds.x_min) *
                        (lower_bounds.y_max-lower_bounds.y_min) *
                        (lower_bounds.z_max-lower_bounds.z_min);
    real_t upper_area = (upper_bounds.x_max-upper_bounds.x_min) *
                        (upper_bounds.y_max-upper_bounds.y_min) *
                        (upper_bounds.z_max-upper_bounds.z_min);
    // return the value
    return lower_area + upper_area;
}

void Scene::intersect_bounding_boxes ( BoundingBox* root, const Ray& r, std::vector<size_t>& target )
{
    // check intersect
    size_t count = 0;
    real_t t, x_t, y_t, z_t;
    // check x_min
    t = (root->bounds.x_min - r.e.x) / r.d.x;
    y_t = r.e.y + t * r.d.y;
    z_t = r.e.z + t * r.d.z;
    if (t>0 && y_t < root->bounds.y_max && y_t > root->bounds.y_min && z_t < root->bounds.z_max && z_t > root->bounds.z_min)
        count ++;
    // check x_max
    t = (root->bounds.x_max - r.e.x) / r.d.x;
    y_t = r.e.y + t * r.d.y;
    z_t = r.e.z + t * r.d.z;
    if (t>0 && y_t < root->bounds.y_max && y_t > root->bounds.y_min && z_t < root->bounds.z_max && z_t > root->bounds.z_min)
        count ++;
    // check y_min
    t = (root->bounds.y_min - r.e.y) / r.d.y;
    x_t = r.e.x + t * r.d.x;
    z_t = r.e.z + t * r.d.z;
    if (t>0 && x_t < root->bounds.x_max && x_t > root->bounds.x_min && z_t < root->bounds.z_max && z_t > root->bounds.z_min)
        count ++;
    // check y_max
    t = (root->bounds.y_max - r.e.y) / r.d.y;
    x_t = r.e.x + t * r.d.x;
    z_t = r.e.z + t * r.d.z;
    if (t>0 && x_t < root->bounds.x_max && x_t > root->bounds.x_min && z_t < root->bounds.z_max && z_t > root->bounds.z_min)
        count ++;
    // check z_min
    t = (root->bounds.z_min - r.e.z) / r.d.z;
    x_t = r.e.x + t * r.d.x;
    y_t = r.e.y + t * r.d.y;
    if (t>0 && x_t < root->bounds.x_max && x_t > root->bounds.x_min && y_t < root->bounds.y_max && y_t > root->bounds.y_min)
        count ++;
    // check z_max
    t = (root->bounds.z_max - r.e.z) / r.d.z;
    x_t = r.e.x + t * r.d.x;
    y_t = r.e.y + t * r.d.y;
    if (t>0 && x_t < root->bounds.x_max && x_t > root->bounds.x_min && y_t < root->bounds.y_max && y_t > root->bounds.y_min)
        count ++;

    // return the value
    if (count > 0 && root->left == NULL && root->right == NULL) {
        target.push_back(root->indices[0]);
    }
    else if (count > 0 && root->left != NULL && root->right != NULL) {
        intersect_bounding_boxes(root->left, r, target);
        intersect_bounding_boxes(root->right, r, target);
    }
}


} /* _462 */
