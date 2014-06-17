#include "p5/physics.hpp"

namespace _462 {

Physics::Physics()
{
    reset();
}

Physics::~Physics()
{
    reset();
}

void Physics::step( real_t dt )
{
    // TODO step the world forward by dt. Need to detect collisions, apply
    // forces, and integrate positions and orientations.
    //
    // Note: put RK4 here, not in any of the physics bodies
    //
    // Must use the functions that you implemented
    //
    // Note, when you change the position/orientation of a physics object,
    // change the position/orientation of the graphical object that represents
    // it
	
	// set the number of substeps in this scope
	size_t num_of_substeps = 10;
	real_t t_sub = dt/ (real_t)(num_of_substeps);
	
	for (size_t i=0; i<num_of_substeps; i++)
	{
	
	// excert force by springs
	for (size_t n=0; n<num_springs(); n++) {
		springs[n]->step(t_sub);
	}
	
	// loop through all spheres
	for (size_t k=0; k<num_spheres(); k++) 
	{
		// excert gravity
		spheres[k]->apply_force(gravity, Vector3::Zero());
		
		// update position and orientation
		spheres[k]->step_position(t_sub, 0.0);
		spheres[k]->step_orientation(t_sub, 0.0);
		
		// check collision with spherebody
		for (size_t n=0; n<num_spheres(); n++) {
			if (k != n) {
				if ( collides(*spheres[k], *spheres[n], collision_damping) ) {
					// TODO: std::cout << "Sphere - Sphere" << std::endl;
				}
			}
		}
		
		// check collision with planebody
		for (size_t n=0; n<num_planes(); n++) {
			if ( collides(*spheres[k], *planes[n], collision_damping) ) {
				// TODO: std::cout << "Sphere - Plane" << std::endl;
			}
		}
		
		// check collision with trianglebody
		for (size_t n=0; n<num_triangles(); n++) {
			if ( collides(*spheres[k], *triangles[n], collision_damping) ) {
				// TODO: std::cout << "Sphere - Triangle" << std::endl;
			}
		}
		
		// clear all the forces and torque
		spheres[k]->clear_force();
	}
	
	}
}

void Physics::add_sphere( SphereBody* b )
{
    spheres.push_back( b );
}

size_t Physics::num_spheres() const
{
    return spheres.size();
}

void Physics::add_plane( PlaneBody* p )
{
    planes.push_back( p );
}

size_t Physics::num_planes() const
{
    return planes.size();
}

void Physics::add_triangle( TriangleBody* t )
{
    triangles.push_back( t );
}

size_t Physics::num_triangles() const
{
    return triangles.size();
}

void Physics::add_spring( Spring* s )
{
    springs.push_back( s );
}

size_t Physics::num_springs() const
{
    return springs.size();
}

void Physics::reset()
{
    for ( SphereList::iterator i = spheres.begin(); i != spheres.end(); i++ ) {
        delete *i;
    }
    for ( PlaneList::iterator i = planes.begin(); i != planes.end(); i++ ) {
        delete *i;
    }
    for ( TriangleList::iterator i = triangles.begin(); i != triangles.end(); i++ ) {
        delete *i;
    }
    for ( SpringList::iterator i = springs.begin(); i != springs.end(); i++ ) {
        delete *i;
    }

    spheres.clear();
    planes.clear();
    triangles.clear();
    springs.clear();
    
    gravity = Vector3::Zero();
	collision_damping = 0.0;
}

}
