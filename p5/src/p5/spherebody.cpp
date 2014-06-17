#include "p5/spherebody.hpp"
#include "math/vector.hpp"
#include "math/matrix.hpp"
#include "scene/sphere.hpp"
#include <iostream>
#include <exception>
#include <algorithm>

namespace _462 {

SphereBody::SphereBody( Sphere* geom )
{
    sphere = geom;
    position = sphere->position;
    radius = sphere->radius;
    orientation = sphere->orientation;
    mass = 0.0;
    velocity = Vector3::Zero();
    angular_velocity = Vector3::Zero();
    force = Vector3::Zero();
    torque = Vector3::Zero();
}

Vector3 SphereBody::step_position( real_t dt, real_t motion_damping )
{
    // Note: This function is here as a hint for an approach to take towards
    // programming RK4, you should add more functions to help you or change the
    // scheme
    // TODO return the delta in position dt in the future
	
	// calculate velocity by RK4
	/*
	Vector3 k1 = velocity * dt;
	Vector3 k2 = (velocity + force/mass * dt/2.0) * dt;
	Vector3 k3 = (velocity + force/mass * dt/2.0) * dt;
	Vector3 k4 = (velocity + force/mass * dt) *dt;
	Vector3 d_pos = 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;
	position += d_pos;
	velocity += force/mass * dt;
	*/
	Vector3 d_pos = RK4_2(velocity, force/mass, dt);
	position += d_pos;
	velocity += RK4_1(force/mass, dt);
	velocity *= 1.0-motion_damping;
	
	// update the actual geometry
	sphere->position = position;

    return d_pos;
}

Vector3 SphereBody::step_orientation( real_t dt, real_t motion_damping )
{
    // Note: This function is here as a hint for an approach to take towards
    // programming RK4, you should add more functions to help you or change the
    // scheme
    // TODO return the delta in orientation dt in the future
    // vec.x = rotation along x axis
    // vec.y = rotation along y axis
    // vec.z = rotation along z axis

	// calculate rotation with RK4
	real_t I = 0.4*mass*radius*radius;	// moment of inertia

	/*
	Vector3 k1 = angular_velocity * dt;
	Vector3 k2 = (angular_velocity + torque/I * dt/2.0) * dt;
	Vector3 k3 = (angular_velocity + torque/I * dt/2.0) * dt;
	Vector3 k4 = (angular_velocity + torque/I * dt) * dt;
	Vector3 d_rot = 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;
	*/
	Vector3 d_rot = RK4_2(angular_velocity, torque/I, dt);
	angular_velocity += RK4_1(torque/I, dt);
	angular_velocity *= 1.0-motion_damping;
	
	// update the actual geometry
	orientation = normalize( Quaternion( orientation*Vector3::UnitX(), d_rot.x ) * orientation );
	orientation = normalize( Quaternion( orientation*Vector3::UnitY(), d_rot.y ) * orientation );
	orientation = normalize( Quaternion( orientation*Vector3::UnitZ(), d_rot.z ) * orientation );
	sphere->orientation = orientation;
	
    return d_rot;
}

void SphereBody::apply_force( const Vector3& f, const Vector3& offset )
{
    // TODO apply force/torque to sphere

    // if offset is zero (less than epsilon), only linear component
    if (length(offset) < 1e-9) {
        // all force goes to linear
        force += f;
        // torque is zero
        torque += Vector3::Zero();
    }
    // if offset is not equal zero, both linear and angular components
    else {
        // linear term is the projection of force on the offset
        force += dot(f, offset) * offset / length(offset) / length(offset);
        // angular term is cross product between force and offset
		// NOTE: order matters
        torque += cross(offset, f);
    }
}

void SphereBody::clear_force()
{
	// clear force and torque
	force = Vector3::Zero();
	torque = Vector3::Zero();
}

}
