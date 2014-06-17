#include "p5/collisions.hpp"

namespace _462 {

bool collides( SphereBody& body1, SphereBody& body2, real_t collision_damping )
{
    // TODO detect collision. If there is one, update velocity

	// use velocity as the zeroth criteria
	// check veocity absolute value, compare with epsilon
	if ( length(body2.velocity-body1.velocity) < 1e-9 ) {
		// collision is not detected
		return false;
	}
	// position and velocity alignment
	if ( dot ( (body1.position-body2.position), (body2.velocity-body1.velocity) ) < 0.0) {
		// collision is not detected
		return false;
	}

	// first check if collison criteria meet
	if ( distance(body1.position, body2.position) < (body1.radius + body2.radius) ) {
		// collision detected, calculation follows project descripton
		Vector3 v1_p = body1.velocity - body2.velocity;
		Vector3 v2_p = Vector3::Zero();
		Vector3 d = (body2.position - body1.position)/length(body2.position - body1.position);
		Vector3 v2_pp = 2*d *(body1.mass/(body1.mass+body2.mass))*dot(v1_p, d);
		// set body2's velocity
		Vector3 v2 = body2.velocity;	// remember the old velocity
		body2.velocity += v2_pp;
		// set body1's velocity
		body1.velocity = (body1.mass*body1.velocity + body2.mass*v2 - body2.mass*body2.velocity)/body1.mass;
		// damp both velocities
		body1.velocity *= (1-collision_damping);
		body2.velocity *= (1-collision_damping);
		// adjust velocity if small
		if (length(body1.velocity) < 1e-9)
			body1.velocity = Vector3::Zero();
		if (length(body2.velocity) < 1e-9)
			body2.velocity = Vector3::Zero();
		// return result
		return true;
	}
	else {
		// collision is not detected
		return false;
	}
}

bool collides( SphereBody& body1, TriangleBody& body2, real_t collision_damping )
{
    // TODO detect collision. If there is one, update velocity
	// NOTE only body1 needs to update its velocity
	
	// first check if collision criteria meet
	// NOTE checking with two stages, first distance, second within triangle
	
	// pre-compute several parameters
	Vector3 a = body1.position - body2.vertices[0];
	// compute normal vector
	Vector3 normal = Vector3::Zero();
	normal = cross(body2.vertices[1]-body2.vertices[0], body2.vertices[2]-body2.vertices[0]);
	normal = normalize(normal);
	// adjust direction of the normal vector to match sign
	if (dot(a, normal) < 0.0) {
		normal = -normal;
	}
	
	// use velocity as the zeroth criteria
	// velocity absolute value
	if ( length(body1.velocity) < 1e-9 ) {
		// collision not detected
		return false;
	}
	// position and velocity alignment
	if ( dot ( normal, body1.velocity ) > 0.0) {
		// collision is not detected
		return false;
	}
	
	// check distance first
	if (length(dot(a, normal)*normal) >= body1.radius) {
		// collision is not detected
		return false;
	}
	
	// check point in the triangle
	Vector3 p = body2.vertices[0] + (a - dot(a, normal)*normal);
	// compute recalculated coordinates
	Vector3 v0 = p - body2.vertices[0];
	Vector3 v1 = body2.vertices[1] - body2.vertices[0];
	Vector3 v2 = body2.vertices[2] - body2.vertices[0];
	// compute helping parameters
	real_t A = dot(v0, v1);
	real_t B = dot(v0, v2);
	real_t C = dot(v1, v1);
	real_t D = dot(v1, v2);
	real_t E = dot(v2, v2);
	// check first parameter
	real_t alpha = (A*E - B*D) / (C*E - D*D);
	if (alpha < 0.0 || alpha > 1.0) {
		// collision is not detected
		return false;
	}
	// check second parameter
	real_t beta = (A*D - B*C) / (D*D - C*E);
	if (beta < 0.0 || beta > 1.0) {
		// collision is not detected
		return false;
	}
	// check two parameters together
	if (alpha + beta > 1.0) {
		// collision is not detected
		return false;
	}
	
	// if reach here, it means collision is detected, calculate velocity
	body1.velocity -= 2*dot(body1.velocity, normal)*normal;
	// damp velocity
	body1.velocity *= (1-collision_damping);
	// adjust velocity if small
	if (length(body1.velocity) < 1e-9)
		body1.velocity = Vector3::Zero();
	// return result
    return true;
}

bool collides( SphereBody& body1, PlaneBody& body2, real_t collision_damping )
{
    // TODO detect collision. If there is one, update velocity
	// NOTE only body1 needs to update its velocity
	
	// pre compute some values
	Vector3 a = body1.position - body2.position;
	Vector3 normal = body2.normal;
	// find the normal align with the position direction
	if (dot(a, normal) < 0.0) {
		normal = -normal;
	}
	real_t d = dot(a, normal);
	
	// use velocity as the zeroth criteria
	// check absolute velocity value
	if ( length(body1.velocity) < 1e-9) {
		// collision not detected
		return false;
	}
	// position and velocity alignment
	if ( dot ( normal, body1.velocity ) > 0.0 ) {
		// collision is not detected
		return false;
	}
	
	// first check if collision criteria meet
	if (d < body1.radius) {
		// collision detected, calculate new velocity
		body1.velocity -= 2*dot(body1.velocity, normal)*normal;
		// damp velocity
		body1.velocity *= (1-collision_damping);
		// adjust velocity if small
		if (length(body1.velocity) < 1e-9)
			body1.velocity = Vector3::Zero();
		return true;
	}
	else {
		// collision is not detected
		return false;
	}
}

}
