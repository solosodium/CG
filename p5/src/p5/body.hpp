#ifndef _462_PHYSICS_BODY_HPP_
#define _462_PHYSICS_BODY_HPP_

#include "math/vector.hpp"
#include "math/quaternion.hpp"
#include <exception>
#include <iostream>

namespace _462 {

class Body
{
public:
    int id;
    int type;
    Vector3 position;
    Quaternion orientation;
    Vector3 velocity;
    Vector3 angular_velocity;

    virtual ~Body() { }
    virtual Vector3 step_position( real_t dt, real_t motion_damping ) = 0;
    virtual Vector3 step_orientation( real_t dt, real_t motion_damping ) = 0;
    virtual void apply_force( const Vector3& f, const Vector3& offset ) = 0;

	// function to compute a RK4
	// this version computes first order (velocity and angular velocity)
	// assume torque and force are consistant
	// ini: initial value
	// grad: gradient of the field
	// step: step
	Vector3 RK4_1 (Vector3 grad, real_t step)
	{
		Vector3 k1 = step * grad;
		Vector3 k2 = step * grad;
		Vector3 k3 = step * grad;
		Vector3 k4 = step * grad;
		// return the final incremental value
		return 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;
	}

	// function to compute a RK4
	// this version computes second order (position and rotation)
	// it calls RK4_1 during computation
	// ini: initial value
	// grad: gradient of the field
	// grad_d: derivative of gradient of the field
	// step: step
	Vector3 RK4_2 (Vector3 grad, Vector3 grad_d, real_t step)
	{
		Vector3 k1 = grad * step;
		Vector3 k2 = (grad + RK4_1(grad_d, step/2.0)) * step;
		Vector3 k3 = (grad + RK4_1(grad_d, step/2.0)) * step;
		Vector3 k4 = (grad + RK4_1(grad_d, step)) * step;
		// return the final incremental value
		return 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;
	}
};

}

#endif
