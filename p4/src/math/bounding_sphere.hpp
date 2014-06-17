/**
 * @file bounding_sphere.cpp
 * @brief Function definitons for a simple bounding sphere class.
 *
 * @author Kristin Siu (kasiu)
 */

#ifndef _462_MATH_BOUNDING_SPHERE_HPP_
#define _462_MATH_BOUNDING_SPHERE_HPP_

#include "math/vector.hpp"
#include "math/math.hpp"

namespace _462 {

class BoundingSphere {
public:
    BoundingSphere(){} // stupid default constructor
    BoundingSphere( Vector3 center, real_t radius );

    /**
     * @brief Returns the amount of intersection between the two spheres, if any
     * @param other Reference to the sphere being tested against
     */
    bool intersect( const BoundingSphere& other ) const;

    Vector3 center;
    real_t radius;
};

} /* _462 */

#endif /* _462_MATH_BOUNDING_SPHERE_HPP */
