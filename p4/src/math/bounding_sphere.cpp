/**
 * @file bounding_sphere.cpp
 * @brief Function definitions for bounding sphers.
 *
 * @author Kristin Siu (kasiu)
 */

#include "math/bounding_sphere.hpp"

namespace _462 {

BoundingSphere::BoundingSphere( Vector3 center, real_t radius ) : center(center), radius(radius) {}

bool BoundingSphere::intersect( const BoundingSphere& other ) const
{
    return squared_distance( center, other.center ) < (( radius + other.radius ) * ( radius + other.radius ) );
}

} /* _462 */

