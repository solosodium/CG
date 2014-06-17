/**
 * @file triangle.hpp
 * @brief Class definition for Triangle.
 *
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_TRIANGLE_HPP_
#define _462_SCENE_TRIANGLE_HPP_

#include "scene/scene.hpp"

namespace _462 {

/**
 * a triangle geometry.
 * Triangles consist of 3 vertices. Each vertex has its own position, normal,
 * texture coordinate, and material. These should all be interpolated to get
 * values in the middle of the triangle.
 * These values are all in local space, so it must still be transformed by
 * the Geometry's position, orientation, and scale.
 */
class Triangle : public Geometry
{
public:

    struct Vertex
    {
        // note that position and normal are in local space
        Vector3 position;
        Vector3 normal;
        Vector2 tex_coord;
        const Material* material;
    };

    // the triangle's vertices, in CCW order
    Vertex vertices[3];

    Triangle();
    virtual ~Triangle();
    virtual void render() const;

    virtual bool intersect(const Ray& r, Intersection* info);
    virtual void get_extended_info(Intersection* info);

    /**
     * TODO: implement get bounds function
     */
    virtual void get_bounds (Bounds& bounds) const;
};

template<class T>

inline T lerp(const Vector3& bary, const T& a, const T& b, const T& c)
{
	return bary.x*a + bary.y*b + bary.z*c;
}

bool triangle_intersect(const Ray& r, const Vector3& v0,
			const Vector3& v1, const Vector3& v2,
			Intersection* info);

Color3 get_tex_color(const Vector3& bary, const Material& mat,
		     const Vector2& v0, const Vector2& v1,
		     const Vector2& v2);

} /* _462 */

#endif /* _462_SCENE_TRIANGLE_HPP_ */
