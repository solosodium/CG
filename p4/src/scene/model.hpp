/**
 * @file model.hpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_MODEL_HPP_
#define _462_SCENE_MODEL_HPP_

#include "scene/scene.hpp"
#include "scene/mesh.hpp"

namespace _462 {

/**
 * A mesh of triangles.
 */
class Model : public Geometry
{
public:

    const Mesh* mesh;
    const Material* material;

    Model();
    virtual ~Model();

    virtual void render() const;

	virtual bool intersect(const Ray& r, Intersection* info);
	virtual void get_extended_info(Intersection* info);

	/**
	 * TODO: implement get bounds function
	 */
	virtual void get_bounds (Bounds& bounds) const;
};


} /* _462 */

#endif /* _462_SCENE_MODEL_HPP_ */

