/**
 * @file heightmap.hpp
 * @brief Backend for the heightmap.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_HEIGHTMAP_HPP_
#define _462_SCENE_HEIGHTMAP_HPP_

#include "math/vector.hpp"
#include "opengl/project.hpp"
#include <vector>

namespace _462 {

class WaterSurface : public Heightmap
{
public:
    /**
     * structure containing information about a wave-emitting point.
     */
    struct WavePoint {
        Vector2 position; // position on surface (between (-1 and 1)
        real_t falloff; // exponential falloff of amplitude
        real_t amplitude; // scaling factor of amplitude
        real_t timerate; // scaling factor of time
        real_t period; // scaling factor of distance
    };

    typedef std::vector<WavePoint> WavePointList;

    /**
     * Construct a new watersurface.
     */
    WaterSurface();

    virtual ~WaterSurface();

    /**
     * Returns the absolute height of the watersurface (in the local
     * coordinate space) for the given (x,z) and time.
     * @param pos The x and z positions. Valid range is (-1,-1) to (1,1).
     * @param time The absolute time.
     * @return The value of y for x, z, and time in the local coordinate space.
     */
    virtual real_t compute_height( const Vector2& pos ) const;

    virtual void update( real_t dt );
#ifdef SOL_OGL
    virtual bool create_gl_data( int resx, int resy );
    virtual void render() const;
#endif /* SOL_OGL */

    // list of all wave-emitting points.
    WavePointList wave_points;

    real_t current_time;

#ifdef SOL_OGL
private:
    typedef std::vector< float > FloatList;
    typedef std::vector< unsigned int > IndexList;
    typedef std::vector< Vector3 > VectorList;

    // the vertex data used for GL rendering
    FloatList vertex_data;
    // the index data used for GL rendering
    IndexList index_data;
    // the positions
    VectorList position_data;
    // the normals
    VectorList normal_data;

    int resx, resy;
#endif /* SOL_OGL */
};

} /* _462 */

#endif /* _462_SCENE_HEIGHTMAP_HPP_ */

