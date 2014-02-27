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
#include "p1/project.hpp"
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

    // list of all wave-emitting points.
    WavePointList wave_points;

    real_t current_time;

};

} /* _462 */

#endif /* _462_SCENE_HEIGHTMAP_HPP_ */
