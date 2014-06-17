/**
 * @file heightmap.cpp
 * @brief Backend for the heightmap.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 * @author Frank Palermo (fpalermo)
 */


#include "scene/heightmap.hpp"

#ifdef SOL_OGL
#include "application/opengl.hpp"
#include <cstring>
#endif /* SOL_OGL */

namespace _462 {

WaterSurface::WaterSurface() : current_time( 0 ) { }

WaterSurface::~WaterSurface() { }

real_t WaterSurface::compute_height( const Vector2& pos ) const
{
    assert( pos.x >= -1.0 && pos.x <= 1.0 );
    assert( pos.y >= -1.0 && pos.y <= 1.0 );

    real_t h = 0;

    for ( size_t i = 0; i < wave_points.size(); ++i ) {
        const WavePoint& p = wave_points[i];
        real_t r = distance( pos, p.position );
        h += p.amplitude * exp( -p.falloff * r ) * sin( p.period * r + p.timerate * current_time );
    }

    return h;
}

#ifdef SOL_OGL
// number of floats per vertex
#define VERTEX_SIZE 8
#endif /* SOL_OGL */

void WaterSurface::update( real_t dt )
{
    current_time += dt;
#ifdef SOL_OGL
    size_t num_vertices = ( resx + 1 ) * ( resy + 1 );
    size_t num_triangles = resx * resy * 2;

    // update height, zero out all normals
    for ( size_t i = 0; i < num_vertices; ++i ) {
        position_data[i].y = compute_height( Vector2( position_data[i].x, position_data[i].z ) );
        normal_data[i] = Vector3::Zero;
    }

    // then sum in all triangle normals
    for ( size_t i = 0; i < num_triangles; ++i ) {
        unsigned int idx[3];
        for ( size_t j = 0; j < 3; ++j ) {
            idx[j] = index_data[i * 3 + j];
        }
        Vector3 pos[3];
        for ( size_t j = 0; j < 3; ++j ) {
            pos[j] = position_data[idx[j]];
        }
        Vector3 normal = normalize( cross( pos[1] - pos[0], pos[2] - pos[0] ) );
        assert( length( cross( pos[1] - pos[0], pos[2] - pos[0] ) ) > 0.0 );
        for ( size_t j = 0; j < 3; ++j ) {
            normal_data[idx[j]] += normal;
        }
    }

    for ( size_t i = 0; i < num_vertices; ++i ) {
        float* vertex = &vertex_data[VERTEX_SIZE * i];
        normalize( normal_data[i] ).to_array( &vertex[2] );
        vertex[6] = position_data[i].y;
    }
#endif /* SOL_OGL */
}

#ifdef SOL_OGL
bool WaterSurface::create_gl_data( int resx, int resy )
{
    if ( resx < 2 || resy < 2 ) {
        return false;
    }

    this->resx = resx;
    this->resy = resy;

    size_t num_vertices = ( resx + 1 ) * ( resy + 1 );
    size_t num_quads = resx * resy;

    // allocate vertex and index data
    vertex_data.resize( num_vertices * VERTEX_SIZE );
    index_data.resize( num_quads * 6 );
    position_data.resize( num_vertices );
    normal_data.resize( num_vertices );

    // build index data
    unsigned int* index = &index_data[0];
    for ( int i = 0; i < resx; ++i ) {
        for ( int j = 0; j < resy; ++j ) {
            unsigned int a = ( i + 0 ) * ( resy + 1 ) + ( j + 0 );
            unsigned int b = ( i + 0 ) * ( resy + 1 ) + ( j + 1 );
            unsigned int c = ( i + 1 ) * ( resy + 1 ) + ( j + 1 );
            unsigned int d = ( i + 1 ) * ( resy + 1 ) + ( j + 0 );

            index[0] = a;
            index[1] = b;
            index[2] = c;
            index[3] = c;
            index[4] = d;
            index[5] = a;
            index += 6;
        }
    }
    assert( index - 1 == &index_data.back() );

    float* vertex = &vertex_data[0];
    size_t idx = 0;
    for ( int i = 0; i <= resx; ++i ) {
        for ( int j = 0; j <= resy; ++j ) {
            // texture coord data
            vertex[0] = float( i ) / resx;
            vertex[1] = float( j ) / resy;
            // normal data
            normal_data[idx] = Vector3::UnitY;
            normal_data[idx].to_array( &vertex[2] );
            // vertex data
            position_data[idx] = Vector3(
                real_t( i * 2 ) / real_t( resx ) - 1.0,
                0.0,
                real_t( j * 2 ) / real_t( resy ) - 1.0
            );
            position_data[idx].to_array( &vertex[5] );
            vertex += VERTEX_SIZE;
            idx++;
        }
    }
    assert( vertex - 1 == &vertex_data.back() );

    return true;
}

void WaterSurface::render() const
{
    assert( index_data.size() > 0 );
    glInterleavedArrays( GL_T2F_N3F_V3F, VERTEX_SIZE * sizeof vertex_data[0], &vertex_data[0] );
    glDrawElements( GL_TRIANGLES, index_data.size(), GL_UNSIGNED_INT, &index_data[0] );
}
#endif /* SOL_OGL */

} /* _462 */

