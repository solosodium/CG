/**
 * @file triangle.cpp
 * @brief Function definitions for the Triangle class.
 *
 * @author Eric Butler (edbutler)
 */

#include "scene/triangle.hpp"
#include "application/opengl.hpp"

namespace _462 {

Triangle::Triangle()
{
    vertices[0].material = 0;
    vertices[1].material = 0;
    vertices[2].material = 0;
}

Triangle::~Triangle() { }

void Triangle::render() const
{
    bool materials_nonnull = true;
    for ( int i = 0; i < 3; ++i )
        materials_nonnull = materials_nonnull && vertices[i].material;

    // this doesn't interpolate materials. Ah well.
    if ( materials_nonnull )
        vertices[0].material->set_gl_state();

    glBegin(GL_TRIANGLES);

    glNormal3dv( &vertices[0].normal.x );
    glTexCoord2dv( &vertices[0].tex_coord.x );
    glVertex3dv( &vertices[0].position.x );

    glNormal3dv( &vertices[1].normal.x );
    glTexCoord2dv( &vertices[1].tex_coord.x );
    glVertex3dv( &vertices[1].position.x);

    glNormal3dv( &vertices[2].normal.x );
    glTexCoord2dv( &vertices[2].tex_coord.x );
    glVertex3dv( &vertices[2].position.x);

    glEnd();

    if ( materials_nonnull )
        vertices[0].material->reset_gl_state();
}

////////////////////////////////////////////////////////////
// Added Function: intersection test for ray and triangle //
////////////////////////////////////////////////////////////

Intersect Triangle::intersect (Ray ray) const
{
	// construct return value
    Intersect result = Intersect();
    result.intersect = false;

    // after transformaton, the center of the object should be (0,0,0)
    Vector3 E = invMat.transform_point(ray.e);        // ray eye coordinate
    Vector3 D = invMat.transform_vector(ray.d);       // ray direction vector

    // calculate all the parameters
    real_t a = vertices[0].position.x - vertices[1].position.x;
    real_t b = vertices[0].position.y - vertices[1].position.y;
    real_t c = vertices[0].position.z - vertices[1].position.z;

    real_t d = vertices[0].position.x - vertices[2].position.x;
    real_t e = vertices[0].position.y - vertices[2].position.y;
    real_t f = vertices[0].position.z - vertices[2].position.z;

    real_t g = D.x;
    real_t h = D.y;
    real_t i = D.z;

    real_t j = vertices[0].position.x - E.x;
    real_t k = vertices[0].position.y - E.y;
    real_t l = vertices[0].position.z - E.z;

    // value cache
    real_t ei_hf = e*i - h*f;
    real_t gf_di = g*f - d*i;
    real_t dh_eg = d*h - e*g;

    real_t M = a*ei_hf + b*gf_di + c*dh_eg;

    // check t intersection
    real_t ak_jb = a*k - j*b;
    real_t jc_al = j*c - a*l;
    real_t bl_kc = b*l - k*c;
    real_t t = -( f*ak_jb + e*jc_al +d*bl_kc ) / M;

    if (t < MIN_T_THRESHOLD) {
        return result;
    }

    // check gamma intersection
    real_t gamma = ( i*ak_jb + h*jc_al + g*bl_kc ) / M;
    if (gamma < 0.0 || gamma > 1.0) {
        return result;
    }

    // check beta intersection
    real_t beta = ( j*ei_hf + k*gf_di + l*dh_eg ) / M;
    if (beta < 0.0 || beta > 1.0-gamma) {
        return result;
    }

    // compute results at this point
    result.intersect = true;
    
    // get normal transformation matrix
    Matrix4 invMat_inv;
    make_transformation_matrix(&invMat_inv, position, orientation, scale);
    // calcualte world position
    Vector3 local_pos = E + t*D;
    result.position = invMat_inv.transform_point(local_pos);

    // calculate world normal
    Vector3 local_norm = normalize(vertices[0].normal + vertices[1].normal + vertices[2].normal);
    result.normal = normalize(normMat * local_norm);

    // calcualte material based on blending

    // calculate ambient color
    result.ambient = vertices[0].material->ambient * (1.0-beta-gamma) + 
                     vertices[1].material->ambient * beta + 
                     vertices[2].material->ambient * gamma;
    // calculate diffuse color
    result.diffuse = vertices[0].material->diffuse * (1.0-beta-gamma) + 
                     vertices[1].material->diffuse * beta + 
                     vertices[2].material->diffuse * gamma;
    // calculate specular color
    result.specular = vertices[0].material->specular * (1.0-beta-gamma) + 
                      vertices[1].material->specular * beta + 
                      vertices[2].material->specular * gamma;
    // calulate refractive index
    result.refractive_index = vertices[0].material->refractive_index * (1.0-beta-gamma) + 
                              vertices[1].material->refractive_index * beta + 
                              vertices[2].material->refractive_index * gamma;

    // calculate texture color at this point
    if (vertices[0].material->texture_filename.empty() &&
        vertices[1].material->texture_filename.empty() &&
        vertices[2].material->texture_filename.empty())
    {
        // if none of the material has texture, return white
        result.texture = Color3::White();
    }
    else 
    {
        // get texture coordinates first
        Vector2 tex = vertices[0].tex_coord * (1.0-beta-gamma) + 
                      vertices[1].tex_coord * beta + 
                      vertices[2].tex_coord * gamma;
        // Bug fix: here need to solve the texture wrap problem
        tex.x = tex.x - floor(tex.x);
        tex.y = tex.y - floor(tex.y);
        
        // get texture color
        // get the width and height
        int width, height;
        // set color accumulator
        Color3 color = Color3::Black();
        // set the first texture color
        vertices[0].material->get_texture_size(&width, &height);
        color += vertices[0].material->get_texture_pixel( (int)(tex.x*(real_t)(width)), (int)(tex.y*(real_t)(height)) ) * (1.0-beta-gamma);
        // set the first texture color
        vertices[1].material->get_texture_size(&width, &height);
        color += vertices[1].material->get_texture_pixel( (int)(tex.x*(real_t)(width)), (int)(tex.y*(real_t)(height)) ) * beta; 
        // set the first texture color
        vertices[2].material->get_texture_size(&width, &height);
        color += vertices[2].material->get_texture_pixel( (int)(tex.x*(real_t)(width)), (int)(tex.y*(real_t)(height)) ) * gamma;

        // set the final texture color
        result.texture = color;
    }

    // return the value
    return result;
}

} /* _462 */
