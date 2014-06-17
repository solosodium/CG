/**
 * @file sphere.cpp
 * @brief Function defnitions for the Sphere class.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#include "scene/sphere.hpp"
#include "application/opengl.hpp"

namespace _462 {

#define SPHERE_NUM_LAT 80
#define SPHERE_NUM_LON 100

#define SPHERE_NUM_VERTICES ( ( SPHERE_NUM_LAT + 1 ) * ( SPHERE_NUM_LON + 1 ) )
#define SPHERE_NUM_INDICES ( 6 * SPHERE_NUM_LAT * SPHERE_NUM_LON )
// index of the x,y sphere where x is lat and y is lon
#define SINDEX(x,y) ((x) * (SPHERE_NUM_LON + 1) + (y))
#define VERTEX_SIZE 8
#define TCOORD_OFFSET 0
#define NORMAL_OFFSET 2
#define VERTEX_OFFSET 5

static unsigned int Indices[SPHERE_NUM_INDICES];
static float Vertices[VERTEX_SIZE * SPHERE_NUM_VERTICES];

static void init_sphere()
{
    static bool initialized = false;
    if ( initialized )
        return;

    for ( int i = 0; i <= SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j <= SPHERE_NUM_LON; j++ ) {
            real_t lat = real_t( i ) / SPHERE_NUM_LAT;
            real_t lon = real_t( j ) / SPHERE_NUM_LON;
            float* vptr = &Vertices[VERTEX_SIZE * SINDEX(i,j)];

            vptr[TCOORD_OFFSET + 0] = lon;
            vptr[TCOORD_OFFSET + 1] = 1-lat;

            lat *= PI;
            lon *= 2 * PI;
            real_t sinlat = sin( lat );

            vptr[NORMAL_OFFSET + 0] = vptr[VERTEX_OFFSET + 0] = sinlat * sin( lon );
            vptr[NORMAL_OFFSET + 1] = vptr[VERTEX_OFFSET + 1] = cos( lat ),
            vptr[NORMAL_OFFSET + 2] = vptr[VERTEX_OFFSET + 2] = sinlat * cos( lon );
        }
    }

    for ( int i = 0; i < SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j < SPHERE_NUM_LON; j++ ) {
            unsigned int* iptr = &Indices[6 * ( SPHERE_NUM_LON * i + j )];

            unsigned int i00 = SINDEX(i,  j  );
            unsigned int i10 = SINDEX(i+1,j  );
            unsigned int i11 = SINDEX(i+1,j+1);
            unsigned int i01 = SINDEX(i,  j+1);

            iptr[0] = i00;
            iptr[1] = i10;
            iptr[2] = i11;
            iptr[3] = i11;
            iptr[4] = i01;
            iptr[5] = i00;
        }
    }

    initialized = true;
}

Sphere::Sphere()
    : radius(0), material(0) {}

Sphere::~Sphere() {}

void Sphere::render() const
{
    // create geometry if we haven't already
    init_sphere();

    if ( material )
        material->set_gl_state();

    // just scale by radius and draw unit sphere
    glPushMatrix();
    glScaled( radius, radius, radius );
    glInterleavedArrays( GL_T2F_N3F_V3F, VERTEX_SIZE * sizeof Vertices[0], Vertices );
    glDrawElements( GL_TRIANGLES, SPHERE_NUM_INDICES, GL_UNSIGNED_INT, Indices );
    glPopMatrix();

    if ( material )
        material->reset_gl_state();
}

//////////////////////////////////////////////////////////
// Added Function: intersection test for ray and sphere //
//////////////////////////////////////////////////////////

Intersect Sphere::intersect(Ray ray) const
{
	// Use the equations in fundamentals of computer graphics p76-77

    // construct return value
    Intersect result = Intersect();
	
	// after transformaton, the center of the object should be (0,0,0)
	Vector3 e = invMat.transform_point(ray.e);		           // ray eye coordinate
	Vector3 d = normalize(invMat.transform_vector(ray.d));    // ray direction vector
	
	// computer the discriminant if it's valid
	real_t discriminant = dot(d,e)*dot(d,e) - dot(d,d)*(dot(e,e)-radius*radius);
	
	// if the discriminant < 0, no intersection
	if (discriminant < 0.0) {
        result.intersect = false;
        return result;
	}
	else {
		// compute the normal matrix
		Matrix4 invMat_inv;
		make_transformation_matrix(&invMat_inv, position, orientation, scale);
		
		// compute the first t value and point
		real_t first_t = (-dot(d,e) + sqrt(discriminant)) / dot(d,d);
        Vector3 first_raw = e + first_t*d;
		Vector3 first = invMat_inv.transform_point(e + first_t*d);
		// compute the second t value and point
		real_t second_t = (-dot(d,e) - sqrt(discriminant)) / dot(d,d);
        Vector3 second_raw = e + second_t*d;
		Vector3 second = invMat_inv.transform_point(e + second_t*d);

        // determine the first or second point
        int CASE = 0;   // 0 & others: no intersection; 1: first; 2: decond
        // actual determination
        if (first_t < MIN_T_THRESHOLD && second_t < MIN_T_THRESHOLD) {
            CASE = 0;
        }
        // if only first t is valid
        else if (first_t > MIN_T_THRESHOLD && second_t < MIN_T_THRESHOLD) {
            CASE = 1;
        }
        // if only second t is valid
        else if (first_t < MIN_T_THRESHOLD && second_t > MIN_T_THRESHOLD) {
            CASE = 2;
        }
        // if both t are valid
        else {
            // if first one closer, return first, else second
            if (distance(ray.e, first) < distance(ray.e, second)) { CASE = 1; }
            else { CASE = 2; }
        }

        // switching between cases
        if (CASE == 1) {
            // compute some parameters
            result.intersect = true;
            result.position = first;
            result.normal = normalize(normMat * (e + first_t*d));
            //result.material = material;
            result.ambient = material->ambient;
            result.diffuse = material->diffuse;
            result.specular = material->specular;
            result.refractive_index = material->refractive_index;
            // compute texture color
            if (material->texture_filename.empty()) {
                result.texture = Color3::White();   // no texture case
            }
            else {
                real_t theta = acos(first_raw.z);
                real_t phi = atan2(first_raw.y, first_raw.x);
                real_t u = phi/2.0/3.141592654;
                real_t v = (3.141592654 - theta)/3.141592654;
                int width, height;
                material->get_texture_size(&width, &height);
                result.texture = material->get_texture_pixel( (int)(u*(real_t)(width)), (int)(v*(real_t)(height)) );
            }
            return result;                  // return first point
        }
        else if (CASE == 2) {
            // compute some parameters
            result.intersect = true;
            result.position = second;
            result.normal = normalize(normMat * (e + second_t*d));
            //result.material = material;
            result.ambient = material->ambient;
            result.diffuse = material->diffuse;
            result.specular = material->specular;
            result.refractive_index = material->refractive_index;
            // compute texture color
            if (material->texture_filename.empty()) {
                result.texture = Color3::White();   // sphere don't have texture
            }
            else {
                real_t theta = acos(second_raw.z);
                real_t phi = atan2(second_raw.y, second_raw.x);
                real_t u = phi/2.0/3.141592654;
                real_t v = (3.141592654 - theta)/3.141592654;
                int width, height;
                material->get_texture_size(&width, &height);
                result.texture = material->get_texture_pixel( (int)(u*(real_t)(width)), (int)(v*(real_t)(height)) );
            }
            return result;                  // return first point
        }
        else {
            // no valid intersection found
            result.intersect = false;
            return result;
        }
	}
}

} /* _462 */

