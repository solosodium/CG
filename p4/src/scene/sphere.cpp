/**
 * @file sphere.cpp
 * @brief Function defnitions for the Sphere class.
 * DO NOT DISTRIBUTE
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

            real_t sinlat = sqrt(1- pow(cos(lat), 2));

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

bool Sphere::intersect(const Ray& r, Intersection* info)
{
    real_t a = squared_length(r.d);
    real_t b = real_t(2)*dot(r.d, r.e);
    real_t c = squared_length(r.e) - radius * radius;
    real_t disc = b*b - real_t(4)*a*c;

    if (disc < real_t(0)) return false;

    int df = b < 0 ? -1 : 1;
    real_t q = (-b + df*sqrt(disc))/real_t(2);

    real_t t0 = q/a;
    real_t t1 = c/q;

    if (t0 > t1) std::swap(t0, t1);

    if (t1 < real_t(0)) return false;

    info->t = t0 > real_t(0) ? t0 : t1;
    info->normal = r.d*(info->t+1e-9) + r.e;
    return true;
}

void Sphere::get_extended_info(Intersection* info)
{
	info->normal = normalize(info->normal);
	info->diffuse = material->diffuse;
	info->ambient = material->ambient;
	info->specular = material->specular;
	info->refractive_index = material->refractive_index;

  int w, h;
  material->get_texture_size(&w, &h);

  if (w == 0 || h == 0) {
    info->texColor = Color3::White();
  } else {
		const Vector3& d = info->normal;

		real_t theta = acos(d.y);
		real_t phi = atan2(d.x, d.z);

		real_t u = phi/(2*PI);
		real_t v = (PI - theta)/PI;
		info->texColor = material->get_texture_pixel(int(u*w)%w, int(v*h)%h);
	}
}

/**
 * TODO: implement get bounds function
 */
void Sphere::get_bounds (Bounds& bounds) const
{
	// sphere bounds more complicated for the rotation and scale
	// first get a bounding box in world space
	Vector3 corners [8];
	corners[0] = invMat_inv.transform_point(Vector3(radius, radius, radius));
	corners[1] = invMat_inv.transform_point(Vector3(radius, radius, -radius));
	corners[2] = invMat_inv.transform_point(Vector3(radius, -radius, radius));
	corners[3] = invMat_inv.transform_point(Vector3(radius, -radius, -radius));
	corners[4] = invMat_inv.transform_point(Vector3(-radius, radius, radius));
	corners[5] = invMat_inv.transform_point(Vector3(-radius, radius, -radius));
	corners[6] = invMat_inv.transform_point(Vector3(-radius, -radius, radius));
	corners[7] = invMat_inv.transform_point(Vector3(-radius, -radius, -radius));
	// loop through thr corners to get the bounds in world space
	for (int k=0; k<8; k++) 
	{
		// fit into either one of the bounds
		if (corners[k].x < bounds.x_min) {
			bounds.x_min = corners[k].x;
		}
		else if (corners[k].x > bounds.x_max) {
			bounds.x_max = corners[k].x;
		}
		if (corners[k].y < bounds.y_min) {
			bounds.y_min = corners[k].y;
		}
		else if (corners[k].y > bounds.y_max) {
			bounds.y_max = corners[k].y;
		}
		if (corners[k].z < bounds.z_min) {
			bounds.z_min = corners[k].z;
		}
		else if (corners[k].z > bounds.z_max) {
			bounds.z_max = corners[k].z;
		}
	}
    // settle down the center values as well
    bounds.x_center = (bounds.x_min + bounds.x_max) / 2.0;
    bounds.y_center = (bounds.y_min + bounds.y_max) / 2.0;
    bounds.z_center = (bounds.z_min + bounds.z_max) / 2.0;
	// MARK: this method will introduce a little extra volume for scaled and rotated sphere
};

} /* _462 */
