/**
 * @file main.hpp
 * @brief OpenGL project application and main function.
 *
 * @author Eric Butler (edbutler)
 */

#include "application/application.hpp"
#include "application/camera_roam.hpp"
#include "p1/project.hpp"
#include "scene/mesh.hpp"
#include "scene/heightmap.hpp"

namespace _462 {

#define KEY_SCREENSHOT SDLK_f

static const int WIDTH = 1024;
static const int HEIGHT = 768;
static const real_t FPS = 60.0;
static const char* TITLE = "15462 Project - OpenGL";

class OpenglApplication : public Application
{
public:

    OpenglApplication() { }
    virtual ~OpenglApplication() { }

    virtual bool initialize();
    virtual void destroy();
    virtual void update( real_t dt );
    virtual void render();
    virtual void handle_event( const SDL_Event& event );

    // the project
    OpenglProject project;
    // the camera
    CameraRoamControl camera_control;

    // the triangle mesh we're rendering
    Mesh mesh;
    // the heightmap we're using
    WaterSurface heightmap;
    // the scene given to the project
    Scene scene;
};

bool OpenglApplication::initialize()
{
    bool rv = true;

    // set camera paremeters
    // aspect ratio is set by the application, so we can ignore it
    Camera* camera = &camera_control.camera;
    camera->fov = PI / 4.0;
    camera->near_clip = 0.1;
    camera->far_clip = 100.0;
    camera->position = Vector3( 0.0, 2.5, 5.5 );
    camera->orientation = Quaternion::Identity;
    camera->pitch( -PI / 5.0 );

    // hacked-in scene (DO NOT CHANGE)

    mesh.filename = "models/pool.obj";
    rv = rv && mesh.load();

    // copy vertices from loaded mesh
    scene.mesh.num_vertices = mesh.vertices.size();
    scene.mesh.vertices = new Vector3[scene.mesh.num_vertices];
    scene.mesh.num_triangles = mesh.triangles.size();
    scene.mesh.triangles = new Triangle[scene.mesh.num_triangles];
    for ( size_t i = 0; i < scene.mesh.num_vertices; ++i ) {
        scene.mesh.vertices[i] = mesh.vertices[i].position;
    }
    for ( size_t i = 0; i < scene.mesh.num_triangles; ++i ) {
        for ( size_t j = 0; j < 3; ++j ) {
            scene.mesh.triangles[i].vertices[j] = mesh.triangles[i].vertices[j];
        }
    }

    scene.mesh_position.position = Vector3( 0.0, -1.0, 0.0 );
    scene.mesh_position.orientation = Quaternion( Vector3::UnitY, PI / 3.6 );
    scene.mesh_position.scale = Vector3( 0.5, 0.5, 0.5 );

    // create a watersurface
    WaterSurface::WavePoint* p;

    heightmap.wave_points.push_back(WaterSurface::WavePoint());
    p = &heightmap.wave_points.back();
    p->position = Vector2(.71,.78);
    p->falloff = 2;
    p->amplitude = .15;
    p->timerate = -2.0*PI;
    p->period = 8*PI;

    heightmap.wave_points.push_back(WaterSurface::WavePoint());
    p = &heightmap.wave_points.back();
    p->position = Vector2(-.21,-.35);
    p->falloff = 2;
    p->amplitude = .2;
    p->timerate = -8.0/3.0*PI;
    p->period = 10*PI;

    scene.heightmap = &heightmap;
    scene.heightmap_position.position = Vector3( 0.0, -0.4, 0.0 );
    scene.heightmap_position.orientation = Quaternion( Vector3::UnitY, PI / 3.6 );
    scene.heightmap_position.scale = Vector3( 1.5, 0.5, 1.5 );

    rv = rv &&  project.initialize( &camera_control.camera, &scene );

    return rv;
}

void OpenglApplication::destroy()
{
    project.destroy();
    delete [] scene.mesh.vertices;
    delete [] scene.mesh.triangles;
}

void OpenglApplication::update( real_t dt )
{
    // update camera
    camera_control.update( dt );
    // update project
    project.update( dt );
}

void OpenglApplication::render()
{
    // adjust camera aspect
    int width, height;
    this->get_dimension( &width, &height );
    camera_control.camera.aspect = real_t( width ) / real_t( height );

    // render
    project.render( &camera_control.camera );
}

void OpenglApplication::handle_event( const SDL_Event& event )
{
    camera_control.handle_event( this, event );

    // open to add more event handlers
    switch ( event.type )
    {
    case SDL_KEYDOWN:
        switch ( event.key.keysym.sym )
        {
        case KEY_SCREENSHOT:
            take_screenshot();
            break;
        default:
            break;
        }
    default:
        break;
    }
}

} /* _462 */

using namespace _462;

static void print_usage( const char* progname )
{
    std::cout << "Usage : " << progname << "\n";
}

int main( int argc, char* argv[] )
{
    // start a new application
    OpenglApplication app;
    if ( argc > 1 ) {
        print_usage( argv[0] );
        return 1;
    }
    return Application::start_application( &app, WIDTH, HEIGHT, FPS, TITLE );
}
