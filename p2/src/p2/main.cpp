
#include "application/application.hpp"
#include "application/camera_roam.hpp"
#include "p2/project.hpp"
#include "application/imageio.hpp"
#include "application/scene_loader.hpp"
#include "application/opengl.hpp"
#include "scene/scene.hpp"
#include <SDL.h>

#define CUBEMAP_PATH ("./cubemaps/")

namespace _462 {

#define KEY_SCREENSHOT SDLK_f
#define KEY_SUBDIVIDE SDLK_y

// pretty sure these are sequential, but use an array just in case
static const GLenum LightConstants[] = {
    GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3,
    GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7
};
static const size_t NUM_GL_LIGHTS = 8;

static const int WIDTH = 1024;
static const int HEIGHT = 768;
static const real_t FPS = 60.0;
static const char* TITLE = "15462 Project - OpenGL";

class GlslApplication : public Application, public SceneRenderer
{
public:

    GlslApplication() { }
    virtual ~GlslApplication() { }

    virtual bool initialize();
    virtual void destroy();
    virtual void update( real_t dt );
    virtual void render();
    virtual void handle_event( const SDL_Event& event );

    virtual void render_scene() const;
    virtual const Camera* get_camera() const {
        return &camera_control.camera;
    }

    // the project
    GlslProject project;
    // the camera
    CameraRoamControl camera_control;
    // the rendered scene
    //Scene scene;
};

bool GlslApplication::initialize()
{
    int width, height;
    this->get_dimension( &width, &height );

    // copy over loaded camera
    camera_control.camera = scene.camera;

    try {
        const Scene::MaterialList& materials = scene.get_materials();
        const Scene::MeshList& meshes = scene.get_meshes();

        // load all textures
        for ( size_t i = 0; i < materials.size(); ++i ) {
            if ( !materials[i]->load() || !materials[i]->create_gl_data() ) {
                std::cout << "Error loading texture, aborting.\n";
                return false;
            }
        }

        // load all meshes
        for ( size_t i = 0; i < meshes.size(); ++i ) {
            if ( !meshes[i]->load() || !meshes[i]->create_gl_data() ) {
                std::cout << "Error loading mesh, aborting.\n";
                return false;
            }
        }

    } catch ( std::bad_alloc const& ) {
        std::cout << "Out of memory error while initializing scene\n.";
        return false;
    }
    project.scene = &scene;
    return project.initialize( this, width, height );
}

void GlslApplication::destroy()
{
    project.destroy();
}

void GlslApplication::update( real_t dt )
{
    // update camera
    camera_control.update( dt );
}

void GlslApplication::render()
{
    // adjust camera aspect
    int width, height;
    this->get_dimension( &width, &height );
    glViewport( 0, 0, width, height );
    camera_control.camera.aspect = real_t( width ) / real_t( height );

    // render
    project.render();
}

void GlslApplication::handle_event( const SDL_Event& event )
{
    camera_control.handle_event( this, event );

    switch ( event.type )
    {
    case SDL_KEYDOWN:
        switch ( event.key.keysym.sym )
        {
        case KEY_SCREENSHOT:
            take_screenshot();
            break;
	case KEY_SUBDIVIDE:
	    project.scene->subdivide_geometries();
	    break;
        default:
            break;
        }
    default:
        break;
    }
}

void GlslApplication::render_scene() const
{
    // backup state so it doesn't mess up students stuff
    glPushAttrib( GL_ALL_ATTRIB_BITS );
    glPushClientAttrib( GL_CLIENT_ALL_ATTRIB_BITS );

    // we set most gl state each time since it all gets popped each frame
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);

    // set camera transform

    const Camera& camera = camera_control.camera;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(camera.get_fov_degrees(),
                   camera.get_aspect_ratio(),
                   camera.get_near_clip(),
                   camera.get_far_clip());

    const Vector3& campos = camera.get_position();
    const Vector3 camref = camera.get_direction() + campos;
    const Vector3& camup = camera.get_up();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    gluLookAt(campos.x, campos.y, campos.z,
              camref.x, camref.y, camref.z,
              camup.x,  camup.y,  camup.z);
    // set light data
    float arr[4];
    arr[3] = 1.0; // w is always 1

    scene.ambient_light.to_array( arr );
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, arr );

    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );

    const Scene::PointLightList& lights = scene.get_lights();

    size_t N_LIGHTS = lights.size();

    for (size_t i = 0; i < NUM_GL_LIGHTS && i < N_LIGHTS; i++)
    {
        const PointLight& light = lights[i];
        glEnable(LightConstants[i]);
        light.color.to_array(arr);
        glLightfv(LightConstants[i], GL_DIFFUSE, arr);
        glLightfv(LightConstants[i], GL_SPECULAR, arr);
        glLightf(LightConstants[i], GL_CONSTANT_ATTENUATION,
		 (GLfloat)light.attenuation.constant);
        glLightf(LightConstants[i], GL_LINEAR_ATTENUATION,
		 (GLfloat)light.attenuation.linear);
        glLightf(LightConstants[i], GL_QUADRATIC_ATTENUATION,
		 (GLfloat)light.attenuation.quadratic);
        light.position.to_array(arr);
        glLightfv(LightConstants[i], GL_POSITION, arr);
    }
    // render each object

    const Scene::GeometryList& geometries = scene.get_geometries();

    for ( size_t i = 0; i < geometries.size(); ++i ) {
        const Geometry& geom = *geometries[i];
        Vector3 axis;
        real_t angle;

        glPushMatrix();

        glTranslated( geom.position.x, geom.position.y, geom.position.z );
        geom.orientation.to_axis_angle( &axis, &angle );
        glRotated( angle * ( 180.0 / PI ), axis.x, axis.y, axis.z );
        glScaled( geom.scale.x, geom.scale.y, geom.scale.z );

        geom.render();

        glPopMatrix();
    }
    glPopMatrix();
    glPopClientAttrib();
    glPopAttrib();
}

} /* _462 */

using namespace _462;

static void print_usage( const char* progname )
{
    std::cout << "Usage : "
	      << progname << " <scene filename> [<cubmap name>]\n";
}

int main( int argc, char* argv[] )
{
    GlslApplication app;
    const char* cubemap_name;
    // parse arguments
    if ( argc < 2 ) {
        print_usage( argv[0] );
        return 2;
    }
    else if (argc == 3)
	cubemap_name = argv[2];
    else
	cubemap_name = "teide";

    // load the given scene
    if ( !load_scene( &app.scene, argv[1] ) ) {
        std::cout << "Error loading scene " << argv[1] << ". Aborting.\n";
        return 1;
    }

    std::cout << "Press 'y' to subdivide." << std::endl;

    const char* face_files[6] = {
	"/posx.png",
	"/negx.png",
	"/posy.png",
	"/negy.png",
	"/posz.png",
	"/negz.png",
    };
    for (int i = 0; i < 6; i++)
	app.project.cubemap_faces[i] = (std::string(CUBEMAP_PATH) +
					cubemap_name) + face_files[i];
    // start a new application
    return Application::start_application( &app, WIDTH, HEIGHT, FPS, TITLE );
}
