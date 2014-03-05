/**
 * @file project.cpp
 * @brief OpenGL project
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "p1/project.hpp"

// use this header to include the OpenGL headers
// DO NOT include gl.h or glu.h directly; it will not compile correctly.
#include "application/opengl.hpp"

// A namespace declaration. All proejct files use this namespace.
// Add this declration (and its closing) to all source/headers you create.
// Note that all #includes should be BEFORE the namespace declaration.
namespace _462 {

// definitions of functions for the OpenglProject class



// constructor, invoked when object is created
OpenglProject::OpenglProject()
{
    // TODO any basic construction or initialization of members
    // Warning: Although members' constructors are automatically called,
    // ints, floats, pointers, and classes with empty contructors all
    // will have uninitialized data!

    
}

// destructor, invoked when object is destroyed
OpenglProject::~OpenglProject()
{
    // TODO any final cleanup of members
    // Warning: Do not throw exceptions or call virtual functions from deconstructors!
    // They will cause undefined behavior (probably a crash, but perhaps worse).

    // For Mesh
    delete [] mesh_normals;
    
    // For Heightmap
    delete [] height_vertices;
    delete [] height_triangles;
    delete [] height_normals;
}

/**
 * Initialize the project, doing any necessary opengl initialization.
 * @param camera An already-initialized camera.
 * @param scene The scene to render.
 * @return true on success, false on error.
 */
bool OpenglProject::initialize( Camera* camera, Scene* scene )
{
    // copy scene
    this->scene = *scene;

    // TODO opengl initialization code and precomputation of mesh/heightmap

    // OpenGL state enabling
    glEnable ( GL_DEPTH_TEST );             // enable OpenGL depth option to handle vertex overlay
    glEnable( GL_NORMALIZE );               // enable OpenGL normalize stuff
    glEnable( GL_LIGHTING );                // enables lighting
    glEnable( GL_LIGHT0 );                  // enables light 0
    // Setup the light
    GLfloat light_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = { 0.0, 20.0, 1.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    // Setup OpenGL buffer clear parameters
    glClearColor (0.0, 0.0, 0.0, 0.0);      // when clear color buffer, set to black
    glClearDepth (1.0);                     // when slear depth buffer, set to 1.0 (not sure about this)

    /**
     * Fix: projection matrix here?
     */
    
    // Initialize the camera with projection matrix
    glMatrixMode( GL_PROJECTION );      // select the projection matrix
    glLoadIdentity();                   // set the matrix to identity (identity)
    // set up the camera using the camera parameters (provided)
    gluPerspective(camera->get_fov_degrees(), camera->get_aspect_ratio(), camera->get_near_clip(), camera->get_far_clip());
    // change back to modelview
    glMatrixMode( GL_MODELVIEW );

    // Initialize Mesh
    this->initializeMesh();
    // Initialize Heightmap
    this->initializeHeightmap();

    return true;
}

/**
 * Clean up the project. Free any memory, etc.
 */
void OpenglProject::destroy()
{
    // TODO any cleanup code, e.g., freeing memory
}

/**
 * Perform an update step. This happens on a regular interval.
 * @param dt The time difference from the previous frame to the current.
 */
void OpenglProject::update( real_t dt )
{
    // update our heightmap
    scene.heightmap->update( dt );

    // TODO any update code, e.g. commputing heightmap mesh positions and normals

    // Update heightmap vertices
    this->calculateVertices();

	// Update Heightmap normals
	this->calculateNormals();
}

/**
 * Clear the screen, then render the mesh using the given camera.
 * @param camera The logical camera to use.
 * @see math/camera.hpp
 */
void OpenglProject::render( const Camera* camera )
{
    // TODO render code

    // Clear the buffers
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// Mesh is drawn here
    this->drawMesh (camera);
    /**
	 * This should be in the update scope
	 */
	// Update Heightmap normals
    //this->calculateNormals();
    
	// Heightmap is drawn here
    this->drawHeightmap (camera);

	// Flush everything
    glFlush();
}

/***************************************/
/* Implementation of Utility Functions */ 
/***************************************/

/**
 * Function to initialize mesh with normal vectors array
 */
void OpenglProject::initializeMesh ()
{ 
    // Data initilaization
    mesh_normals = new Vector3 [scene.mesh.num_vertices];       // normal for each vertex (3 coordinates)   
    // initialize the array to zero first
    for (int k=0; k<(int)scene.mesh.num_vertices; k++) {
        mesh_normals[k] = Vector3::Zero;
    }
    // compute normal array
    for (int k=0; k<(int)scene.mesh.num_triangles; k++)
    {
        // get the indices of three points of a triangle
        unsigned int* index = scene.mesh.triangles[k].vertices;
        // compute the normal 
        Vector3 normal = normalize(cross(scene.mesh.vertices[index[1]] - scene.mesh.vertices[index[0]], scene.mesh.vertices[index[2]] - scene.mesh.vertices[index[0]]));
        // assign the normal vector to appropraite vertices
        for (int n=0; n<3; n++) {
            mesh_normals[index[n]] += normal;
        }
    }
	// normalize all the nomral vectors
	for (int k=0; k<(int)scene.mesh.num_vertices; k++) {
        mesh_normals[k] = normalize(mesh_normals[k]);
    }
}

/**
 * Function to initialize Heightmap 
 */
void OpenglProject::initializeHeightmap ()
{
    // Data initilaization
    height_size = 64;                                                       // height map side size
    height_vertices = new Vector3 [height_size*height_size];                // heightmap vertices array
    height_triangles = new Triangle [(height_size-1)*(height_size-1)*2];    // heightmap triangles array
    height_normals = new Vector3 [height_size*height_size];                 // heightmap normal array
    // initialize vertices
    this->calculateVertices();
    // initialize triangle mesh and normals
    for (int k=0; k<(int)(height_size-1); k++) {
        for (int n=0; n<(int)(height_size-1); n++)
        {
            // first triangle
            height_triangles[(k*(height_size-1)+n)*2 + 0].vertices[0] = k*height_size+n;
            height_triangles[(k*(height_size-1)+n)*2 + 0].vertices[1] = (k+1)*height_size+n;
            height_triangles[(k*(height_size-1)+n)*2 + 0].vertices[2] = (k+1)*height_size+n+1;
            // second triangle
            height_triangles[(k*(height_size-1)+n)*2 + 1].vertices[0] = k*height_size+n;
            height_triangles[(k*(height_size-1)+n)*2 + 1].vertices[1] = (k+1)*height_size+n+1;
            height_triangles[(k*(height_size-1)+n)*2 + 1].vertices[2] = k*height_size+n+1;
        }
    }
    // initialize normals
    this->calculateNormals();
}

/**
 * Function to calculate vertices for heightmap
 */
void OpenglProject::calculateVertices ()
{
    // initialize vertices
    for (int k=0; k<(int)height_size; k++) {
        for (int n=0; n<(int)height_size; n++) 
        {
            height_vertices[k*height_size + n].x = -1.0 + 2.0*(double)(k)/(double)(height_size);      // x
            height_vertices[k*height_size + n].z = -1.0 + 2.0*(double)(n)/(double)(height_size);      // z
            Vector2 loc (height_vertices[k*height_size + n].x, height_vertices[k*height_size + n].z);
            height_vertices[k*height_size + n].y = scene.heightmap->compute_height(loc);              // y
        }
    }
}

/**
 * Function to calculate normals for heightmap
 */
void OpenglProject::calculateNormals ()
{
    // initialize normals to zero
    for (int k=0; k<(int)(height_size*height_size); k++) {
        height_normals[k] = Vector3::Zero;
    }
    // compute the normals based on the triangle array
    for (int k=0; k<(int)(height_size-1); k++) {
        for (int n=0; n<(int)(height_size-1); n++)
        {
            // first triangle index 
            unsigned int* t1 = height_triangles[(k*(height_size-1)+n)*2 + 0].vertices;
            // second triangle
            unsigned int* t2 = height_triangles[(k*(height_size-1)+n)*2 + 1].vertices;
            // compute the normal 
            Vector3 t1_normal = normalize(cross(height_vertices[t1[2]] - height_vertices[t1[0]], height_vertices[t1[1]] - height_vertices[t1[0]]));
            Vector3 t2_normal = normalize(cross(height_vertices[t2[2]] - height_vertices[t2[0]], height_vertices[t2[1]] - height_vertices[t2[0]]));
            // assign the normal vector to appropraite vertices in t1 and t2
            for (int m=0; m<3; m++) {
                height_normals[t1[m]] += t1_normal;
                height_normals[t2[m]] += t2_normal;
            }
        }
    }
	// normalize all the nomral vectors
	for (int k=0; k<(int)(height_size*height_size); k++) {
        height_normals[k] = normalize(height_normals[k]);
    }
}

/**
 * Function to draw the Mesh
 */
void OpenglProject::drawMesh (const Camera* camera)
{
    // start drawing
		glLoadIdentity();       // Clears the matrix
		gluLookAt(camera->get_position().x, camera->get_position().y, camera->get_position().z, camera->get_position().x + camera->get_direction().x, camera->get_position().y + camera->get_direction().y, camera->get_position().z + camera->get_direction().z, camera->get_up().x, camera->get_up().y, camera->get_up().z);
        
		glPushMatrix();
		
		// transform the object
        glTranslatef(scene.mesh_position.position.x, scene.mesh_position.position.y, scene.mesh_position.position.z); 
        glRotatef(scene.mesh_position.orientation.w, scene.mesh_position.orientation.x, scene.mesh_position.orientation.y, scene.mesh_position.orientation.z); 
        glScalef(scene.mesh_position.scale.x, scene.mesh_position.scale.y, scene.mesh_position.scale.z);
        // enable
        glEnableClientState( GL_VERTEX_ARRAY ); // enable using vertex array
        glEnableClientState( GL_NORMAL_ARRAY ); // enable using normal array
			
			// draw in between
            glColor3f(1.0f, 0.2f, 0.1f);           // set the color of the object
            GLfloat diffuse_intensity[] = {1.0f, 0.2f, 0.0f, 1.0f};     // set the diffuse material
            GLfloat specular_intensity[] = {0.1f, 0.1f, 0.1f, 1.0f};    // set the specular material
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_intensity);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_intensity);
            glVertexPointer(3, GL_DOUBLE, 0, scene.mesh.vertices);
            glNormalPointer(GL_DOUBLE, 0, mesh_normals);
            glDrawElements(GL_TRIANGLES, scene.mesh.num_triangles*3, GL_UNSIGNED_INT, scene.mesh.triangles);
			
		// disable
        glDisableClientState( GL_VERTEX_ARRAY ); // disable using vertex array
        glDisableClientState( GL_NORMAL_ARRAY ); // disable using normal array
		
		glPopMatrix();
}

/**
 * Function to draw the Heightmap
 */
void OpenglProject::drawHeightmap (const Camera* camera)
{
    // start drawing
		glLoadIdentity();       // Clears the matrix
		gluLookAt(camera->get_position().x, camera->get_position().y, camera->get_position().z, camera->get_position().x + camera->get_direction().x, camera->get_position().y + camera->get_direction().y, camera->get_position().z + camera->get_direction().z, camera->get_up().x, camera->get_up().y, camera->get_up().z);
        
		// transform the object

		glPushMatrix();

        glTranslatef(scene.heightmap_position.position.x, scene.heightmap_position.position.y, scene.heightmap_position.position.z); 
        glRotatef(scene.heightmap_position.orientation.w, scene.heightmap_position.orientation.x, scene.heightmap_position.orientation.y, scene.heightmap_position.orientation.z); 
        glScalef(scene.heightmap_position.scale.x, scene.heightmap_position.scale.y, scene.heightmap_position.scale.z);
        // enable
        glEnableClientState( GL_VERTEX_ARRAY ); // enable using vertex array
        glEnableClientState( GL_NORMAL_ARRAY ); // enable using normal array
			
            // draw in between
            glColor3f(0.0f, 0.3f, 1.0f);
            GLfloat diffuse_intensity[] = {0.0f, 0.3f, 1.0f, 1.0f};     // set the diffuse material
            GLfloat specular_intensity[] = {0.0f, 0.3f, 1.0f, 1.0f};    // set the specular material
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_intensity);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_intensity);
            glVertexPointer(3, GL_DOUBLE, 0, height_vertices);
            glNormalPointer(GL_DOUBLE, 0, height_normals);
            glDrawElements(GL_TRIANGLES, (height_size-1)*(height_size-1)*6, GL_UNSIGNED_INT, height_triangles);
		// disable
        glDisableClientState( GL_VERTEX_ARRAY ); // disable using vertex array
        glDisableClientState( GL_NORMAL_ARRAY ); // disable using normal array

		glPopMatrix();
}

} /* _462 */
