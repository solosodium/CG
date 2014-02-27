#include "scene/mesh.hpp"

// Added libraries
#include <math.h>
#include "math/vector.hpp"

namespace _462 {

// Define edge data structure (winged)
struct MyEdge
{
	// contructor take in edge and left triangle info
	MyEdge (int start_in, int end_in, int left_in, Vector3 norm_left_in)
	: right(-1), tri_count(1), create(-1)
	{
		start = start_in;
		end = end_in;
		left = left_in;
		norm_left = norm_left_in;
	}
	// member data (before creating new vertex)
	int start;			// start of the edge
	int end;			// end of the edge
	int left;			// left point (if there)
	int right;			// right point (if there)
	int tri_count;		// number of adjacent triangles
	
	// member data (after creating new vertex)
	int create;			// the created vertex index
	
	// this member data is for adaptive subdivision
	Vector3 norm_left;	// normal vector of left triangle
	Vector3 norm_right;	// normal vector of right triangle
};

// A customed type for the Edge list
typedef std::vector <MyEdge> EdgeList;

// Define a triangle pointing to three edges
struct MyTriangle
{
	MyEdge * one;		// first edge
	MyEdge * two;		// second edge
	MyEdge * three;		// third edge
};

bool Mesh::subdivide()
{
    /*
      You should implement loop subdivision here.

      Triangles are stored in an std::vector<MeshTriangle> in 'triangles'.
      Vertices are stored in an std::vector<MeshVertex> in 'vertices'.

      Check mesh.hpp for the Mesh class definition.
     */

	// NOTE: CONTROL PARAMETER
	bool ADAPTIVE = false;			// control for adaptive subdivision
	float THRESHOLD = 0.999;		// threshold for adaptive subdivision
	
	/************************************************************/
	/************** Stage 1: Preprocess Data ********************/
	/************************************************************/
	
	/*
		What I want to do here is do use the original data structure 
		as much as possible. At the same time, I also need to build 
		an efficient adjacency data structure for fast traversal.
	*/
	
	// Get the information about the existing vertices and triangles data
	int vert_ori_size = vertices.size();	// the original size of the vertices
	int tri_ori_size = triangles.size();	// the original size of the triangles
	
	// Define the edge list data structure
	EdgeList edge_list;						// instantiation of EdgeList
	edge_list.reserve(triangles.size()*3);	// reserve the memory (large enough)
	
	// Define the triangle list as a 1D array
	MyTriangle* tri_list = new MyTriangle [tri_ori_size];
	
	// Travel through all the triangles ONCE to build the winged edges (edge_list)
	// This also build the tri_list which points to all its three edges
	for (int k=0; k<tri_ori_size; k++)
	{
		// get the three points from a triangle
		int first = (int) triangles[k].vertices[0];
		int second = (int) triangles[k].vertices[1];
		int third = (int) triangles[k].vertices[2];
		
		// set up booleans for the flags for the three edges
		bool first_exist = false;
		bool second_exist = false;
		bool third_exist = false;
        
        // compute the normal vector of the triangle
        Vector3 tri_norm = normalize(vertices[first].normal
                                     + vertices[second].normal
                                     + vertices[third].normal);
		
		// loop through edge_list to find possible hits (reversed)
		for (int n=0; n<(int)edge_list.size(); n++)
		{
			// search for reverse edges,
            // if no reverse edge, means this edge does not exist
			if (edge_list[n].start == second && edge_list[n].end == first)
            {
				// First edge reverse found
				first_exist = true;					// set the flag
				edge_list[n].right = third;			// initialize right point
				edge_list[n].tri_count = 2;			// two triangles full
				// set the normal
				edge_list[n].norm_right = tri_norm;
				// set pointer in triangle list
				tri_list[k].one = &edge_list[n];
			}
			else if (edge_list[n].start == third && edge_list[n].end == second)
            {
				// Second edge reverse found
				second_exist = true;			// set the flag
				edge_list[n].right = first;		// initialize right point
				edge_list[n].tri_count = 2;		// two triangles full
				// set the normal
				edge_list[n].norm_right = tri_norm;
				// set pointer in triangle list
				tri_list[k].two = &edge_list[n];
			}
			else if (edge_list[n].start == first && edge_list[n].end == third)
            {
				// Third edge reverse found
				third_exist = true;				// set the flag
				edge_list[n].right = second;	// initialize right point
				edge_list[n].tri_count = 2;		// two triangles full
				// set the normal
				edge_list[n].norm_right = tri_norm;
				// set pointer in triangle list
				tri_list[k].three = &edge_list[n];
			}
			
			// if all the three reverse are found, no need for looking
			if (first_exist && second_exist && third_exist) {
				break;
			}
		}
		
		// if reverse edge is not found, need to add new edge
		if (!first_exist) {
			// add first edge
			MyEdge edge (first, second, third, tri_norm);
			edge_list.push_back(edge);
			// set pointer in triangle list
			tri_list[k].one = &edge_list.back();
		}
		if (!second_exist) {
			// add second edge
			MyEdge edge (second, third, first, tri_norm);
			edge_list.push_back(edge);
			// set pointer in triangle list
			tri_list[k].two = &edge_list.back();
		}
		if (!third_exist) {
			// add third edge
			MyEdge edge (third, first, second, tri_norm);
			edge_list.push_back(edge);
			// set pointer in triangle list
			tri_list[k].three = &edge_list.back();
		}	
	}
	
	// TEST: output the vertices, triangles and edges information
	std::cout << "Old vertices = "
              << vert_ori_size << "; "
              << "Old traingles = "
              << tri_ori_size << std::endl;
	
	// The number of the edge is the number of the new vertices.
	// Build the new vertices adjacency data structure upon this.
	
	// Vertices neighbor data structure (2D array)
    // the x dimension of vertext adjacency array is original vertex size (even)
	int vert_adj_dimen_x = vert_ori_size;
    // the y dimension of vertext adjacency array is new vertext size (odd+even)
	int vert_adj_dimen_y = vert_ori_size + edge_list.size();
	
	// if is not for adaptive subdivision, the vert_adj_dimen_y can be
    // reduced to only edge_list.size() the reason is that if all edges
    // are splited, even vertices only neighbor with odd ones
	
	// create the adjacency 2D array
	bool ** vert_adj = new bool * [vert_adj_dimen_x];
	// loop through the 2D array to 
	for (int k=0; k<vert_adj_dimen_x; k++) 
	{
		vert_adj[k] = new bool[vert_adj_dimen_y];
		// initialize all with values to false
		for (int n=0; n<vert_adj_dimen_y; n++) {
			vert_adj[k][n] = false;
		}
	}
	
	/************************************************************/
	/************ Stage 2: Compute Odd Vertices *****************/
	/************************************************************/
	
	/*
		This will loop through all the edges, and do:
		1. Compute the position of all the odd vertices
		2. Update and form new triangles based on the created vertices
		3. Set the adjacency matrix for vertices
	*/
	
	// Loop through all the available edges, compute odd vertices
	for (int k=0; k<(int)edge_list.size(); k++)
	{	
		// create the new vertex
		MeshVertex vertex;
		
		// if this edge is a boundary edge
		if (edge_list[k].tri_count == 1)
		{
			// calculate the position
			vertex.position = vertices[edge_list[k].start].position * 0.5
                              + vertices[edge_list[k].end].position * 0.5;
			// push this vertex to the vertex array
			vertices.push_back(vertex);
			// set the index of this vertex in the vertices vector at edge
			edge_list[k].create = vertices.size() - 1;
		}
		// if this edge is an inner edge
		else if (edge_list[k].tri_count == 2)
		{
			// if the edge has two winged triangles, we need adaptive subdivision
			
			// for adaptive subdivison
			if (ADAPTIVE)
			{
				// the creation of the new vertex is valid if the dot product is less than threshold
				if (dot(edge_list[k].norm_left, edge_list[k].norm_right) < THRESHOLD)
				{
					// calculate the position
					vertex.position = vertices[edge_list[k].start].position * 0.375
                                      + vertices[edge_list[k].end].position * 0.375
                                      + vertices[edge_list[k].left].position * 0.125
                                      + vertices[edge_list[k].right].position * 0.125;
					// push this vertex to the vertex array
					vertices.push_back(vertex);
					// set the index of this vertex in the vertices vector at edge
					edge_list[k].create = vertices.size() - 1;
				}
			}
			
			// for normal loop subdivision
			else
			{
				// calculate the position
                vertex.position = vertices[edge_list[k].start].position * 0.375
                                  + vertices[edge_list[k].end].position * 0.375
                                  + vertices[edge_list[k].left].position * 0.125
                                  + vertices[edge_list[k].right].position * 0.125;
				// push this vertex to the vertex array
				vertices.push_back(vertex);
				// set the index of this vertex in the vertices vector at edge
				edge_list[k].create = vertices.size() - 1;
			}
		}
	}
	
	// Loop through all the triangles to form new triangles and adjacent vertices
	for (int k=0; k<tri_ori_size; k++)
	{
		// Each triangle in here will be divided into four triangles,
		// so I select one of them to be stored in the old triangle place,
		// and add three more to append the triangles vector
		
		// code for adaptive subdivision
		if (ADAPTIVE)
		{
			// This is the most difficult part for adaptive subdivision
			// there are three cases: one new vertex, two new vertices, three new vertices
			// the case no new vertex is omitted as there are nothing to change
			
			// Case 0: no new vertex
			
			if (tri_list[k].one->create == -1
                && tri_list[k].two->create == -1
                && tri_list[k].three->create == -1)
			{
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][triangles[k].vertices[1]]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[1]][triangles[k].vertices[0]]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[2]][triangles[k].vertices[1]]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[0]]
                = true;
			}
			
			// Case 1: one new vertex
			
			// the first edge has new vertex
			else if (tri_list[k].one->create != -1
                     && tri_list[k].two->create == -1
                     && tri_list[k].three->create == -1)
			{
				// add new triangle
				MeshTriangle tri;
				tri.vertices[0] = triangles[k].vertices[0];
				tri.vertices[1] = tri_list[k].one->create;
				tri.vertices[2] = triangles[k].vertices[2];
				triangles.push_back(tri);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[0]]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[1]]
                = true;
				
				// update old triangle
				triangles[k].vertices[0] = tri_list[k].one->create;
			}
			// the second edge has new vertex
			else if (tri_list[k].one->create == -1
                     && tri_list[k].two->create != -1
                     && tri_list[k].three->create == -1)
			{
				// add new triangle
				MeshTriangle tri;
				tri.vertices[0] = triangles[k].vertices[1];
				tri.vertices[1] = tri_list[k].two->create;
				tri.vertices[2] = triangles[k].vertices[0];
				triangles.push_back(tri);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[1]]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[0]]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[0]]
                = true;
				
				// update old triangle
				triangles[k].vertices[1] = tri_list[k].two->create;
			}
			// the third edge has new vertex
			else if (tri_list[k].one->create == -1
                     && tri_list[k].two->create == -1
                     && tri_list[k].three->create != -1)
			{
				// add new triangle
				MeshTriangle tri;
				tri.vertices[0] = triangles[k].vertices[2];
				tri.vertices[1] = tri_list[k].three->create;
				tri.vertices[2] = triangles[k].vertices[1];
				triangles.push_back(tri);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[1]]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[0]]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[1]]
                = true;
				
				// update old triangle
				triangles[k].vertices[2] = tri_list[k].three->create;
			}
			
			// Case 2: two new vertices
			
			// first and second edge have new vertices
			else if (tri_list[k].one->create != -1
                     && tri_list[k].two->create != -1
                     && tri_list[k].three->create == -1)
			{
				// add new triangle
				MeshTriangle tri_1;
				tri_1.vertices[0] = triangles[k].vertices[0];
				tri_1.vertices[1] = tri_list[k].one->create;
				tri_1.vertices[2] = tri_list[k].two->create;
				triangles.push_back(tri_1);
				
				MeshTriangle tri_2;
				tri_2.vertices[0] = tri_list[k].one->create;
				tri_2.vertices[1] = triangles[k].vertices[1];
				tri_2.vertices[2] = tri_list[k].two->create;
				triangles.push_back(tri_2);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[0]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[1]][tri_list[k].two->create]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[0]]
                = true;
				
				// update old triangle
				triangles[k].vertices[1] = tri_list[k].two->create;
			}
			
			// second and third edge have new vertices
			else if (tri_list[k].one->create == -1
                     && tri_list[k].two->create != -1
                     && tri_list[k].three->create != -1)
			{
				// add new triangle
				MeshTriangle tri_1;
				tri_1.vertices[0] = triangles[k].vertices[1];
				tri_1.vertices[1] = tri_list[k].two->create;
				tri_1.vertices[2] = tri_list[k].three->create;
				triangles.push_back(tri_1);
				
				MeshTriangle tri_2;
				tri_2.vertices[0] = tri_list[k].two->create;
				tri_2.vertices[1] = triangles[k].vertices[2];
				tri_2.vertices[2] = tri_list[k].three->create;
				triangles.push_back(tri_2);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[0]][triangles[k].vertices[1]]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[1]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[0]]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[2]][tri_list[k].three->create]
                = true;
				
				// update old triangle
				triangles[k].vertices[2] = tri_list[k].three->create;
			}
			
			// third and first edge have new vertices
			else if (tri_list[k].one->create != -1
                     && tri_list[k].two->create == -1
                     && tri_list[k].three->create != -1)
			{
				// add new triangle
				MeshTriangle tri_1;
				tri_1.vertices[0] = triangles[k].vertices[2];
				tri_1.vertices[1] = tri_list[k].three->create;
				tri_1.vertices[2] = tri_list[k].one->create;
				triangles.push_back(tri_1);
				
				MeshTriangle tri_2;
				tri_2.vertices[0] = tri_list[k].three->create;
				tri_2.vertices[1] = triangles[k].vertices[0];
				tri_2.vertices[2] = tri_list[k].one->create;
				triangles.push_back(tri_2);
				
				// update adjacency table
				vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[0]][tri_list[k].one->create]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[1]][triangles[k].vertices[2]]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[2]][tri_list[k].three->create]
                = vert_adj[triangles[k].vertices[2]][triangles[k].vertices[1]]
                = true;
				
				// update old triangle
				triangles[k].vertices[0] = tri_list[k].one->create;
			}
			
			// Case 3: three new vertices
			
			else if (tri_list[k].one->create != -1
                     && tri_list[k].two->create != -1
                     && tri_list[k].three->create != -1)
			{
				// create three new triangles first (we want to use the old data)
				MeshTriangle tri_1;
				tri_1.vertices[0] = triangles[k].vertices[0];
				tri_1.vertices[1] = tri_list[k].one->create;
				tri_1.vertices[2] = tri_list[k].three->create;
				triangles.push_back(tri_1);

				MeshTriangle tri_2;
				tri_2.vertices[0] = triangles[k].vertices[1];
				tri_2.vertices[1] = tri_list[k].two->create;
				tri_2.vertices[2] = tri_list[k].one->create;
				triangles.push_back(tri_2);

				MeshTriangle tri_3;
				tri_3.vertices[0] = triangles[k].vertices[2];
				tri_3.vertices[1] = tri_list[k].three->create;
				tri_3.vertices[2] = tri_list[k].two->create;
				triangles.push_back(tri_3);

				// The vertices adjacency data can be updated here
				vert_adj[triangles[k].vertices[0]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
                = true;
				vert_adj[triangles[k].vertices[1]][tri_list[k].one->create]
                = vert_adj[triangles[k].vertices[1]][tri_list[k].two->create]
                = true;
				vert_adj[triangles[k].vertices[2]][tri_list[k].two->create]
                = vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
                = true;

				// update the existing triangle
				triangles[k].vertices[0] = tri_list[k].one->create;
				triangles[k].vertices[1] = tri_list[k].two->create;
				triangles[k].vertices[2] = tri_list[k].three->create;
			}
		}
		
		// code for loop subdivision
		else
		{
			// create three new triangles first (we want to use the old data)
			MeshTriangle tri_1;
			tri_1.vertices[0] = triangles[k].vertices[0];
			tri_1.vertices[1] = tri_list[k].one->create;
			tri_1.vertices[2] = tri_list[k].three->create;
			triangles.push_back(tri_1);

			MeshTriangle tri_2;
			tri_2.vertices[0] = triangles[k].vertices[1];
			tri_2.vertices[1] = tri_list[k].two->create;
			tri_2.vertices[2] = tri_list[k].one->create;
			triangles.push_back(tri_2);

			MeshTriangle tri_3;
			tri_3.vertices[0] = triangles[k].vertices[2];
			tri_3.vertices[1] = tri_list[k].three->create;
			tri_3.vertices[2] = tri_list[k].two->create;
			triangles.push_back(tri_3);

			// The vertices adjacency data can be updated here

			// first point
			vert_adj[triangles[k].vertices[0]][tri_list[k].one->create]
            = vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
            = true;
			// second point
			vert_adj[triangles[k].vertices[1]][tri_list[k].one->create]
            = vert_adj[triangles[k].vertices[1]][tri_list[k].two->create]
            = true;
			// third point
			vert_adj[triangles[k].vertices[2]][tri_list[k].two->create]
            = vert_adj[triangles[k].vertices[0]][tri_list[k].three->create]
            = true;

			// update the existing triangle
			triangles[k].vertices[0] = tri_list[k].one->create;
			triangles[k].vertices[1] = tri_list[k].two->create;
			triangles[k].vertices[2] = tri_list[k].three->create;
		}
	}
	
	// TEST: output the vertices, triangles and edges information
	std::cout << "New vertices = " << vertices.size() << "; "
              << "New Triangles = " << triangles.size() << std::endl;
	
	/************************************************************/
	/************** Stage 3: Move Even Vertices *****************/
	/************************************************************/
	
	/*
		In this section, I need to go through the adjacency 2D array 
		to compute the new position of the even vertices (the original 
		vertices, old part in the vertices array) 
	*/
	
	// Loop through vert_adj (x dimension)
	for (int k=0; k<vert_adj_dimen_x; k++)
	{
		// prepare the data needed after the traversal
		Vector3 sum(0.0, 0.0, 0.0);					// initialize as 0
		int N = 0;									// initialize as 0
		const float pi = 3.141592653589793238;		// the PI constant
		// Loop through vert_adj (y dimension)
		for (int n=0; n<vert_adj_dimen_y; n++)
		{
			// if n is a neighbor of k
			if (vert_adj[k][n]) {
				// add the position to the sum
				sum += vertices[n].position;
				// add the number to N
				N += 1;
			}
		}
		
		// if the number of neighbor < 2 => error
		if (N < 2) {
			// This case may trigger for adaptive subdivision
			// Do nothing when this happens
			//std::cout << "Less than 2 neighbors" << std::endl;
		}
		// if the number of neighbor == 2 => boundary case
		else if (N == 2) {
			vertices[k].position = 0.75 * vertices[k].position + 0.125 * sum;
		}
		else if (N > 2) {
			// After looping, start computing
			float Nf = (float) N;
			float beta = 1/Nf * (0.625 -(0.375+0.25*cos(2*pi/Nf))*(0.375+0.25*cos(2*pi/Nf)));
			// set the vertex to a new position
			vertices[k].position = (1 - beta*Nf)*vertices[k].position+beta*sum;
		}
	}
	
	// I guess these are the Jackpot lines to make everything show
	has_tcoords = false;
	has_normals = false;
	create_gl_data();
	
	// Memory operation
	delete [] tri_list;			// delete triangle data structure
	for (int k=0; k<vert_adj_dimen_x; k++) {
		delete [] vert_adj[k];	// delete first level of vertex adjacency
	}
	delete [] vert_adj;			// delete second level of vertex adjacency
	
    //std::cout << "Subdivision not yet implemented" << std::endl;
    return true;
}

} /* _462 */
