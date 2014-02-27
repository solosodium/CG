varying vec3 norm;
varying vec3 cam_dir;
varying vec3 color;

// TODO: declare the cubemap variable [samplerCube]
uniform samplerCube cubemap;

void main(void)
{
	norm = gl_NormalMatrix * gl_Normal;
	cam_dir = (gl_ModelViewMatrix*gl_Vertex).xyz;
	gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;
	vec4 col = gl_Color;
	if (all(equal(gl_Color, vec4(1.0)))) col = vec4(0.0);
	gl_FrontColor = col;
	gl_BackColor = col;
}
