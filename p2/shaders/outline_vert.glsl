// Varying parameters
varying float value_l;      // dot product of light and normal
varying float value_c;      // dot product of camera and normal

void main(void)
{
    // get normal vector
	vec3 normal = gl_NormalMatrix * gl_Normal;
	
    // get camera direction
    vec3 camera = (gl_ModelViewMatrix*gl_Vertex).xyz;
    
    // get light direction
    vec3 light = gl_LightSource[0].halfVector.xyz;
    
    // compute the dot product value of normal and light
    value_l = dot(normalize(normal), normalize(light));
    
    // compute the dot product value of normal and camera
    value_c = dot(normalize(normal), normalize(camera));
    
    // compute the position
    gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;
}
