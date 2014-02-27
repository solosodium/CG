// set by the vertex shader
varying vec2 tex_coord;
// set by the C program
uniform float green;

void main(void)
{
    // render a color based on the texture coordinate
    gl_FragColor = vec4( tex_coord.s, 0.8, tex_coord.t, 1.0 );
}

