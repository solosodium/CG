varying vec3 normal;

void main(void)
{
    vec3 n_normal = normalize( normal );

    // convert normal values from domain [-1.0 - 1.0] to [0.0, 1.0]
    // y = 0.5x + 0.5
    vec3 c_normal = 0.5 * ( n_normal + vec3( 1.0, 1.0, 1.0 ) );

    gl_FragColor = vec4( c_normal, 0.0 );
}

