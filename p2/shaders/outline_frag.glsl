// Varying parameters
varying float value_l;      // dot product of light and normal
varying float value_c;      // dot product of camera and normal

void main(void)
{
    // if the dot product of camera and normal close to zero
    // it should have the outliner shading effect
    // the result will look like a toon shader
    
    // the outline width
    float width = 0.4;
    
    // check the close to zero condition
    if (value_c < width && value_c > -width) {
        // the outline dark edge should change with the light
        gl_FragColor = vec4(0.1, 0.1, 0.1, 1.0) * value_l;
    }
    // if not close to zero
    else {
        // should set as pure color (make toon effect)
        gl_FragColor = vec4(1.0, 1.0, 0.0, 1.0);
    }
}
