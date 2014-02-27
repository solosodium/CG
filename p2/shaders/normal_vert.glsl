varying vec3 normal;

void main(void)
{
    gl_Position = ftransform();
    normal = normalize( gl_Normal );
}
