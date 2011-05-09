uniform sampler2D tone0;
uniform sampler2D tone1;
uniform sampler2D tone2;
uniform sampler2D tone3;
uniform sampler2D tone4;
uniform sampler2D tone5;
uniform sampler2D alph;
varying float lightIntensity;
varying vec4 _origTexCoord;

void main()
{
	gl_Position = ftransform();	
	//gl_Normal*=gl_NormalMatrix;
	//normalize(gl_Normal);
	vec3 lightoff = gl_ModelViewMatrix*gl_Vertex - gl_LightSource[0].position;
	lightIntensity = dot(normalize(lightoff),gl_Normal);
	gl_TexCoord[0] =gl_TextureMatrix[0] * gl_MultiTexCoord0;
	_origTexCoord = gl_MultiTexCoord0;
}
