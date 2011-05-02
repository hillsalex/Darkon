uniform sampler2D tone0;
uniform sampler2D tone1;
uniform sampler2D tone2;
uniform sampler2D tone3;
uniform sampler2D tone4;
uniform sampler2D tone5;
varying float lightIntensity;

void main()
{
	vec4 color_0 = texture2DProj(tone0,gl_TexCoord[0]);
	vec4 color_1 = texture2DProj(tone1,gl_TexCoord[0]);
	vec4 color_2 = texture2DProj(tone2,gl_TexCoord[0]);
	vec4 color_3 = texture2DProj(tone3,gl_TexCoord[0]);
	vec4 color_4 = texture2DProj(tone4,gl_TexCoord[0]);
	vec4 color_5 = texture2DProj(tone5,gl_TexCoord[0]);
	
	float p0 = 0.0;
	float p1 = 0.08;
	float p2 = 0.15;
	float p3 = 0.3;
	float p4 = 0.4;
	float p5 = 0.8;
	float p6 = 1.0;
	


	if(lightIntensity<p1)
		gl_FragColor = mix(color_5,color_4,lightIntensity/p1);
	else if(lightIntensity<p2)
		gl_FragColor = mix(color_4,color_3,(lightIntensity-p1)/(p2-p1));
	else if(lightIntensity<p3)
		gl_FragColor = mix(color_3,color_2,(lightIntensity-p2)/(p3-p2));
	else if(lightIntensity<p4)
		gl_FragColor = mix(color_2,color_1,(lightIntensity-p3)/(p4-p3));
	else if(lightIntensity<p5)
		gl_FragColor = mix(color_1,color_0,(lightIntensity-p4)/(p5-p4));
	else if (lightIntensity<p6)
		gl_FragColor = mix(color_0,vec4(1.0,1.0,1.0,1.0),(lightIntensity-p5)/(p6-p5));
	else
		gl_FragColor = vec4(1.0,1.0,1.0,1.0);

	/*if(8*gl_FragColor.x - 3.5 <= 0)
		gl_FragColor = vec4(0.0,0.0,0.0,0.0);
	else
		gl_FragColor = vec4(1.0,1.0,1.0,1.0);*/


}
