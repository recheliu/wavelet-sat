"/*\n"
"\n"
"This is the shader program for direct volume rendering\n"
"\n"
"*/\n"
"	uniform sampler3D t3dVolume;	// the texture hold the depths of each knots\n"
"	uniform sampler2D t2dPrevLayer;	// the texture hold the depths of each knots\n"
"	uniform sampler2D t2dsDepth;	// the texture hold the depths of each knots\n"
"	uniform float fThicknessGain;\n"
"	uniform float fWindowWidth;\n"
"	uniform float fWindowHeight;\n"
"\n"
"	#include \"clip_frag_func.frag.h\"\n"
"	#include \"tf1d_frag_func.frag.h\"\n"
"\n"
"void \n"
"main()\n"
"{\n"
"	vec4 v4FragCoord = gl_FragCoord;\n"
"	v4FragCoord.x /= fWindowWidth;\n"
"	v4FragCoord.y /= fWindowHeight;\n"
"	v4FragCoord.w = 1.0 / v4FragCoord.w;\n"
"\n"
"				// read the previous layer\n"
"	vec4 v4PrevFragData = texture2D(t2dPrevLayer, v4FragCoord.xy);\n"
"	float fPrevDepth = v4PrevFragData.r;\n"
"	float fPrevV_normalized = v4PrevFragData.g;\n"
"\n"
"						// calcualte the distance between this fragment and the previous fragment in the object space \n"
"	vec4 v4PrevCoord_ndc = v4FragCoord;\n"
"	v4PrevCoord_ndc.z = fPrevDepth;\n"
"	v4PrevCoord_ndc.xyz = v4PrevCoord_ndc.xyz * 2.0 - 1.0;\n"
"	v4PrevCoord_ndc.w = 1.0;\n"
"	vec4 v4PrevCoord_obj = gl_ModelViewProjectionMatrixInverse * v4PrevCoord_ndc;\n"
"	v4PrevCoord_obj = v4PrevCoord_obj / v4PrevCoord_obj.w;\n"
"	vec4 v4PrevCoord_eye = gl_ModelViewMatrix * v4PrevCoord_obj;\n"
"\n"
"	vec4 v4Coord_ndc = v4FragCoord;\n"
"	v4Coord_ndc.xyz = v4Coord_ndc.xyz * 2.0 - 1.0; \n"
"	v4Coord_ndc.w = 1.0;\n"
"	vec4 v4Coord_obj = gl_ModelViewProjectionMatrixInverse * v4Coord_ndc;\n"
"	v4Coord_obj = v4Coord_obj / v4Coord_obj.w;\n"
"	vec4 v4Coord_eye = gl_ModelViewMatrix * v4Coord_obj;\n"
"\n"
"	float fThickness_obj;\n"
"	fThickness_obj = distance(v4Coord_obj, v4PrevCoord_obj);\n"
"\n"
"	fThickness_obj *= fThicknessGain;\n"
"\n"
"	float fThicknessRatio = FAdjustThickness(v4FragCoord.z, fPrevDepth, v4FragCoord.xy);\n"
"	fThickness_obj *= fThicknessRatio;\n"
"\n"
"				// read the value \n"
"	float fV_normalized = texture3D(t3dVolume, gl_TexCoord[0].xyz).x;\n"
"\n"
"				// convert the value into color via the transfer function\n"
"	vec4 v4Color = F4GetColorFrom1DTf(fV_normalized);\n"
"\n"
"	v4Color.a = 1.0 - exp(-v4Color.a * fThickness_obj);\n"
"\n"
"	float fBackgroundDepth = texture2D(t2dsDepth, v4FragCoord.xy).r;\n"
"	if( v4FragCoord.z > fBackgroundDepth )\n"
"		v4Color.a = 0.0;\n"
"\n"
"	vec4 v4Data = vec4(\n"
"		v4FragCoord.z, \n"
"		fV_normalized,\n"
"		0.0, \n"
"		1.0);\n"
"\n"
"	gl_FragData[0] = v4Color;\n"
"	gl_FragData[1] = v4Data;\n"
"}\n"
"\n"
"/*\n"
"\n"
"$Log: not supported by cvs2svn $\n"
"Revision 1.3  2010/09/10 14:08:40  leeten\n"
"\n"
"[09/10/2010]\n"
"1. [MOD] Treat the sampler t2dsDepth as a regular 2D texture.\n"
"\n"
"Revision 1.2  2010/05/15 16:06:23  leeten\n"
"\n"
"[05/15/2010]\n"
"1. [MOD] Change the ray integral equation from A = 1 - (1 - A)^D to A = 1 - exp(-AD).\n"
"\n"
"Revision 1.1  2010/04/06 20:01:15  leeten\n"
"\n"
"[04/06/2010]\n"
"1. [1ST] First time checkin.\n"
"\n"
"\n"
"*/\n"
""