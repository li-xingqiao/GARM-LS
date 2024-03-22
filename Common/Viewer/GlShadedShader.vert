R"(
#version 450 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec4 aColor;

out vec4 vertColor;

layout (std140, binding = 0) uniform PassConstants
{
	mat4 uProjView;			// 0
	vec4 uViewPos;			// 64
	vec3 uAmbientStrength;	// 80
	vec3 uLightStrength;	// 96
	vec3 uLightDir;			// 112
	float uTotalTime;		// 128
	float uDeltaTime;		// 132
};

void main()
{
	gl_Position = uProjView * vec4(aPos, 1.0);
	vertColor = aColor;
}
)"
