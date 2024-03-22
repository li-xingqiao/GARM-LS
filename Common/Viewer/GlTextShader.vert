R"(
#version 450 core

layout (location = 0) in vec4 aPos;
layout (location = 1) in vec4 aTexCoord;
layout (location = 2) in vec4 aColor;

out vec4 vertColor;
out vec2 vertTexCoord;

void main()
{
	vec2 uv = vec2(gl_VertexID & 1, (gl_VertexID >> 1) &1); 
	gl_Position = vec4(aPos.x + (aPos.z * uv.x), aPos.y - (aPos.w * uv.y), 0, 1);
	vertColor = aColor;
	vertTexCoord = vec2(aTexCoord.x + (aTexCoord.z * uv.x), aTexCoord.y + (aTexCoord.w * uv.y));
}
)"
