R"(
#version 450 core

in vec4 vertColor;

out vec4 fragColor;

void main()
{
	fragColor = vec4(vertColor);
}
)"
