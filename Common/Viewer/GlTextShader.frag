R"(
#version 450 core

in vec4 vertColor;
in vec2 vertTexCoord;

out vec4 fragColor;

uniform sampler2D uText;

void main()
{
	fragColor = vec4(vertColor.rgb, vertColor.a * texture(uText, vertTexCoord).r);
}
)"
