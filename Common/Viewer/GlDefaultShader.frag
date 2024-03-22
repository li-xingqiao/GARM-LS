R"(
#version 450 core

in vec3 vertPos;
in vec3 vertNormal;
in float vertHeat;

out vec4 fragColor;

uniform vec4 uDiffuseAlbedo;
uniform vec3 uFresnelR0;
uniform float uRoughness;
uniform uint uEnableColorMap;

vec3 diffuseColor;

layout (std140, binding = 0) uniform PassConstants
{
	mat4 uProjView;			// 0
	vec3 uViewPos;			// 64
	vec3 uAmbientStrength;	// 80
	vec3 uLightStrength;	// 96
	vec3 uLightDir;			// 112
	float uTotalTime;		// 128
	float uDeltaTime;		// 132
};

vec3 calcSchlickFresnel(const vec3 R0, const vec3 normal, const vec3 lightDir)
{
	float cosIncidentAngle = clamp(dot(normal, lightDir), 0.0, 1.0);
	float f0 = 1.0 - cosIncidentAngle;
	return R0 + (1.0 - R0) * (f0 * f0 * f0 * f0 * f0);
}

vec3 calcDiffuseAndSpecular(const vec3 lightStrength, const vec3 lightDir, const vec3 normal, const vec3 viewDir)
{
	vec3 halfDir = normalize(lightDir + viewDir);
	float m = (1.0 - uRoughness) * 256.0;
	float roughnessFactor = (m + 8.0) * pow(max(dot(halfDir, normal), 0.0), m) / 8.0;
	vec3 fresnelFactor = calcSchlickFresnel(uFresnelR0, halfDir, lightDir);
	vec3 specularAlbedo = fresnelFactor * roughnessFactor;
	// Apply LDR.
	specularAlbedo = specularAlbedo / (specularAlbedo + 1.0);
	return (diffuseColor + specularAlbedo) * lightStrength;
}

vec3 ComputeDirectionalLight(vec3 normal, const vec3 viewDir)
{
	if (dot(normal, viewDir) < 0) normal = -normal;
	vec3 lightStrength = uLightStrength * max(dot(uLightDir, normal), 0.0);
	return calcDiffuseAndSpecular(lightStrength, uLightDir, normal, viewDir);
}

float jetBase(const float val)
{
	if (val <= 0.125) return 0.0;
	else if (val <= 0.375) return 4.0 * val - 0.5;
	else if (val <= 0.625) return 1.0;
	else if (val <= 0.875) return 3.5 - 4.0 * val;
	else return 0.0;
}

vec3 getJet(const float heat)
{
	return vec3(jetBase(heat - 0.25), jetBase(heat), jetBase(heat + 0.25));
}

void main()
{
	diffuseColor = uEnableColorMap == 0 ? uDiffuseAlbedo.rgb : getJet(vertHeat);
	vec3 ambient = uAmbientStrength * diffuseColor;
	vec3 diffspec = ComputeDirectionalLight(normalize(vertNormal), normalize(uViewPos - vertPos));
	fragColor = vec4(ambient + diffspec, uDiffuseAlbedo.a);
}
)"
