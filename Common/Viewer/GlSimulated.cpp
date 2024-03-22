#include "GlSimulated.h"

#include "IO.h"

#include <format>
#include <iostream>

#include <cstdlib>

namespace PhysX {
GlSimulated::GlSimulated(GlProgram *const program, const std::string &outputDir, const uint endFrame, const int dim, const YAML::Node &node) :
	GlRenderItem(program)
{
	// Read name.
	if (!node["name"]) reportError("unnamed object");
	const std::string name = node["name"].as<std::string>();

	// Read data mode.
	if (!node["data_mode"]) reportError("missing data mode");
	const std::string dataMode = node["data_mode"].as<std::string>();

	// Read primitive type.
	if (!node["primitive_type"]) reportError("missing primitive type");
	const std::string primitiveType = node["primitive_type"].as<std::string>();
	if (strToMode.find(primitiveType) == strToMode.end()) reportError("invalid primive type");
	_mode = strToMode[primitiveType];

	// Read indexed.
	if (node["indexed"]) _indexed = node["indexed"].as<bool>();
	else _indexed = false;

	// Read color map.
	if (node["color_map"] && node["color_map"]["enabled"])
		_enableColorMap = node["color_map"]["enabled"].as<bool>();
	else _enableColorMap = false;
	bool normalizedHeat = false;
	if (node["color_map"] && node["color_map"]["normalized"])
		normalizedHeat = node["color_map"]["normalized"].as<bool>();

	// Read material.
	_diffuseAlbedo = Vector4f(0.5f, 0.5f, 0.5f, 1.0f);
	_fresnelR0 = Vector3f(0.02041f, 0.02041f, 0.02041f);
	_roughness = dim > 2 ? 0.75f : 1.03125f;
	if (node["material"]) {
		if (node["material"]["diffuse_albedo"])
			_diffuseAlbedo = node["material"]["diffuse_albedo"].as<Vector4f>();
		if (node["material"]["fresnel_r0"])
			_fresnelR0 = node["material"]["fresnel_r0"].as<Vector3f>();
		if (node["material"]["roughness"])
			_roughness = node["material"]["roughness"].as<float>();
	}

	// Initialize meshes.
	std::vector<Vector3f> positions;
	std::vector<Vector3f> normals;
	std::vector<float> heats;
	std::vector<uint> indices;
	if (dataMode == "static") {
		_vtxFrameOffset.push_back(0);
		if (_indexed) _idxFrameOffset.push_back(0);
		loadMesh(
			std::format("{}/frames/0/{}.mesh", outputDir, name),
			dim,
			&positions,
			&normals,
			_enableColorMap ? &heats : nullptr,
			_indexed ? &indices : nullptr);
	}
	else if (dataMode == "dynamic") {
		_vtxFrameOffset.push_back(0);
		if (_indexed) _idxFrameOffset.push_back(0);
		for (uint frame = 0; frame < endFrame; frame++) {
			loadMesh(
				std::format("{}/frames/{}/{}.mesh", outputDir, frame, name),
				dim,
				&positions,
				&normals,
				_enableColorMap ? &heats : nullptr,
				_indexed ? &indices : nullptr);
		}
	}
	else if (dataMode == "semi-dynamic") {
		_vtxFrameOffset.push_back(0);
		if (_indexed) _idxFrameOffset.push_back(0);
		for (uint frame = 0; frame < endFrame; frame++) {
			loadMesh(
				std::format("{}/frames/{}/{}.mesh", outputDir, frame, name),
				dim,
				&positions,
				&normals,
				_enableColorMap ? &heats : nullptr,
				_indexed && !frame ? &indices : nullptr);
		}
	}
	else reportError("invalid data mode");

	for (size_t frame = 1; frame < _vtxFrameOffset.size(); frame++) {
		_vtxFrameOffset[frame] += _vtxFrameOffset[frame - 1];
	}
	if (_indexed) {
		for (size_t frame = 1; frame < _idxFrameOffset.size(); frame++) {
			_idxFrameOffset[frame] += _idxFrameOffset[frame - 1];
		}
	}

	// Normalized heats.
	if (_enableColorMap && !normalizedHeat) {
		const auto minmax = std::minmax_element(heats.begin(), heats.end());
		float minimum = *minmax.first;
		float maximum = *minmax.second;
		if (minimum == maximum) minimum -= 1, maximum += 1;
		for (auto &x : heats) x = (x - minimum) / (maximum - minimum);
	}

	// Create vertex buffers.
	auto sizeBase = positions.size() * sizeof(float);
	glCreateBuffers(1, &_vbo);
	glNamedBufferData(_vbo, (_enableColorMap ? 7 : 6) * sizeBase, nullptr, GL_STATIC_DRAW);		
	// Attribute 0.
	glNamedBufferSubData(_vbo, 0, 3 * sizeBase, positions.data());
	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 3 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 0, 0);
	glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glEnableVertexArrayAttrib(_vao, 0);
	// Attribute 1.
	glNamedBufferSubData(_vbo, 3 * sizeBase, 3 * sizeBase, normals.data());
	glVertexArrayVertexBuffer(_vao, 1, _vbo, 3 * sizeBase, 3 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 1, 1);
	glVertexArrayAttribFormat(_vao, 1, 3, GL_FLOAT, GL_FALSE, 0);
	glEnableVertexArrayAttrib(_vao, 1);
	// Attribute 2.
	if (_enableColorMap) {
		glNamedBufferSubData(_vbo, 6 * sizeBase, sizeBase, heats.data());
		glVertexArrayVertexBuffer(_vao, 2, _vbo, 6 * sizeBase, sizeof(float));
		glVertexArrayAttribBinding(_vao, 2, 2);
		glVertexArrayAttribFormat(_vao, 2, 1, GL_FLOAT, GL_FALSE, 0);
		glEnableVertexArrayAttrib(_vao, 2);
	}

	// Create index buffers.
	if (_indexed) {
		glCreateBuffers(1, &_ebo);
		glNamedBufferStorage(_ebo, indices.size() * sizeof(indices[0]), indices.data(), 0);
		glVertexArrayElementBuffer(_vao, _ebo);
	}
}

void GlSimulated::beginDraw()
{
	// Set parameters according to current frame.
	if (_indexed) {
		const size_t vtxFrame = _currentFrame < _vtxFrameOffset.size() - 1 ? _currentFrame : 0;
		const size_t idxFrame = _currentFrame < _idxFrameOffset.size() - 1 ? _currentFrame : 0;
		_count = _idxFrameOffset[idxFrame + 1] - _idxFrameOffset[idxFrame];
		_indices = reinterpret_cast<const void *>(_idxFrameOffset[idxFrame] * sizeof(uint));
		_baseVertex = _vtxFrameOffset[vtxFrame];
	}
	else {
		const size_t vtxFrame = _currentFrame < _vtxFrameOffset.size() - 1 ? _currentFrame : 0;
		_first = _vtxFrameOffset[vtxFrame];
		_count = _vtxFrameOffset[vtxFrame + 1] - _vtxFrameOffset[vtxFrame];
	}
	// Set uniforms.
	_program->setUniform("uDiffuseAlbedo", _diffuseAlbedo);
	_program->setUniform("uFresnelR0", _fresnelR0);
	_program->setUniform("uRoughness", _roughness);
	_program->setUniform("uEnableColorMap", uint(_enableColorMap));
}

void GlSimulated::loadMesh(
	const std::string &fileName,
	const int dim,
	std::vector<Vector3f> *const positions,
	std::vector<Vector3f> *const normals,
	std::vector<float> *const heats,
	std::vector<uint> *const indices)
{
	std::ifstream fin(fileName, std::ios::binary);
	if (!fin) {
		std::cerr << std::format("Error: [GlSimulated] failed to open {}.", fileName) << std::endl;
		std::exit(-1);
	}
	uint vtxCnt;
	IO::readValue(fin, vtxCnt);
	_vtxFrameOffset.push_back(vtxCnt);
	if (positions) {
		positions->resize(positions->size() + vtxCnt, Vector3f::Zero().eval());
		if (dim > 2)
			IO::readArray(fin, positions->data() + positions->size() - vtxCnt, vtxCnt);
		else {
			for (uint i = 0; i < vtxCnt; i++) {
				IO::read(fin, positions->data() + positions->size() - vtxCnt + i, sizeof(Vector2f));
			}
		}
	}
	if (normals) {
		normals->resize(normals->size() + vtxCnt, Vector3f::Zero().eval());
		if (dim > 2)
			IO::readArray(fin, normals->data() + normals->size() - vtxCnt, vtxCnt);
		else {
			for (uint i = 0; i < vtxCnt; i++)
				normals->operator[](normals->size() - vtxCnt + i).z() = 1.0f;
		}
	}
	if (heats) {
		heats->resize(heats->size() + vtxCnt);
		IO::readArray(fin, heats->data() + heats->size() - vtxCnt, vtxCnt);
	}
	if (indices) {
		uint idxCnt;
		IO::readValue(fin, idxCnt);
		_idxFrameOffset.push_back(idxCnt);
		indices->resize(indices->size() + idxCnt);
		IO::readArray(fin, indices->data() + indices->size() - idxCnt, idxCnt);
	}
}

void GlSimulated::reportError(const std::string &msg) const
{
	std::cerr << std::format("Error: [GlSimulated] encountered {}.", msg) << std::endl;
	std::exit(-1);
}

}
