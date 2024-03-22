#pragma once

#include "Camera.h"
#include "Constants.h"
#include "GlAxes.h"
#include "GlText.h"
#include "StepTimer.h"
#include "Types.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <memory>
#include <numbers>
#include <string>
#include <unordered_map>
#include <vector>

namespace PhysX {

class GlApp
{
protected:

	enum class RenderLayer : uint { Opaque, Transparency, Text, Count };

	struct PassConstants
	{
		Matrix4f projView;			// 0
		Vector3f viewPos;			// 64
		float _pad1 = 0.0f;
		Vector3f ambientStrength;	// 80
		float _pad2 = 0.0f;
		Vector3f lightStrength;		// 96
		float _pad3 = 0.0f;
		Vector3f lightDir;			// 112
		float totalTime;			// 124
		float deltaTime;			// 128
	};

	static inline GlApp *_this = nullptr;

	static constexpr int _kMinimalWidth = 640;
	static constexpr int _kMinimalHeight = 480;
	static constexpr float _kSavedLightPhi = 0.0f;
	static constexpr float _kSavedLightTheta = 0.25f * float(kPi);

	static constexpr GLchar _kShadedVsCode[] =
#include "GlShadedShader.vert"
		;
	static constexpr GLchar _kShadedFsCode[] =
#include "GlShadedShader.frag"
		;
	static constexpr GLchar _kDefaultVsCode[] =
#include "GlDefaultShader.vert"
		;
	static constexpr GLchar _kDefaultFsCode[] =
#include "GlDefaultShader.frag"
		;
	static constexpr GLchar _kTextVsCode[] =
#include "GlTextShader.vert"
		;
	static constexpr GLchar _kTextFsCode[] =
#include "GlTextShader.frag"
		;

	GLuint _uboPassConstants = 0;

	bool _paused = false;
	bool _enableInfo = false;
	bool _enableSrgb = false;
	bool _enableMsaa = false;
	bool _enableWireframe = false;

	GLFWwindow *_window;

	const int _savedWidth;
	const int _savedHeight;
	float _xScale = 1;
	float _yScale = 1;
	float _scale = 1;
	const std::string _savedTitle;


	int _width = 0;
	int _height = 0;
	Vector3f _bgColor = Vector3f(176, 196, 222) / 255;
	Vector2d _lastMousePos = Vector2d::Zero();

	int _dim = 3;
	float _radius = 5;
	std::unique_ptr<Camera> _camera;
	float _lightPhi = _kSavedLightPhi;
	float _lightTheta = _kSavedLightTheta;

	StepTimer _timer;
	uint _framesPerSecond = 0;

	std::unordered_map<std::string, std::unique_ptr<GlProgram>> _programs;

	std::vector<std::unique_ptr<GlRenderItem>> _ritems;
	std::vector<GlRenderItem *> _ritemLayers[size_t(RenderLayer::Count)];
	GlAxes *_axes = nullptr;
	GlText *_text = nullptr;

	PassConstants _passConstants;

public:

	GlApp(const int width, const int height, const std::string &caption);

	GlApp(const GlApp &rhs) = delete;
	GlApp &operator=(const GlApp &rhs) = delete;

	virtual ~GlApp()
	{
		glDeleteBuffers(1, &_uboPassConstants);
		glfwTerminate();
	}

	void run();

protected:

	void initialize();

	virtual void resize(const int width, const int height);
	virtual void initGlStates() const;
	virtual void setCallbacks() const;
	virtual void initPrograms();
	virtual void initUniformBuffers();
	virtual void buildRenderItems();
	virtual void initCamera();

	virtual void processInput(const double dt);
	virtual void update(const double dt);
	virtual void clearBuffers() const;
	virtual void draw() const;

	void drawRenderLayer(const RenderLayer layer) const
	{
		for (auto ritem : _ritemLayers[size_t(layer)])
			if (ritem->isVisible()) ritem->draw();
	}

	void updateFrameRate();
	virtual void updateText();
	virtual void updateUniforms();

	virtual void updateGlStates()
	{
		_enableSrgb ? glEnable(GL_FRAMEBUFFER_SRGB) : glDisable(GL_FRAMEBUFFER_SRGB);
		_enableMsaa ? glEnable(GL_MULTISAMPLE) : glDisable(GL_MULTISAMPLE);
		_enableWireframe ? glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) : glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	static void framebufferSizeCallback(GLFWwindow *window, int width, int height);
	static void windowContentScaleCallback(GLFWwindow *window, float xScale, float yScale);
	static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
	static void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods);
	static void cursorPosCallback(GLFWwindow *window, double xpos, double ypos);
	static void scrollCallback(GLFWwindow *window, double xoffset, double yoffset);
	static void windowIconifyCallback(GLFWwindow *window, int iconified);

	static void APIENTRY debugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, const void *userParam);
};

}
