#include "GlApp.h"

#include "OrbitCamera.h"
#include "OrthoCamera.h"

#include <algorithm>
#include <iostream>

#include <cstdlib>

namespace PhysX {

GlApp::GlApp(const int width, const int height, const std::string &title) :
	_savedWidth(width),
	_savedHeight(height),
	_savedTitle(title)
{
	_this = this;
	// Initialize GLFW.
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA
#ifdef _DEBUG
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
#endif
	// Create window.
	glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
	glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE); // high DPI support
	_window = glfwCreateWindow(_savedWidth, _savedHeight, _savedTitle.c_str(), nullptr, nullptr);
	if (!_window) {
		std::cerr << "Error: [GLApp] Failed to create GLFW window." << std::endl;
		std::exit(-1);
	}
	// Resize components to current window size, get High-DPI working.
	glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);
	glfwGetWindowContentScale(_window, &_xScale, &_yScale);
	_scale = std::min(_xScale, _yScale);
	glfwShowWindow(_window);
	glfwMakeContextCurrent(_window);
	// Initialize GLAD.
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cerr << "Error: [GLApp] Failed to initialize GLAD." << std::endl;
		std::exit(-1);
	}
	glfwSetWindowSizeLimits(_window, _kMinimalWidth, _kMinimalHeight, GLFW_DONT_CARE, GLFW_DONT_CARE);
}

void GlApp::run()
{
	initialize();
	_timer.reset();
	while (!glfwWindowShouldClose(_window)) {
		_timer.tick();
		if (!_paused) {
			const double dt = _timer.deltaTime().count();
			processInput(dt);
			update(dt);
			clearBuffers();
			draw();
			// Swap buffers and poll IO events.
			glfwSwapBuffers(_window);
			glfwPollEvents();
		}
		else
			glfwWaitEvents();
	}
}

void GlApp::initialize()
{
	initGlStates();
	setCallbacks();
	initPrograms();
	initUniformBuffers();
	buildRenderItems();
	initCamera();

	int width;
	int height;
	glfwGetFramebufferSize(_window, &width, &height);
	resize(width, height);
}

void GlApp::resize(const int width, const int height)
{
	_width = width;
	_height = height;
	glViewport(0, 0, width, height);
	_text->resize(width, height);
	_camera->setAspectRatio(float(width) / height);
}

void GlApp::initGlStates() const
{
#ifdef _DEBUG
	glEnable(GL_DEBUG_OUTPUT);
	glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
	glDebugMessageCallback(debugMessageCallback, nullptr);
#endif
	glPointSize(4.8f);
	glLineWidth(1.5f);
}

void GlApp::setCallbacks() const
{
	glfwSetFramebufferSizeCallback(_window, framebufferSizeCallback);
	glfwSetWindowContentScaleCallback(_window, windowContentScaleCallback);
	glfwSetKeyCallback(_window, keyCallback);
	glfwSetMouseButtonCallback(_window, mouseButtonCallback);
	glfwSetCursorPosCallback(_window, cursorPosCallback);
	glfwSetScrollCallback(_window, scrollCallback);
	glfwSetWindowIconifyCallback(_window, windowIconifyCallback);
}

void GlApp::initPrograms()
{
	_programs["shaded"] = std::make_unique<GlProgram>(_kShadedVsCode, _kShadedFsCode);
	_programs["default"] = std::make_unique<GlProgram>(_kDefaultVsCode, _kDefaultFsCode);
	_programs["text"] = std::make_unique<GlProgram>(_kTextVsCode, _kTextFsCode);
}

void GlApp::initUniformBuffers()
{
	glCreateBuffers(1, &_uboPassConstants);
	glNamedBufferStorage(_uboPassConstants, sizeof(PassConstants), nullptr, GL_DYNAMIC_STORAGE_BIT);
	glBindBufferBase(GL_UNIFORM_BUFFER, 0, _uboPassConstants);
}

void GlApp::buildRenderItems()
{
	_ritems.push_back(std::make_unique<GlAxes>(_programs["shaded"].get(), _dim, 2.0f));
	_ritemLayers[size_t(RenderLayer::Opaque)].push_back(_ritems.back().get());
	_axes = dynamic_cast<GlAxes *>(_ritems.back().get());

	_ritems.push_back(std::make_unique<GlText>(_programs["text"].get()));
	_ritemLayers[size_t(RenderLayer::Text)].push_back(_ritems.back().get());
	_text = dynamic_cast<GlText *>(_ritems.back().get());
}

void GlApp::initCamera()
{
	if (_dim > 2) {
		_camera = std::make_unique<OrbitCamera>(
			0.25f * float(kPi), // fovy
			.05f, // zNear
			500.0f, // zFar
			_radius, // radius
			0.0f, // phi
			0.5f * float(kPi), // theta
			Vector3f::Zero()); // target
	}
	else {
		_camera = std::make_unique<OrthoCamera>(
			Vector3f::Unit(2) * 15.0f, // pos
			1 / _radius, // yScale
			1 / 150.0f); // zScale
	}
}

void GlApp::processInput(const double dt)
{
	if (glfwGetKey(_window, GLFW_KEY_LEFT) == GLFW_PRESS)
		_lightPhi -= float(dt);
	if (glfwGetKey(_window, GLFW_KEY_RIGHT) == GLFW_PRESS)
		_lightPhi += float(dt);
	if (glfwGetKey(_window, GLFW_KEY_UP) == GLFW_PRESS)
		_lightTheta -= float(dt);
	if (glfwGetKey(_window, GLFW_KEY_DOWN) == GLFW_PRESS)
		_lightTheta += float(dt);
	_lightPhi = std::fmod(_lightPhi, 2.0f * float(kPi));
	if (_lightPhi < 0.0f) _lightPhi += 2.0f * float(kPi);
	_lightTheta = std::clamp(_lightTheta, 0.1f, 0.5f * float(kPi));
}

void GlApp::update(const double dt)
{
	updateFrameRate();
	updateText();
	updateUniforms();
	updateGlStates();
}

void GlApp::clearBuffers() const
{
	glClearColor(_bgColor[0], _bgColor[1], _bgColor[2], 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void GlApp::draw() const
{
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	drawRenderLayer(RenderLayer::Opaque);
	glDepthMask(GL_FALSE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	drawRenderLayer(RenderLayer::Transparency);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	drawRenderLayer(RenderLayer::Text);
	glDepthMask(GL_TRUE);
}

void GlApp::updateFrameRate()
{
	static uint framesCnt = 0;
	static double lastTotTime = _timer.totalTime().count();
	double currTotTime = _timer.totalTime().count();

	framesCnt++;
	if (currTotTime - lastTotTime >= 1.0) {
		_framesPerSecond = uint(framesCnt / (currTotTime - lastTotTime) + 0.5);
		framesCnt = 0;
		lastTotTime = currTotTime;
	}
}

void GlApp::updateText()
{
	_text->reset();
	if (_enableInfo) {
		_text->set(
			std::format(
				"  Rendering info (F1):  {}\n"
				"Gamma correction (F2):  {}\n"
				"   Multisampling (F3):  {}\n"
				"  Wireframe mode (F4):  {}",
				_enableInfo ? "on" : "off",
				_enableSrgb ? "on" : "off",
				_enableMsaa ? "on" : "off",
				_enableWireframe ? "on" : "off"),
			Vector2f(10.24f / _width * _scale, 7.68f / _height * _scale),
			Vector2f(0.75f * _scale, 0.75f * _scale),
			Vector4f(0, 0, 0, 1));
		_text->set(
			std::format("FPS: {:>3}", _framesPerSecond),
			Vector2f(1.0f - 10.24f / _width * _scale, 7.68f / _height * _scale),
			Vector2f(0.75f * _scale, 0.75f * _scale),
			Vector4f(0, 0, 0, 1),
			GlText::Alignment::Right);
		_text->set(
			std::format("Light direction: ({:.2f}, {:.2f})", _lightPhi, _lightTheta),
			Vector2f(10.24f / _width * _scale, 1.0f - (24.0f + 7.68f) / _height * _scale),
			Vector2f(0.75f * _scale, 0.75f * _scale),
			Vector4f(0, 0, 0, 1));
	}
}

void GlApp::updateUniforms()
{
	_camera->update();
	// Construct pass constants.
	_passConstants.projView = _camera->projView();
	_passConstants.viewPos = _camera->pos();
	_passConstants.ambientStrength = Vector3f(0.25f, 0.25f, 0.35f);
	_passConstants.lightStrength = Vector3f(1.0f, 1.0f, 0.9f);
	_passConstants.lightDir = OrbitCamera::sphericalToCartesian(1.0f, _lightPhi, _lightTheta);
	_passConstants.totalTime = float(_timer.totalTime().count());
	_passConstants.deltaTime = float(_timer.deltaTime().count());
	// Upload to Graphics Memory.
	glNamedBufferSubData(_uboPassConstants, 0, sizeof(_passConstants), &_passConstants);
}

void GlApp::framebufferSizeCallback(GLFWwindow *window, int width, int height)
{
	_this->resize(width, height);
}

void GlApp::windowContentScaleCallback(GLFWwindow *window, float xScale, float yScale)
{
	_this->_xScale = xScale;
	_this->_yScale = yScale;
	_this->_scale = std::min(xScale, yScale);
}

void GlApp::keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch (key) {
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, true);
			break;
		case GLFW_KEY_SPACE:
			glfwSetWindowSize(window, int(_this->_savedWidth * _this->_xScale), int(_this->_savedHeight * _this->_yScale));
			_this->_camera->reset();
			_this->_lightPhi = _kSavedLightPhi;
			_this->_lightTheta = _kSavedLightTheta;
			break;
		case GLFW_KEY_F1:
			_this->_enableInfo = !_this->_enableInfo;
			break;
		case GLFW_KEY_F2:
			_this->_enableSrgb = !_this->_enableSrgb;
			break;
		case GLFW_KEY_F3:
			_this->_enableMsaa = !_this->_enableMsaa;
			break;
		case GLFW_KEY_F4:
			_this->_enableWireframe = !_this->_enableWireframe;
			break;
		case GLFW_KEY_0:
			_this->_axes->flipVisibility();
			break;
		}
	}
}

void GlApp::mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
	if (action == GLFW_PRESS) {
		glfwGetCursorPos(window, &_this->_lastMousePos.x(), &_this->_lastMousePos.y());
	}
}

void GlApp::cursorPosCallback(GLFWwindow *window, double xpos, double ypos)
{
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		_this->_camera->rotate(float(xpos - _this->_lastMousePos.x()), float(ypos - _this->_lastMousePos.y()));
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
		_this->_camera->translate(float(xpos - _this->_lastMousePos.x()), float(ypos - _this->_lastMousePos.y()));
	}
	_this->_lastMousePos << xpos, ypos;
}

void GlApp::scrollCallback(GLFWwindow *window, double xoffset, double yoffset)
{
	_this->_camera->scale(float(yoffset));
}

void GlApp::windowIconifyCallback(GLFWwindow *window, int iconified)
{
	if (iconified) {
		_this->_paused = true;
		_this->_timer.stop();
	}
	else {
		_this->_paused = false;
		_this->_timer.start();
	}
}

void APIENTRY GlApp::debugMessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *message, const void *userParam)
{
	std::cerr << "Debug message (" << id << "): " << message << std::endl;

	switch (source) {
	case GL_DEBUG_SOURCE_API:				std::cerr << "Source: API";					break;
	case GL_DEBUG_SOURCE_WINDOW_SYSTEM:		std::cerr << "Source: Window System";		break;
	case GL_DEBUG_SOURCE_SHADER_COMPILER:	std::cerr << "Source: Shader Compiler";		break;
	case GL_DEBUG_SOURCE_THIRD_PARTY:		std::cerr << "Source: Third Party";			break;
	case GL_DEBUG_SOURCE_APPLICATION:		std::cerr << "Source: Application";			break;
	case GL_DEBUG_SOURCE_OTHER:				std::cerr << "Source: Other";				break;
	}
	std::cerr << std::endl;

	switch (type) {
	case GL_DEBUG_TYPE_ERROR:				std::cerr << "Type: Error";					break;
	case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:	std::cerr << "Type: Deprecated Behaviour";	break;
	case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:	std::cerr << "Type: Undefined Behaviour";	break;
	case GL_DEBUG_TYPE_PORTABILITY:			std::cerr << "Type: Portability";			break;
	case GL_DEBUG_TYPE_PERFORMANCE:			std::cerr << "Type: Performance";			break;
	case GL_DEBUG_TYPE_MARKER:				std::cerr << "Type: Marker";				break;
	case GL_DEBUG_TYPE_PUSH_GROUP:			std::cerr << "Type: Push Group";			break;
	case GL_DEBUG_TYPE_POP_GROUP:			std::cerr << "Type: Pop Group";				break;
	case GL_DEBUG_TYPE_OTHER:				std::cerr << "Type: Other";					break;
	}
	std::cerr << std::endl;

	switch (severity) {
	case GL_DEBUG_SEVERITY_HIGH:			std::cerr << "Severity: high";				break;
	case GL_DEBUG_SEVERITY_MEDIUM:			std::cerr << "Severity: medium";			break;
	case GL_DEBUG_SEVERITY_LOW:				std::cerr << "Severity: low";				break;
	case GL_DEBUG_SEVERITY_NOTIFICATION:	std::cerr << "Severity: notification";		break;
	}
	std::cerr << std::endl << std::endl;
}

}
