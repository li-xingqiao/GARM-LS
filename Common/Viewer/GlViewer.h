#pragma once

#include "GlSimulated.h"

#include "GlApp.h"
#include "Yaml.h"

#include <iostream>

#include <cstdio>

namespace PhysX {

class GlViewer : public GlApp
{
protected:

	static inline GlViewer *_this = nullptr;

	const std::string _outputDir;
	const std::string _recordingCommand;
	std::string _recordingDir;
	YAML::Node _root;

	uint _endFrame;
	uint _frameRate;

	bool _enableFrameNum = true;
	bool _enableRecording = false;
	bool _playing = false;
	double _currentFrame = 0;

	std::vector<GlSimulated *> _simulatedObjects;

public:

	GlViewer(const std::string &outputDir, const uint frameRate);

	GlViewer(const GlViewer &rhs) = delete;
	GlViewer &operator=(const GlViewer &rhs) = delete;
	virtual ~GlViewer() = default;

protected:

	virtual void setCallbacks() const override;
	virtual void buildRenderItems() override;

	virtual void update(const double dt) override;
	virtual void updateText() override;

	void capture(const std::string &fileName) const;
	void generateVideo() const { std::system(std::vformat(_recordingCommand, std::make_format_args(_recordingDir)).c_str()); }

	static void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
};

}
