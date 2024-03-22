#pragma once

#include "Types.h"

#include <glad/glad.h>

#include <string>

namespace PhysX {

class GlProgram final
{
protected:

	GLuint _program;

public:

	GlProgram(const GLchar *vsCode, const GLchar *fsCode, const GLint vsLength = -1, const GLint fsLength = -1);

	GlProgram(const GlProgram &rhs) = delete;
	GlProgram &operator=(const GlProgram &rhs) = delete;
	virtual ~GlProgram() = default;

	void use() const { glUseProgram(_program); }

	void setUniform(const std::string &name, const GLfloat value) const { glProgramUniform1f(_program, glGetUniformLocation(_program, name.c_str()), value); }
	void setUniform(const std::string &name, const Vector2f &value) const { glProgramUniform2f(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1]); }
	void setUniform(const std::string &name, const Vector3f &value) const { glProgramUniform3f(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2]); }
	void setUniform(const std::string &name, const Vector4f &value) const { glProgramUniform4f(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2], value[3]); }

	void setUniform(const std::string &name, const GLint value) const { glProgramUniform1i(_program, glGetUniformLocation(_program, name.c_str()), value); }
	void setUniform(const std::string &name, const Vector2i &value) const { glProgramUniform2i(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1]); }
	void setUniform(const std::string &name, const Vector3i &value) const { glProgramUniform3i(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2]); }
	void setUniform(const std::string &name, const Vector4i &value) const { glProgramUniform4i(_program, glGetUniformLocation(_program, name.c_str()), value[0], value[1], value[2], value[3]); }

	void setUniform(const std::string &name, const GLuint value) const { glProgramUniform1ui(_program, glGetUniformLocation(_program, name.c_str()), value); }

	void setUniform(const std::string &name, const Matrix2f &value) const { glProgramUniformMatrix2fv(_program, glGetUniformLocation(_program, name.c_str()), 1, GL_FALSE, value.data()); }
	void setUniform(const std::string &name, const Matrix3f &value) const { glProgramUniformMatrix3fv(_program, glGetUniformLocation(_program, name.c_str()), 1, GL_FALSE, value.data()); }
	void setUniform(const std::string &name, const Matrix4f &value) const { glProgramUniformMatrix4fv(_program, glGetUniformLocation(_program, name.c_str()), 1, GL_FALSE, value.data()); }

protected:

	void checkCompileErrors(GLuint object, const std::string &type) const;
};

}
