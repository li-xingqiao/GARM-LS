#pragma once
#pragma warning (disable : 4996)

#include "Types.h"

#include <yaml-cpp/yaml.h>

namespace YAML {

template <typename Derived, int Rows>
struct convert<Eigen::Matrix<Derived, Rows, 1>>
{
	static Node encode(const Eigen::Matrix<Derived, Rows, 1> &rhs)
	{
		Node node;
		for (int i = 0; i < rhs.size(); i++) {
			node.push_back(rhs[i]);
		}
		node.SetStyle(EmitterStyle::Flow);
		return node;
	}

	static bool decode(const Node &node, Eigen::Matrix<Derived, Rows, 1> &rhs)
	{
		if (!node.IsSequence() || node.size() != rhs.size()) {
			return false;
		}
		for (int i = 0; i < rhs.size(); i++) {
			rhs[i] = node[i].as<Derived>();
		}
		return true;
	}
};

}
