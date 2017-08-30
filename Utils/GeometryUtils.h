#pragma once
#include <cmath>
#include <tuple>

namespace utils
{
	template <typename T>
	inline double_t dotProduct(T a, T b)
	{
		return a.x * b.x + a.y * b.y;
	}

	template <typename T>
	inline double_t dotProduct(T a, T b, T c)
	{
		return dotProduct(b - a, c - a);
	}

	template <typename T>
	inline double_t dotProduct(T a)
	{
		return dotProduct(std::get<0>(a), std::get<1>(a), std::get<2>(a));
	}

	template <typename T>
	inline double_t crossProduct(T a, T b)
	{
		return a.x * b.y - a.y * b.x;
	}

	template <typename T>
	inline double_t crossProduct(T a, T b, T c)
	{
		return crossProduct(b - a, c - a);
	}

	template <typename T>
	inline double_t crossProduct(T a)
	{
		return crossProduct(std::get<0>(a), std::get<1>(a), std::get<2>(a));
	}

	template <typename T>
	inline double_t length(T a)
	{
		return sqrt(a.x * a.x + a.y * a.y);
	}

	template <typename T>
	inline double_t distance(T a, T b)
	{
		return length(a - b);
	}

	template <typename T>
	inline T nomalize(T a)
	{
		return a / length(a);
	}
}
