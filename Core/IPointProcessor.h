#pragma once
#include <vector>

namespace leaf
{
	template <typename DataType>
	class IPointProccessor
	{
		// Query functions.
		virtual DataType& getNearest(const DataType &a) = 0;
		virtual std::vector<uint32_t>& getRelativeNeighborByIndex(const DataType &a, double_t controlRadius) = 0;

		// Edit functions.
		virtual void insert(const DataType &a) = 0;
		virtual void rebuild(const std::vector<DataType> &a) = 0;
	};
}
