#pragma once
#include <vector>

namespace leaf
{
	template <typename DataType, typename RawType>
	class IPointProccessor
	{
		// Query functions.
		virtual DataType& getNearest(const DataType &a) = 0;
		virtual std::vector<uint32_t>& getRelativeNeighborByIndex(const DataType &a, double_t controlRadius) = 0;

		// Edit functions.
		virtual void insert(const DataType &a) = 0;

		// We could using rawData() to get all points, modify it, and rebuild the data structure.
		virtual RawType& rawData() = 0;
		virtual void rebuild() = 0;
	};
}
