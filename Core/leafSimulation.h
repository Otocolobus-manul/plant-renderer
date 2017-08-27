#pragma once
#include "IPointProcessor.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <type_traits>

namespace leaf
{
	template <typename DataType, class PointProcessor>
	class LeafSimulation
	{
		static_assert(std::is_base_of<IPointProccessor, PointProcessor>::value, "PointProcessor must be derived from IPointProcessor");

	private:
		PointProcessor<DataType> *auxin, *veinNodes;
		
		// maintain the connection relationship of vein nodes.
		boost::adjacency_list<boost::listS, boost::vecS, boost::directedS> veinTopology;
		
		// maintain the mininum distance from every vein node to the base of the leaf.
		vector<double> minDistToBase;

	public:
		// expand the vein for one step.
		// Parameters:
		// auxinGenerates: numbers of new auxin nodes to be generated.
		// stepSize:       how far should new vein nodes goes.
		// auxinRadius:    the control radius of every auxin node.
		void veinGrow(int auxinGenerates, float stepSize, float auxinRadius)
		{

		}
	};
}
