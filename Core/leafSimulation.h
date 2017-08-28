#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <algorithm>
#include <cmath>
#include <vector>
#include <random>
#include <tuple>
#include <type_traits>

#include "IPointProcessor.h"
#include "../Utils/GeometryUtils.h"

namespace leaf
{
	template <typename DataType, class PointProcessor>
	class LeafSimulation
	{
		static_assert(std::is_base_of<IPointProccessor, PointProcessor>::value, "PointProcessor must be derived from IPointProcessor");

	private:
		// The shape of the leaf.
		// Describes a polygon, and the index 0 point should be the base.
		std::vector<DataType> shape;
		
		// Decompose the leaf polygon into triangles.
		std::vector<std::tuple<double_t, double_t, double_t>> triangles;

		PointProcessor *auxin, *veinNodes;
		
		// maintain the connection relationship of vein nodes.
		boost::adjacency_list<boost::listS, boost::vecS, boost::directedS> veinTopology;
		
		// maintain the mininum distance from every vein node to the base of the leaf.
		vector<double> minDistToBase;

	public:
		LeafSimulation()
		{
			auxin = new PointProcessor<DataType>();
			veinNodes = new PointProcessor<DataType>();
		}

		~LeafSimulation()
		{
			delete auxin;
			delete veinNodes;
		}

		// expand the vein for one step.
		// Parameters:
		// auxinGenerates:   numbers of new auxin nodes to be generated.
		// stepSize:         how far should new vein nodes goes.
		// auxinRadius:      the control radius of every auxin node.
		// auxinAttenuation: the attenuation coefficient of the auxin. We use a exponent function as 
		//                   attenuation function, i.e. auxin_strenth = e ^ (-distance * coefficient).
		// bs:               the minimum distance the new auxin nodes should keep from other auxin nodes.
		// bv:               the minimum distance the new auxin nodes should keep from other vein nodes.
		void veinGrow(uint32_t auxinGenerates, double_t stepSize, double_t auxinRadius, 
			double_t auxinAttenuation, double_t bs, double_t bv)
		{
			// Generate auxins.
			// First calculate the area of the polygon triangles 
			// so that random auxin nodes could be distributed evenly.
			static std::vector<double_t> areas;
			areas.clear();
			areas.reserve(triangles.size());
			for (auto i = triangles.begin(); i != triangles.end(); ++i)
			{
				areas.push_back(abs(utils::crossProduct(*i)) + *(areas.end() - 1));
			}
			for (auto i = areas.begin(); i != areas.end(); ++i)
				*i /= *(areas.end() - 1);
			
			// Generate auxins by picking a point randomly in the polygon and checking 
			// whether the point is legal.
			// Because we don't want a dead loop, the total number of attemps is limited.
			static std::uniform_real_distribution<double_t> unif(0.0, 1.0);
			static std::default_random_engine re;
			for (uint32_t generated = 0, i = 0; generated < auxinGenerates && i < auxinGenerates * 5; ++i)
			{
				// We first randomly pick a triangle.
				auto randTriangle = unif(re);
				auto picked = std::lower_bound(areas.begin(), areas.end(), randTriangle);
				// Then generate a point inside the triangle.
				double_t r1 = sqrt(unif(re)), r2 = unif(re);
				DataType p = std::get<0>(picked) * (1 - r1) + 
					         std::get<1>(picked) * (r1 * (1 - r2)) + 
					         std::get<2>(picked) * (r1 * r2);
				// Check the legality of the generated point.
				if (utils::distance(auxin->getNearest(p), p) > bs &&
					utils::distance(veinNodes->getNearest(p), p) > bv)
				{
					auxin->insert(p);
					generated++;
				}
			}

			// For every vein node, get all auxin nodes that controls it.
			static std::vector<std::vector<DataType>> controlledBy;
			controlledBy.resize(veinNodes->rawData().size());
			for (auto i = controlledBy.begin(); i != controlledBy.end(); ++i)
				i->clear();
			for (auto i = auxin->rawData().begin(); i != auxin->rawData().end(); ++i)
			{
				auto controls = veinNodes->getRelativeNeighborByIndex(i);
				for (auto j = controls.begin(); j != controls.end(); ++j)
					controlledBy[*j].push_back(*i);
			}

		}
	};
}
