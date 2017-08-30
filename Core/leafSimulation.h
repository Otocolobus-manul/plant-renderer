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
#include "../Utils/defines.h"

namespace leaf
{
	template <typename DataType, class PointProcessor>
	class LeafSimulation
	{
		static_assert(std::is_base_of<IPointProccessor, PointProcessor>::value, "PointProcessor must be derived from IPointProcessor");

	private:
		// Data structure for saving information of a node on the boundry of the leaf.
		struct edgeNode
		{
			// node position
			DataType pos;
			// density of auxin on the point
			double_t auxinDensity;
			// Auxin node index in WrappedProcessor *auxin.
			// When leaf grows, we should also modify data used in veinGrow() function.
			int32_t auxinNode;
		};
		// information of the edge of the leaf.
		// Describes a polygon, and the first should be the base.
		std::vector<edgeNode> edge;
		
		// So I decided to put this into veinGrow() function.
		// The partial sum of auxin density on edge, being used by vein growing function 
		// to generate random auxin nodes on the edge.
		// TODO: this could be maintained by BST if it constitutes a performance bottleneck.
		// std::vector<double_t> edgeAuxinPartialSum;

		// Decompose the leaf polygon into triangles.
		std::vector<std::tuple<double_t, double_t, double_t>> triangles;

		// Wrapped point processor class.
		// So we can get both unarranged data and arranged data.
		class WrappedProcessor : PointProcessor
		{
		public:
			std::vector<DataType> rawData;
			virtual void insert(const DataType &a)
			{
				PointProcessor::insert(a);
				rawData.push_back(a);
			}
			void rebuild()
			{
				PointProcessor::rebuild(rawData);
			}
		};
		WrappedProcessor *auxin, *veinNodes;
		
		// maintain the connection relationship of vein nodes.
		boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS> veinTopology;
		
		// For generating closed patterns, serveral veins may grow to the same auxin source. 
		// We use this vector to save the correspondence from auxin node index to vein node index.
		std::vector<int32_t> auxinNode2VeinIndex;

		// maintain the mininum distance from every vein node to the base of the leaf.
		vector<double_t> minDistToBase;

	public:
		LeafSimulation()
		{
			auxin = new WrappedProcessor<DataType>();
			veinNodes = new WrappedProcessor<DataType>();
		}

		~LeafSimulation()
		{
			delete auxin;
			delete veinNodes;
		}

		void readyForSimulation()
		{
			auxinNode2VeinIndex.clear();
		}

		// expand the vein for one step.
		// Parameters:
		// edgeAuxins:       numbers of new auxin nodes generates on the edge of the leaf.
		// innerAuxins:      numbers of new auxin nodes generates on the internal of the leaf.
		// auxinGenerates:   numbers of new auxin nodes to be generated.
		// stepSize:         how far should new vein nodes goes.
		// auxinRadius:      the control radius of every auxin node.
		// auxinAttenuation: the attenuation coefficient of the auxin. We use a exponent function as 
		//                   attenuation function, i.e. auxin_strenth = e ^ (-distance * coefficient).
		// bs:               the minimum distance the new auxin nodes should keep from other auxin nodes.
		// bv:               the minimum distance the new auxin nodes should keep from other vein nodes.
		void veinGrow(
			uint32_t edgeAuxins,
			uint32_t innerAuxins, 
			double_t stepSize, 
			double_t auxinRadius, 
			double_t auxinAttenuation, 
			double_t bs, 
			double_t bv)
		{
			// Generate edge auxins.
			// Please refer to the following internal auxin generation process.
			static std::vector<double_t> areas;
			// TODO: The following code should be deleted in the next submittion.
			/*areas.clear();
			areas.reserve(edgeAuxin.size());
			int32_t index = 0;
			FOR_EACH(i, edgeAuxin)
			{
				auto area = (*i + *(i + 1)) * utils::distance(edge[index], edge[index + 1]);
				areas.push_back(area + (i == edgeAuxin.begin() ? 0 : *(areas.end() - 1)));
				index++;
			}
			FOR_EACH(i, edgeAuxin)
				*i /= *(areas.end() - 1);
			static std::uniform_real_distribution<double_t> unif(0.0, 1.0);
			static std::default_random_engine re;
			for (uint32_t generated = 0, i = 0; generated < edgeAuxins && i < edgeAuxins * 5; ++i)
			{
				auto randSegment = unif(re);
				auto picked = std::lower_bound(areas.begin(), areas.end(), randSegment);
				double_t r = unif(re), low_ = edgeAuxin[picked], high_ = edgeAuxin[picked + 1];
				double_t low = 2 * low_ / (low_ + high_), delta = 2 * (high_ - low_) / (low_ + high_);
				r = low * r + 0.5 * delta * r * r;
				DataType p = shape[picked + 1] * r + shape[picked] * (1 - r);
				if (utils::distance(auxin->getNearest(p), p) > bs &&
					utils::distance(veinNodes->getNearest(p), p) > bv)
				{
					auxin->insert(p);
					generated++;
				}
			}*/

			// Generate internal auxins.
			// First calculate the area of the polygon triangles 
			// so that random auxin nodes could be distributed evenly.
			areas.clear();
			areas.reserve(triangles.size());
			FOR_EACH(i, triangles)
			{
				areas.push_back(abs(utils::crossProduct(*i)) + (i == triangles.begin() ? 0 : *(areas.end() - 1)));
			}
			FOR_EACH(i, areas)
				*i /= *(areas.end() - 1);
			
			// Generate auxins by picking a point randomly in the polygon and checking 
			// whether the point is legal.
			// Because we don't want a dead loop, the total number of attemps is limited.
			for (uint32_t generated = 0, i = 0; generated < innerAuxins && i < innerAuxins * 5; ++i)
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
			auxinNode2VeinIndex.resize(auxin->rawData.size(), -1);

			// For every vein node, get all auxin nodes that controls it.
			static std::vector<std::vector<std::tuple<uint32_t, DataType>>> controlledBy;
			controlledBy.resize(veinNodes->rawData.size());
			for (auto i = controlledBy.begin(); i != controlledBy.end(); ++i)
				i->clear();
			uint32_t index = 0;
			FOR_EACH(i, auxin->rawData)
			{
				auto controls = veinNodes->getRelativeNeighborByIndex(*i);
				if (auxinNode2VeinIndex[index] == -1)
					for (auto j = controls.begin(); j != controls.end(); ++j)
						controlledBy[*j].push_back(std::make_tuple(index, *i));
				else
				{
					// We should consider the situation, that the segments of veins rather than
					// nodes of veins intersect with the forbidden area.
					// So first we find all the vein nodes that connect to the one that coincides
					// with current auxin node.
					static std::vector<DataType> neighbor;
					auto iter = boost::in_edges(auxinNode2VeinIndex[index], veinTopology);
					for (; iter.first != iter.second; ++iter.first)
						neighbor.push_back(veinNodes->rawData[boost::get(boost::vertex_index, veinTopology, *iter)]);
					for (auto j = controls.begin(); j != controls.end(); ++j)
					{
						//                         . <- auxin node
						//                        /|
						//                       / |
						//                      /  |
						//                     /   |
						// added vein nodes-> .    . <- node to be checked
						// Then, we check the angle. If it's less than pi, then the auxin node 
						// can't control the checking node.
						bool flag = true;
						for (auto k = neighbor.begin(); k != neighbor.end(); ++k)
							if (utils::dotProduct(*i, *j, *k) > 0)
							{
								flag = false;
								break;
							}
						if (flag)
							controlledBy[*j].push_back(std::make_tuple(indes, *i));
					}
				}	
			}
			
			// Calculate the growing direction and grow.
			auto veinSize = veinNodes->rawData.size();
			for (int i = 0; i < veinSize; ++i)
			{
				DataType &curVein = veinNodes->rawData[i];
				DataType direction;
				double_t nearest = 1e+8;
				uint32_t nearestIndex = 0;
				DataType nearestP;
				FOR_EACH(controlledBy, j)
				{
					auto delta = std::get<1>(*j) - curVein;
					double_t distance = utils::length(delta);
					delta /= distance;
					direction += delta * exp(-distance * auxinAttenuation);
					if (distance < nearest)
					{
						nearest = distance;
						nearestIndex = std::get<0>(*j);
						nearestP = std::get<1>(*j);
					}
				}
				// if an auxin node is reachable in one step, then just reach it. 
				if (nearest < stepSize)
				{
					// Several veins may grow to the same auxin node.
					int32_t &veinIndex = auxinNode2VeinIndex[nearestIndex];
					if (veinIndex == -1)
					{
						veinIndex = veinNodes->rawData.size();
						veinNodes->insert(nearestP);
					}
					boost::add_edge(i, veinIndex, veinTopology);
				}
				direction = utils::normalize(direction) * stepSize;
				boost::add_edge(i, veinNodes->rawData.size(), veinTopology);
				veinNodes->insert(curVein + direction);
			}
		}
	};
}
