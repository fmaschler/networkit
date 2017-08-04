/* GCE.h
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#ifndef GCE_H_
#define GCE_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "../auxiliary/SetIntersector.h"


namespace NetworKit {


/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GCE: public NetworKit::SelectiveCommunityDetector {

public:

	GCE(const Graph& G, std::string objective);


	std::map<node, std::set<node> >  run(const std::set<node>& seeds) override;

	/**
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	std::set<node> expandSeed(node s);

	NetworKit::Partition runPartition(const std::set<node>& seeds);
	 /**
	  * Not Implemented yet, returns empty Partition!
	 * @param seed Seed nodes for which communities are to be found.

	 * @return Partition that makes up the best communities found around nodes @a seed.
	 */

protected:

    std::string objective;    // name of objective function


};

} /* namespace NetworKit */
#endif
