/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include "../community/Conductance.h"
#include "../auxiliary/Parallel.h"
#include <cmath>
#include <vector>

namespace NetworKit {

PageRankNibble::PageRankNibble(Graph& g, double alpha, double epsilon): SelectiveCommunityDetector(g), alpha(alpha), epsilon(epsilon) {
	graphVolume = 0.0;
	G.forNodes([&](node u) {
		graphVolume += G.volume(u);
	});
	TRACE("Graph Volume: ", graphVolume);
}

PageRankNibble::~PageRankNibble() {

}


std::pair<std::set<node>, double> PageRankNibble::bestSweepSet(std::vector<std::pair<node, double>>& pr) {
	TRACE("Support size: ", pr.size());


	// order vertices
	TRACE("Before sorting");
	auto comp([&](const std::pair<node, double>& a, const std::pair<node, double>& b) {
		return (a.second / G.weightedDegree(a.first)) > (b.second / G.weightedDegree(b.first));
	});
	Aux::Parallel::sort(pr.begin(), pr.end(), comp);
	TRACE("After sorting");

	for (std::vector<std::pair<node, double>>::iterator it = pr.begin(); it != pr.end(); it++) {
		TRACE("(", it->first, ", ", it->second, ")");
	}

	// find best sweep set w.r.t. conductance
	double bestCond = std::numeric_limits<double>::max();
	double cut = 0.0;
	double volume = 0.0;
	index bestSweepSetIndex = 0;
	std::unordered_map<node, bool> withinSweepSet;
	std::vector<node> currentSweepSet;

	for (std::vector<std::pair<node, double>>::iterator it = pr.begin(); it != pr.end(); it++) {
		// update sweep set
		node v = it->first;
		G.forNeighborsOf(v, [&](node neigh) {
			if (withinSweepSet.find(neigh) == withinSweepSet.end()) {
				cut += G.weight(v, neigh);
			} else {
				cut -= G.weight(v, neigh);
			}
		});
		volume += G.volume(v);
		currentSweepSet.push_back(v);
		withinSweepSet[v] = true;

		// compute conductance
		double cond = cut / fmin(volume, graphVolume - volume);

		std::stringstream debug;

		debug << "Current vertex: " << v << "; Current sweep set conductance: " << cond << std::endl;
		debug << "Current cut weight: " << cut << "; Current volume: " << volume << std::endl;
		debug << "Total graph volume: " << graphVolume << std::endl;

		TRACE(debug.str());

		if (cond < bestCond) {
			bestCond = cond;
			bestSweepSetIndex = currentSweepSet.size();
		}
	}

	std::set<node> bestSweepSet;

	for (index j = 0; j < bestSweepSetIndex; j++) {
		bestSweepSet.insert(currentSweepSet[j]);
	}
	return std::pair<std::set<node>, double>(bestSweepSet, bestCond);
}


std::set<node> PageRankNibble::expandSeed(node seed) {
	DEBUG("APR(G, ", alpha, ", ", epsilon, ")");
	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<std::pair<node, double>> pr = apr.run(seed);
	std::pair<std::set<node>, double> set_cond = bestSweepSet(pr);
	seed_cond.push_back(std::pair<node, double>(seed, set_cond.second));
	return set_cond.first;
}

std::map<node, std::set<node> > PageRankNibble::run(std::set<unsigned int>& seeds) {
	std::map<node, std::set<node> > result;
	seed_cond = std::vector<std::pair<node, double>>();
	for (auto seed : seeds) {
		auto community = expandSeed(seed);
		result[seed] = community;
	}
	return result;
}

Partition PageRankNibble::runPartition(std::set<unsigned int>& seeds) {
	auto result = run(seeds);
	Partition partition(G.upperNodeIdBound());
	partition.allToOnePartition();
	// sort conductance ascending
	auto comp([&](const std::pair<node, double>& a, const std::pair<node, double>& b) {
		return (a.second) < (b.second);
	});
	Aux::Parallel::sort(seed_cond.begin(), seed_cond.end(), comp);
	for (std::vector<std::pair<node, double>>::iterator it = seed_cond.begin(); it != seed_cond.end(); it++) {
		node seed = it->first;
		partition.toSingleton(seed);
		index id = partition[seed];
		auto cluster = result[seed];
		for (auto entry: cluster) {
			// only move unassigned nodes from first partition
			// leave partition with lower conductance
			if (partition[entry] == 0) {
				partition.moveToSubset(id, entry);
			}
		}
	}
	return partition;
}

} /* namespace NetworKit */
