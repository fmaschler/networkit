/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include "../auxiliary/Parallel.h"
#include <vector>
#include <algorithm>
#include <unordered_set>

namespace NetworKit {

PageRankNibble::PageRankNibble(const Graph& g, double alpha, double epsilon): SelectiveCommunityDetector(g), alpha(alpha), epsilon(epsilon) {

}

PageRankNibble::~PageRankNibble() {

}


std::pair<std::set<node>, double> PageRankNibble::bestSweepSet(std::vector<std::pair<node, double>>& pr) {
	TRACE("Finding best sweep set. Support size: ",  pr.size());
	// order vertices
	TRACE("Before sorting");
	for (size_t i = 0; i < pr.size(); i++) {
		pr[i].second = pr[i].second / G.volume(pr[i].first);
	}
	auto comp([&](const std::pair<node, double>& a, const std::pair<node, double>& b) {
		return a.second > b.second;
	});
	Aux::Parallel::sort(pr.begin(), pr.end(), comp);
	TRACE("After sorting");

	#ifndef NDEBUG
	for (auto it = pr.begin(); it != pr.end(); it++) {
		TRACE("(", it->first, ", ", it->second, ")");
	}
	#endif

	// find best sweep set w.r.t. conductance
	double bestCond = std::numeric_limits<double>::max();
	double cut = 0.0;
	double volume = 0.0;
	index bestSweepSetIndex = 0;
	std::unordered_set<node> withinSweepSet;
	std::vector<node> currentSweepSet;

	// generate total volume.
	double totalVolume = G.totalEdgeWeight() * 2;

	for (auto it = pr.begin(); it != pr.end(); it++) {
		// update sweep set
		node v = it->first;
		double wDegree = 0.0;
		G.forNeighborsOf(v, [&](node, node neigh, edgeweight w) {
			wDegree += w;
			if (withinSweepSet.find(neigh) == withinSweepSet.end()) {
				cut += w;
			} else {
				cut -= w;
			}
		});
		volume += wDegree;
		currentSweepSet.push_back(v);
		withinSweepSet.insert(v);

		// compute conductance
		double cond = cut / std::min(volume, totalVolume - volume);

		if ((cond < bestCond) && (currentSweepSet.size() < G.numberOfNodes())) {
			bestCond = cond;
			bestSweepSetIndex = currentSweepSet.size();
		}
	}

	DEBUG("Best conductance: ", bestCond, "\n");

	std::set<node> bestSweepSet(currentSweepSet.begin(), currentSweepSet.begin() + bestSweepSetIndex);
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

std::map<node, std::set<node> > PageRankNibble::run(const std::set<node>& seeds) {
	std::map<node, std::set<node> > result;
	seed_cond = std::vector<std::pair<node, double>>();
	for (auto seed : seeds) {
		auto community = expandSeed(seed);
		result[seed] = community;
	}
	return result;
}

Partition PageRankNibble::runPartition(const std::set<node>& seeds) {
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
		index id_old = partition[seed];
		partition.toSingleton(seed);
		index id = partition[seed];
		std::set<node> cluster = result[seed];
		for (auto entry: cluster) {
			// only move unassigned nodes from first partition
			// leave partition with better conductance
			if (partition[entry] == 0) {
				partition.moveToSubset(id, entry);
			}
		}
		// revert partitions that only consist of the seed
		if (partition.getMembers(id).size() == 1) {
			INFO("Revert singleton parition at ", seed);
			partition.moveToSubset(id_old, seed);
		}
	}
	return partition;
}

} /* namespace NetworKit */
