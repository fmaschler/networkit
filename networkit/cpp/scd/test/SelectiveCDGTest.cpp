#include "SelectiveCDGTest.h"

#include "../PageRankNibble.h"
#include "../../community/Modularity.h"
#include "../../community/Conductance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(SCDGTest2, testPageRankNibble) {
	METISGraphReader reader;
	Graph G = reader.read("input/hep-th.graph");
	// parameters
	node seed = 50;
	std::set<unsigned int> seeds = {(unsigned int) seed};
	double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);

	PageRankNibble prn(G, alpha, epsilon);

	// run PageRank-Nibble and partition the graph accordingly
	DEBUG("Call PageRank-Nibble(", seed, ")");
	auto partition = prn.runPartition(seeds);

	EXPECT_GT(partition.numberOfSubsets(), 1u);
	int cluster_size = partition.subsetSizeMap()[partition[seed]];
	EXPECT_GT(cluster_size, 0u);

	// evaluate result
	Conductance conductance;
	double targetCond = 0.4;
	double cond = conductance.getQuality(partition, G);
	EXPECT_LT(cond, targetCond);
	INFO("Conductance of PR-Nibble: ", cond, "; seed partition size: ", cluster_size, "; number of partitions: ", partition.numberOfSubsets());
}

TEST_F(SCDGTest2, testWeightedPageRankNibble) {
	METISGraphReader reader;
	Graph wG = reader.read("input/lesmis.graph");
	// parameters
	node seed = 50;
	std::set<unsigned int> seeds = {(unsigned int) seed};
	double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);

	Graph G = wG.toUnweighted();
	PageRankNibble prn(G, alpha, epsilon);
	PageRankNibble wPrn(wG, alpha, epsilon);

    double cond = 0;
    double wCond = 0;
	// run PageRank-Nibble and partition the graph accordingly
	DEBUG("Call PageRank-Nibble(", seed, ")");
    auto partition = prn.runPartition(seeds);
    auto wPartition = wPrn.runPartition(seeds);

    EXPECT_GT(partition.numberOfSubsets(), 1u);
    EXPECT_GT(wPartition.numberOfSubsets(), 1u);
    int size = partition.subsetSizeMap()[partition[seed]];
    int wSize = wPartition.subsetSizeMap()[wPartition[seed]];
    EXPECT_GT(size, 0u);
    EXPECT_GT(wSize, 0u);

    // evaluate result
    Conductance conductance;
    cond = conductance.getQuality(partition, G);
    wCond = conductance.getQuality(wPartition, wG);
	EXPECT_LT(wCond, cond);
	INFO("Conductance of weighted PR-Nibble: ", wCond, "; unweighted: ", cond);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
