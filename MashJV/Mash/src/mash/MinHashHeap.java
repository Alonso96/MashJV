package mash;

import java.util.ArrayList;
import java.util.AbstractMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

import com.google.common.hash.BloomFilter;
import opennlp.tools.util.HashList;

public class MinHashHeap {
	private boolean use64;
	
	private HashSet<Boolean> hashes;
	PriorityQueue<Boolean> hashesQueue;
	HashSet<Boolean> hashesPending;
	PriorityQueue<Boolean> hashesQueuePending;
	
	long cardinalityMaximum;
	long multiplicityMinimum;
	
	long multiplicitySum;
	
    BloomFilter<Integer> bloomFilter;
    
    long kmersTotal;
    long kmersUsed;
	public MinHashHeap(boolean use64New, long cardinalityMaximumNew, long multiplicityMinimumNew  /* =1 default*/, long memoryBoundBytes /*= 0default*/) {
		use64=use64New;
		hashes.add(use64New);
		hashesQueue.add(use64New);
		hashesPending.add(use64New);
		hashesQueuePending.add(use64New);
		cardinalityMaximum = cardinalityMaximumNew;
		multiplicityMinimum =multiplicityMinimumNew;
		multiplicitySum=0;
		if(memoryBoundBytes ==0) bloomFilter=null;
		else {
			/*
			 * bloom_parameters bloomParams;
		
		bloomParams.projected_element_count = 1000000000;//(uint64_t)parameters.genomeSize * 10l; // TODO: error rate based on platform and coverage
		bloomParams.false_positive_probability = 0;//parameters.bloomError;
		bloomParams.maximum_size = memoryBoundBytes * 8l;
		bloomParams.compute_optimal_parameters();
		*/
			kmersTotal=0;
			kmersUsed=0;
			
		}
	}
	public MinHashHeap(boolean use64New, long cardinalityMaximumNew) {
		long memoryBoundBytes;
		use64=use64New;
		hashes.add(use64New);
		hashesQueue.add(use64New);
		hashesPending.add(use64New);
		hashesQueuePending.add(use64New);
		cardinalityMaximum = cardinalityMaximumNew;
		multiplicityMinimum =1;
		multiplicitySum=0;
		memoryBoundBytes= 0;
		if(memoryBoundBytes ==0) bloomFilter=null;
		else {
			/*
			 * bloom_parameters bloomParams;
		
		bloomParams.projected_element_count = 1000000000;//(uint64_t)parameters.genomeSize * 10l; // TODO: error rate based on platform and coverage
		bloomParams.false_positive_probability = 0;//parameters.bloomError;
		bloomParams.maximum_size = memoryBoundBytes * 8l;
		bloomParams.compute_optimal_parameters();
		*/
			kmersTotal=0;
			kmersUsed=0;
			
		}
	}
	
	public void computeStats() {
		ArrayList<Integer> counts = new ArrayList<Integer>();
		hashes.contains(o)
	}
	public void clear() {}
	public double estimateMultiplicity{}
	public double etimateSetSize() {}
	public void toCounts(ArrayList<Long> counts) {}
	public void toHashList(Hashlist hashList)
}
