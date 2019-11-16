package mash;

import java.util.HashMap;
import java.util.Map;

import opennlp.tools.util.HashList;

public class HashSet extends HashList{
	
	private static final long serialVersionUID = 1L;
	private boolean use64;
	HashMap<Integer,Integer>hashes32;
	HashMap<Long,Integer>hashes64;
	
	public HashSet(boolean use64New) {
		use64 = use64New;
	
	}
	
	    
	
	public void toHashList(HashList hashList) {
		if(use64) {
			for(Map.Entry<Long,Integer> i : hashes64.entrySet()) {
				hashList.put(i.getValue());
			}
		}
	}
	
	
}
