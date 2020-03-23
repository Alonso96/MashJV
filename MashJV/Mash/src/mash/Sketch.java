/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mash;
import java.io.File;
import mash.Hash_U;
import opennlp.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Queue;
import java.util.TreeMap;

import javax.swing.text.StyledEditorKit;

import com.sun.java.swing.plaf.windows.resources.windows;
/**
 *
 * @author alfonso
 */
public class Sketch {
    final static String capnpHeader = "Cap'n Proto";
    final static Integer capnpHeaderLength = capnpHeader.length();
    final static String suffixSketch = ".msh";
    final static String suffixSketchWindowed = ".msw";
    final static String alphabetNucleotide = "ACGT";
    final static String alphabetProtein = "ACDEFGHIKLMNPQRSTVWY";
    private ArrayList<Reference> references = new ArrayList<Reference>();
    private Map<Long,ArrayList<Locus>> lociByHash = new HashMap<Long,ArrayList<Locus>>();
    Map<Long,ArrayList<PositionHash>> LociByHash_map = new HashMap<Long,ArrayList <PositionHash>>();
    ArrayList<ArrayList<PositionHash>> positionHashesByReference = new ArrayList<ArrayList<PositionHash>>();
    Parameters parameters;
    double kmerSpace;
    String file;
    private Map<String,Integer> referenceIndecesById = new HashMap<String,Integer>();
    
    
   // Hashtable<long, long> hash_t = new Hashtable<>();
    class Parameters{
        public Parameters(){
            parallelism =1;
            kmerSize =0;
            alphabetSize =0;
            preserveCase=false;
            use64 = false;
            seed = 0;
            error =0;
            warning =0;
            minHashesPerWindow=0;
            windowSize=0;
            windowed =false;
            concatenated=false;
            noncanonical =false;
            reads = false;
            memoryBound =0;
            minCov=1;
            targetCov=0;
            genomeSize= 0;
            //memseT(alphabet,0,256);
            
        }
        public Parameters(Parameters other){
            parallelism =other.parallelism;
            kmerSize =other.kmerSize;
            alphabetSize =other.alphabetSize;
            preserveCase=other.preserveCase;
            use64 = other.use64;
            seed = other.seed;
            error =other.error;
            warning =other.warning;
            minHashesPerWindow=other.minHashesPerWindow;
            windowSize=other.windowSize;
            windowed =other.windowed;
            concatenated=other.concatenated;
            noncanonical =other.noncanonical;
            reads = other.reads;
            memoryBound =other.memoryBound;
            minCov=other.minCov;
            targetCov=other.targetCov;
            genomeSize= other.genomeSize;
        
        }
        int parallelism;
        int kmerSize;
       boolean[] alphabet = new boolean[256]; //[256]; //controllare
       
        int alphabetSize;
        boolean preserveCase;
        boolean use64;
        int seed;
        double error;
        double warning;
        long minHashesPerWindow;
        long windowSize;
        boolean windowed;
        boolean concatenated;
        boolean noncanonical;
        boolean reads;
        long memoryBound;
        int minCov;
        double targetCov;
        long genomeSize;
    }// Parameters
    
    class PositionHash{
        Hashtable hash ;
        int position;
        public PositionHash(int positionNew, Hashtable hashNew){
            position=positionNew;
            hash = hashNew;
        }

    } //PositionHash
    class Locus{
        public Locus(int sequenceNew,int positionNew){
             sequence = sequenceNew;
             position = positionNew;
        }
        int sequence;
        int position;
    }
    class Reference{
        String name;
        String comment;
        long length;
       // HashList hashesSorted;
        ArrayList<Integer> counts = new ArrayList<Integer>();
    }
    class SketchInput{
        ArrayList<String> fileNames = new ArrayList<String>();
        char seq;
        long length;
        String name;
        String comment;
        Parameters parameters;
        public SketchInput(ArrayList<String> fileNamesNew,char seqNew, long lenghNew,String nameNew, String commentNew,Parameters parametersNew){
            fileNames= fileNamesNew;
            seq = seqNew;
            length = lenghNew;
            name = nameNew;
            comment = commentNew;
            parameters = parametersNew;
        } //costruttore
    } //SketchInput
    class SketchOutput{
        Reference references;
        ArrayList<PositionHash> positionHashesByReference;
    }
    public void getAlphabetAsString(String alphabet){
    	StringBuffer newString = new StringBuffer(alphabet); 

    	for (int i=0;i<256;i++) {
    		if(parameters.alphabet[i]) {
    			 newString.insert(1,i); //da verificare
    			 alphabet= newString.toString();
    		}
    	}
    }
    
    public Integer getAlphabetSize(){
        return parameters.alphabetSize;
    }
    
    public boolean getConcatenated(){
        return Parameters.concatenated;
    }
    
   public float getError(){
       return (float) parameters.error;
   }
   
   public Integer getHashCount(){
       return lociByHash.size();
   }
   
   public Integer getHashSeed(){
       return parameters.seed;
   }
   
   public ArrayList<Locus> getLociByHash(long hash){
       return lociByHash.get(hash);
   }
   
   public float getMinHashesPerWindow(){
       return parameters.minHashesPerWindow;
   }
   
   public int getMinKmerSize(long reference){
	   return (int) Math.ceil(Math.log(references.get((int)reference).length * (1 - parameters.warning) / parameters.warning) / Math.log(parameters.alphabetSize));
   }
	   public boolean getPreserveCase(){
        return parameters.preserveCase;
    }
	   
    public double getRandomKmerChance(long reference){
    	return 1. / (kmerSpace / references.get((int)reference).length + 1.);
    }
    
    public Reference getReference(long index){
            return references.get((int) index);
    }
    
    public long getReferenceCount()  {
        return references.size();
    }
    
    public  void getReferenceHistogram(long index, HashMap<Integer,Long> histogram){
    	Reference reference = references.get((int) index);
    	histogram.clear();
    	for(long i=0;i<reference.counts.size();i++) {
    		int count = reference.counts.get((int)i);
    		if(Collections.frequency((Collection<?>) histogram, count)==0) {
    			histogram.put(count, (long) 1);
    		}else
    		{
    			histogram.put(count, histogram.get(count)+1);
    		}
    	}
    	
    }//reference Histogram
    public long getReferenceIndex(String id){
    	if(Collections.frequency((Collection<?>) referenceIndecesById, id) ==1) {
    		return referenceIndecesById.get(id);
    	}else return -1;
    } 
    public int getKmerSize()  {
        return parameters.kmerSize;
    }
    public double getKmerSpace() {
        return kmerSpace;
    }
    public boolean getUse64(){
        return parameters.use64;
    }
    public long getWindowSize() {
        return parameters.windowSize;
    }
    public boolean getNoncanonical(){
        return parameters.noncanonical;
    }
    public boolean hasHashCounts(){
        return references.size() > 0 && references.get(0).counts.size() > 0;
    }
  /*  public boolean hasLociByHash(Hashtable hash) {
        return Collections.frequency((Collection<?>) lociByHash, hash);
    }
    */
    public int initFromFiles(ArrayList<String> files, Parameters  parametersNew, int verbosity , boolean enforceParameters, boolean contain){ //argDefaultint verbosity = 0, bool enforceParameters = false, bool contain = false
    	parameters = parametersNew;
    	for(int i=0; i<files.size();i++) {
    		boolean isSketch = hasSuffix(files.get(i),parameters.windowed ? suffixSketchWindowed:suffixSketch);
    		if(isSketch) {
    			// init header to check params
    			//
    			Sketch sketchTest = null;
    			sketchTest.initParametersFromCapnp(files.get(i));
    			
            	if ( i == 0 && ! enforceParameters )
            	{
            		initParametersFromCapnp(files.get(i));
            	}
            	
                String alphabet=null ;
                String alphabetTest = null;
                
                getAlphabetAsString(alphabet);
                sketchTest.getAlphabetAsString(alphabetTest);
                
                if ( alphabet != alphabetTest )
                {
                	System.err.println( "\n WARNING: The sketch file " + files.get(i) + " has different alphabet ("+ alphabetTest + ") than the current alphabet ("+ alphabet + "). This file will be skipped.");
                	continue;
                }
                if(sketchTest.getHashSeed()!= parameters.seed) {
                	System.err.println("\n WARNING: The sketch "+files.get(i)+ "has a seed size("+ sketchTest.getHashSeed() +") than the current alphabet ("+ alphabet +"). This file will be skipped");
                	continue;
                }
                if(sketchTest.getKmerSize()!= parameters.kmerSize) {
                	System.err.println("\nWARNING: The sketch " + files.get(i)+ " has a kmer size " +sketchTest.getKmerSize() +") that does not match the current kmer size (" + parameters.kmerSize +"). This file will be skipped.");
                	//iniziare da qui/
                	continue;
                }
                if(! contain && sketchTest.getMinHashesPerWindow() < parameters.minHashesPerWindow){
                	System.err.println("\n WARNING: The sketch file "+ files.get(i)+ " has a target sketch size (" +sketchTest.getMinHashesPerWindow() +") that is smaller than the current sketch size ("+parameters.minHashesPerWindow+"). This file will be skipped.");
                	continue;
                }
                if(sketchTest.getNoncanonical()!= parameters.noncanonical) {
                	System.err.println("\nWARNING The sketch file " + files.get(i) +"is" +(sketchTest.getNoncanonical() ? "noncanonical" : "canonical") + ", which is incompatible with the current setting. This file will be skipped.");
                	
                }
                if (sketchTest.getMinHashesPerWindow() > parameters.minHashesPerWindow) {
                	System.err.println("\n WARNING: The sketch file "+ files.get(i)+ " has a target sketch size (" +sketchTest.getMinHashesPerWindow() +") that is larger than the current sketch size ("+parameters.minHashesPerWindow+"). This file will be skipped.");
                }
                ArrayList<String> file = new ArrayList<String>();
                file.add(files.get(i));
                //threadPool.runWhenThreadAvailable(new SketchInput(file, 0, 0, "", "", parameters), loadCapnp);
                
                
    		}//ifIsSketch 
    		else {
    			File inStream =null;
    			if(files.get(i).equals("-")) {
    				if (verbosity>0) System.err.println("Sketching from stdin...");
    				inStream= new File(System.in.toString());
    				
    			}else {
    				if(verbosity>0) {
    					System.err.println("Sketching "+ files.get(i)+"...");
    				}
    				try{
    					inStream =new File(files.get(i).toString());
    				}catch(Exception e) {
    					System.err.println("ERROR: could not open "+ files.get(i)+ "for reading");
    					e.printStackTrace();
    					System.exit(1);
    					
    				}
    			}
    			if(parameters.concatenated) {
    				if(!files.get(i).equals("-")) {
    					inStream= null;
    				}
    				ArrayList<String> file = new ArrayList<String>();
    				file.add(files.get(i));
    				sketchFile(new SketchInput(file,'0',0,"","",parameters));
    				}else
    				{
    				//	ThreadPool<SketchInput, SketchOutput> threadPool;
						if(!sketchFileBySequence(inStream, null)) {
    						System.err.println("\nERROR: reading" + files.get(i));
    						System.exit(1);
    					}
						inStream=null;
						
    				}
    			}
    		} 
    		createIndex();
    		return 0;
    
    }
    public void initFromReads(ArrayList<String> files,Parameters parametersNew){
    	parameters = parametersNew;
    	useThreadOutput(sketchFile(new SketchInput(files, (char) 0, 0, "", "", parameters)));
    	createIndex();
    	
    }
    long initParametersFromCapnp(String file){ // con capnproto da vedere meglio
    }
    public void setReferenceName(int i, String name) {
            references.get(i).name = name;
    }
    public void setReferenceComment(int i, String comment) {
        references.get(i).comment = comment;
    }
    public boolean sketchFileBySequence(File file, ThreadPool<SketchInput, SketchOutput>  threadPool) { // da capire meglio cosa fa nel codice C
    	
    	
    }
    public void useThreadOutput(SketchOutput output){}
    public void warnKmerSize(long  lengthMax, String lengthMaxName, double randomChance, int kMin, int warningCount){}
    public boolean writeToFile(){} 
    int writeToCapnp(String file){} // capnproto
    private void createIndex(){
    	
    	for (int i= 0 ; i< references.size(); i++) {
    		referenceIndecesById.put(references.get(i).name,i);
    	}
    	for ( int i = 0; i < positionHashesByReference.size(); i++ )
        {
            for ( int j = 0; j < positionHashesByReference.get(i).size(); j++ )
            {
                PositionHash positionHash = positionHashesByReference.get(i).get(j);
                
                lociByHash.get(positionHash.hash).add(new Locus(i, positionHash.position));
            }
        }
        
        kmerSpace = Math.pow(parameters.alphabetSize, parameters.kmerSize);
    }
    private void addMinHashes(MinHashHeap minHashHeap, char [] seq, long length, Parameters  parameters) {
    	int kmerSize = parameters.kmerSize;
    	long mins = parameters.minHashesPerWindow;
    	boolean noncanonical = parameters.noncanonical;
    	for( long i =0; i< length; i++) {
    		if(! parameters.preserveCase && seq[(int) i]  > 96 && seq[(int)i]< 123)
    			
    			seq[(int)i] -= 32;
    	} // fine for
    	char seqRev [] = null;
    	if(!noncanonical) {
    		seqRev = new char [(int) length];
    		reverseComplement(seq, seqRev, length);
    	}
    	long j =0;
    	for(long i=0; i<length -kmerSize +1 ; i++) {
    		//repeatedly skip kmers with bad characters
    		boolean bad = false;
    		for(; j < i + kmerSize && i + kmerSize <= length; j++ ) {
    			if(!parameters.alphabet[seq[(int)j]]) {
    				i = j++;
    				bad= true;
    				break;
    			}
    		}
    		if(bad) continue;
    		if(i+kmerSize > length) break;
    		
    		final String kmer_fwd = seq.toString().substring((int) (seq.length-i));
    	
    		final String kmer_rev = seqRev.toString().substring((int) (seqRev.length + length - i -kmerSize));
    		final String kmer = (noncanonical || memcmp(kmer_fwd.substring(0,kmerSize).getBytes(), kmer_rev.substring(0,kmerSize).getBytes()) <= 0) ? kmer_fwd : kmer_rev;
    		boolean filter = false;
    	} // end for
    	if (!noncanonical) seqRev = null;
    	
    }
    private void getMinHashPositions(ArrayList<PositionHash>  loci, char [] seq, int length, Parameters  parameters, int verbosity ) { //verbosity = 0 default
    
    	int kmerSize = parameters.kmerSize;
    	    int mins = (int) parameters.minHashesPerWindow;
    	    int windowSize = (int) parameters.windowSize;
    	    
    	    int nextValidKmer = 0;
    	    
    	    if ( windowSize > length - kmerSize + 1 )
    	    {
    	        windowSize = length - kmerSize + 1;
    	    }
    	    
    	    if ( verbosity > 1 ) System.out.println(seq);
    	    // Associate positions with flags so they can be marked as min-hashes
    	    // at any point while the window is moved across them
    	    //
    	   class CandidateLocus
    	    {
    	        public CandidateLocus(int positionNew) {
    	            
    	            position=positionNew;
    	            isMinmer=false;
    	            }
    	        
    	        int position;
    	        boolean isMinmer;
    	    }
    	   NavigableMap<Long, Deque<CandidateLocus>> candidatesByHash = new TreeMap<Long,Deque<CandidateLocus>>(); //Uso treemap per poter accedere all'ulimo elemento inserito
    	   Deque<CandidateLocus> dq;
    	  // Iterator it = dq.iterator();	
    	   Queue<Map<Long,Iterator<Deque<CandidateLocus>>>> windowQueue;
    	   //Entry<String, Integer> lastEntry = map.lastEntry();
    	   HashMap<Long,Iterator<Deque<CandidateLocus>>> maxMinmer = (HashMap<Long, Iterator<Deque<CandidateLocus>>>) candidatesByHash.lastEntry();
    	   HashMap<Long,Iterator<Deque<CandidateLocus>>> newCandidates;
    	   int unique = 0;
    	   for ( int i = 0; i < length - kmerSize + 1; i++ )
    	    {
    	        // Increment the next valid kmer if needed. Invalid kmers must still be
    	        // processed to keep the queue filled, but will be associated with a
    	        // dummy iterator. (Currently disabled to allow all kmers; see below)
    	        //
    	        if ( i >= nextValidKmer )
    	        {
    	            for ( int j = i; j < i + kmerSize; j++ )
    	            {
    	                char c = seq[j];
    	                
    	                if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' )
    	                {
    	                    // Uncomment to skip invalid kmers
    	                    //
    	                    //nextValidKmer = j + 1;
    	                    
    	                    break;
    	                }
    	            }
    	        }
    	        if ( i < nextValidKmer && verbosity > 1 )
    	        {
    	            System.out.println( "  [");
    	        
    	            for ( int j = i; j < i + kmerSize; j++ )
    	            {
    	            	System.out.println(seq[j]);
    	            }
    	            
    	            System.out.println( "  ]");
    	        }
    	        if(i >= nextValidKmer) {
    	        	Hash_U h = new Hash_U();
    	        	long hash =  h.getHash(seq.toString().substring((int) (seq.length-i)), kmerSize, parameters.seed, parameters.use64);
    	        	 if ( verbosity > 1 )
    	             {
    	                 System.out.println("  ");
    	             
    	                 for ( int j = i; j < i + kmerSize; j++ )
    	                 {
    	                	 System.out.println(seq[j]);
    	                 }
    	             
    	                 System.out.println("  "+ i + '\t' +hash);
    	             }
    	             
    	        }
    	    }//fine for
    	   
    } //verbosity =0 default
    
    private void getMinHashPositions(ArrayList<PositionHash>  loci, char [] seq, int length, Parameters  parameters) { //overload per parametro opzionale
    	getMinHashPositions(loci, seq, length, parameters,0);
    }
    private boolean hasSuffix(String whole, String suffix) {} 
    private SketchOutput loadCapnp(SketchInput  input) {} // capnproto da vedere meglio
    private void reverseComplement( char[] seq, char[]  seqRev, long length) {} 
    private void setAlphabetFromString(Parameters parameters, char  characters) {}
    void setMinHashesForReference(Reference  reference, MinHashHeap  hashes) {}
    private SketchOutput sketchFile (SketchInput  input) {
    	
    }
    private SketchOutput sketchSequence(SketchInput  input) {
    	
    }
    
    public static int memcmp(final byte[] a, final byte[] b) {
        final int length = Math.min(a.length, b.length);
        if (a == b) {  // Do this after accessing a.length and b.length
          return 0;    // in order to NPE if either a or b is null.
        }
        for (int i = 0; i < length; i++) {
          if (a[i] != b[i]) {
            return (a[i] & 0xFF) - (b[i] & 0xFF);  // "promote" to unsigned.
          }
        }
        return a.length - b.length;
      }
      
    
}
