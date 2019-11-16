/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mash;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.math.*;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

import mash.Sketch.Parameters;
/**
 *
 * @author alfonso
 */
public class CommandContain extends Command {
    public CommandContain(){
        super();
        
        name = "within";
        summary = "Estimate the containment of query sequences within references.";
        description = "Estimate the containment of each query file (or sequence with -i) in the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or mash sketch files (.msh) with matching k-mer sizes. Query files can also be files of file names (see -l). The score is the number of intersecting min-hashes divided by the query set size. The output format is [score, error-bound, reference-ID, query-ID].";
        argumentString = "<reference> <query> [<query>] ...";
    
        addOption("list", new Option(Type.Boolean, "l", "Input", "List input. Each query file contains a list of sequence files, one per line. The reference file is not affected.", "",0,0));
        addOption("errorThreshold",new Option(Type.Number, "e", "Output", "Error bound threshold for reporting scores values. Error bounds can generally be increased by increasing the sketch size of the reference.", "0.05",0,0));
        useOption("help");
        useSketchOptions();
        
    }
    class ContainInput{
        Sketch sketchRef;
        Sketch sketchQuery;
        long indexRef;
        long indexQuery;
        long pairCount;
        String nameRef;
        Parameters parameters;
        public ContainInput(Sketch sketchRefNew, Sketch sketchQueryNew,long indexRefNew,long indexQueryNew, long pairCountNew, Parameters parametersNew){
            sketchRef= sketchRefNew;
            sketchQuery=sketchQueryNew;
            indexRef= indexRefNew;
            indexQuery= indexQueryNew;
            pairCount= pairCountNew;
            parameters= parametersNew;
        }
    }
    
    class ContainOutput{
        Sketch sketchRef;
        Sketch sketchQuery;
        long indexRef;
        long indexQuery;
        long pairCount;
        PairOutput pairs;
        public ContainOutput(Sketch sketchRefNew,Sketch sketchQueryNew, long indexRefNew,long indexQueryNew, long pairCountNew){
            sketchRef= sketchRefNew;
            sketchQuery= sketchQueryNew;
            indexRef = indexRefNew;
            indexQuery = indexQueryNew;
            pairCount = pairCountNew;
            pairs = new PairOutput(pairCount);
        
        }
       // pairs = new PairOutput(pairCount);
      
        }
     // ContainOutput
        class PairOutput{
            double score;
            double error;
            long pairCount; //da rivedere
            //costruttore forse
            public PairOutput(long pairCount){
                this.pairCount = pairCount;
            }
        }

      public ContainOutput contain(ContainInput input){
           Sketch sketchRef = input.sketchRef;
           Sketch sketchQuery = input.sketchQuery;
           ContainOutput output = new ContainOutput(input.sketchRef,input.sketchQuery,input.indexRef,input.indexQuery,input.pairCount);
           long i = input.indexQuery;
           long j = input.indexRef;
           for(long k=0;k<input.pairCount && i< sketchQuery.getReferenceCount();k++){
               output.pairs[k].score = containSketches(sketchRef.getReference(j).hashesSorted),sketchQuery.getReference(i).hashesSorted,output.pairs[k].error);
               j++;
               if(j == sketchRef.getReferenceCount()){
                   j=0;
                   i++;
               }
           
           }
           return output;
            
        } //contain       
                
    public int run() throws FileNotFoundException{
        if ( arguments.size() < 2 || options.get("help").active ){ //ereditati da Command{
        print();
        return 0;
        }
        int threads = (int) options.get("threads").getArgumentAsNumber();
        boolean list = options.get("list").active;
        Parameters parameters;
        parameters.error = options.get("errorThreshold").getArgumentAsNumber();
        if(sketchParameterSetup(parameters,(Command) this)) return 1;
        Sketch sketchRef;
        String fileReference = arguments.get(0); 
        boolean isSketch = hasSuffix(fileReference,suffixShetch);
        if(isSketch){
            if(options.get("kmer").active){
                System.err.println("ERROR: The option"+ options.get("kmer").identifier+ " cannot be used when a sketch is provided; it is inherited from the sketch ");
                return 1;
            }
            if(options.get("noncanonical").active){
                System.err.println("ERROR: The option"+ options.get("noncanonical").identifier+ " cannot be used when a sketch is provided; it is inherited from the sketch ");
                return 1;
            }
            if(options.get("protein").active){
                System.err.println("ERROR: The option"+ options.get("protein").identifier+ " cannot be used when a sketch is provided; it is inherited from the sketch ");
                return 1;
            }
            if(options.get("alphabet").active){
                System.err.println("ERROR: The option"+ options.get("alphabet").identifier+ " cannot be used when a sketch is provided; it is inherited from the sketch ");
            }
            else
                System.err.println("Sketching " + fileReference +" (provide sketch file made with \\\"mash sketch\\\" to skip)..." );
            
            ArrayList<String> refArgVector = new ArrayList<>();
            refArgVector.add(fileReference);
            sketchRef.initFromFiles(refArgVector,parameters);
            
            if(isSketch){
                parameters.minHashesPerWindow = sketchRef.getMinHashesPerWindow();
                parameters.kmerSize = sketchRef.getKmerSize();
                parameters.noncanonical = sketchRef.getNoncanonical();
                parameters.preserveCase = sketchRef.getPreserveCase();
                parameters.seed = sketchRef.getHashSeed();
                String alphabet = null;
                sketchRef.getAlphabetAsString(alphabet);
                setAlphabetFromString(parameters, alphabet.c_str());
            }else
                System.err.println("done.\n");
        }
        ThreadPoolExecutor threadPool =  (ThreadPoolExecutor) Executors.newFixedThreadPool(threads);
        threadPool.submit(()->{
           contain(ContainInput input); //far vedere al prof 
        });
        ArrayList<String> queryFiles = new ArrayList<String>();
        for ( int i = 1; i < arguments.size(); i++ )
        {
             if ( list )
            {
                splitFile(arguments.get(i), queryFiles);
            }
            else
            {
            queryFiles.add(arguments.get(i));
            }
        }
    
        Sketch sketchQuery;
        sketchQuery.initFromFiles(queryFiles, parameters, 0, true, true);

        long pairCount = sketchRef.getReferenceCount() * sketchQuery.getReferenceCount();
        long pairsPerThread = pairCount / parameters.parallelism;

            if ( pairsPerThread == 0 )
            {
                pairsPerThread = 1;
            }

           final long maxPairsPerThread = 0x1000;

            if ( pairsPerThread > maxPairsPerThread )
            {
                pairsPerThread = maxPairsPerThread;
            }

            long iFloor = pairsPerThread / sketchRef.getReferenceCount();
            long iMod = pairsPerThread % sketchRef.getReferenceCount();

            for ( long i = 0, j = 0; i < sketchQuery.getReferenceCount(); i += iFloor, j += iMod )
            {
                if ( j >= sketchRef.getReferenceCount() )
                {
                    if ( i == sketchQuery.getReferenceCount() - 1 )
                    {
                        break;
                    }

                    i++;
                    j -= sketchRef.getReferenceCount();
                }

                        threadPool.execute((Runnable)new ContainInput(sketchRef, sketchQuery, j, i, pairsPerThread, parameters));

                        while ( threadPool.awaitTermination(iMod, TimeUnit.DAYS) )
                        {
                                writeOutput(threadPool.popOutputWhenAvailable(), parameters.error);
                        }
            }

            while (! threadPool.isTerminated())
            {
                writeOutput(threadPool.popOutputWhenAvailable(), parameters.error);
            }
    
        return 0;
    
    }
    private void writeOutput(ContainOutput output, float error){
        long i = output.indexQuery;
        long j = output.indexRef;
        for(long k =0;k<output.pairCount && i< output.sketchQuery.getReferenceCount();k++){
            PairOutput pair = output.pairs[k];
            if(pair.error <= error){
                System.out.println(pair.score +"\t" + pair.error+ "\t"+ output.sketchRef.getReference(j).name+ "\t"+ output.sketchQuery.getReference(i).name);
            }
            j++;
            if(j == output.sketchRef.getReferenceCount()){
                j=0;
                i++;
            }
        }
        output=null;
    }// writeOutput
    public double containSketches(ArrayList hashesSortedRef,ArrayList hashesSortedQuery,double errorToSet){
        int common=0;
        int denom = hashesSortedRef.size() < hashesSortedQuery.size() ?
                hashesSortedRef.size() :
                hashesSortedQuery.size();
        int i=0;
        int j=0;
        for(int steps =0;steps < denom && i < hashesSortedRef.size(); steps++){
            if(hashLessThan(hashesSortedRef.get(i), hashesSortedQuery.get(j), hashesSortedRef.get64())){
                i++;
                steps --;
            }
            else if(hassLessThan(hashesSortedQuery.get(j),hashesSortedRef.get(i),hashesSortedRef.get64()))
                j++;
            
            else {
                i++;
                j++;
                common++;
            }
        }
        errorToSet = 1. /Math.sqrt(j);
        return (double) common /j;
    } // fine containSketches
} 
