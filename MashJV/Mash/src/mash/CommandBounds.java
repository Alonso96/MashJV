package mash;
import ch.visnet.jgsl.Jgsl;
import ch.visnet.jgsl.sf.

public class CommandBounds extends Command {
	
    public CommandBounds(){
    name = "bounds";
    summary = "Print a table of Mash error bounds.";
    description = "Print a table of Mash error bounds for various sketch sizes and Mash distances based on a given k-mer size and desired confidence. Note that these calculations assume sequences are much larger than the sketch size, and will overestimate error bounds if this is not the case.";
    argumentString = "";
    
    useOption("help");
    addOption("kmer", new Option(Option.Type.Integer, "k", "", "k-mer size.", "21", 1, 32));
    addOption("prob", new Option(Option.Type.Number, "p", "", "Mash distance estimates will be within the given error bounds with this probability.", "0.99", 0, 1));
   
    }
    
    public int run(){
    	ch.visnet.jgsl.Jgsl.init();
    	
        if (options.get("help").active){
            print();
            return 0;
        }
    final  int sketchSizeCount = 9;
    final  double [] sketchSizes = new double[]{100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    final double [] dists = new double[]{0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    int k = getOption("kmer").getArgumentAsNumber();
    double q2 = (1.0 - getOption("prob").getArgumentAsNumber()) / 2.0;
    System.out.println("Parameters (run with -h for details):");
    System.out.println("  k:  " + k);
    System.out.println("  p:  "+getOption("prob").getArgumentAsNumber());
    for(int cont =0;cont<2;cont++){
        if(cont !=0){
            System.out.println("\tScreen distance");
        } 
        else 
             System.out.println("\tMash Distance");
        System.out.println("Sketch");
        for (int i =0;i<distCount;i++){
            System.out.println("\t"+ dists[i]);
        }
        for(int i=0;i<sketchSizeCount;i++){
            int s = (int) sketchSizes[i];
            System.out.println(s);
        
             for ( int j = 0; j < distCount; j++ ){
                double m2j;
                
        if ( cont!=0 ){
                    //m2j = exp(-k * dists[j]);
                    m2j = Math.pow( (1.0 - dists[j]), k); // binomial model
                    }else{
                        m2j = 1.0 / (2.0 * Math.exp(k * dists[j]) - 1.0);
            }
            
            int x = 0;
            
                        while ( x < s ){
                            Jgsl.init();
                            double cdfx = cdf(binomial(s,m2j),x);
                            double cdfx =gsl_cdf_binomial_P(x,m2j,s);
                            if(cdfx > q2) break;
                            x++;
                        } // fine while
                        double je = ((double) x) / s;
                        double j2m;
                        if(cont !=0) j2m = 1.0 -Math.pow(je,1. /k );
                        else j2m = -1.0 / k * Math.log(2.0 * je / (1.0 + je));
                        System.out.println("\t"+ (j2m - dists[j]));
                        
                }
        }
    
    }
    return 0;
    }
}



