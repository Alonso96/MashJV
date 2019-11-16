package mash;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import javafx.util.Pair;

public class Command  {
    enum Type{Boolean,Number,Integer,Size,File,Stringa};
        class Option {
            
            Type type;
            public String category;
            public String identifier;
            public String description;
            public String argument;
            public String argumentDefault;
            float argumentAsNumber;
            float argumentMin;
            float argumentMax;
            boolean active;
            boolean changed;
        public Option(Type type,String identifiernew,String categorynew,String descriptionnew, String argumentDefaultnew,float argumentMinnew,float argumentMaxnew){
            this.type = type;
            identifier = identifiernew;
            category = categorynew;
            description = descriptionnew;
            argumentDefault = argumentDefaultnew;
            argumentMin = argumentMinnew;
            argumentMax = argumentMaxnew;
            this.active = false;
        
    
   
    }
    
    public float getArgumentAsNumber(){
        return this.argumentAsNumber;
    }
    public void setArgument(String argumentnew){
    	this.argument = argumentnew;
		long a = (long) argumentAsNumber;
		if(type == Type.Number || type == Type.Integer) { // primo if
			if (argument.length() == 0) {
				argumentAsNumber = 0;
				return;
			}
			boolean failed = false;
			try {
				argumentAsNumber = Float.parseFloat(argument);
				if(argumentMin != argumentMax && (argumentAsNumber < argumentMin || argumentAsNumber > argumentMax)) {
					failed = true;
				}
				else if(type == Type.Integer && (a != argumentAsNumber))
					failed = true;
			}catch(Exception e) {
				
				failed= true;
				e.printStackTrace();
			}
			if (failed) {
				System.err.println("ERROR: Argument must be an integer");
				if(argumentMin != argumentMax) {
					System.err.println("between "+ argumentMin + "and " + argumentMax);
				}
				System.err.println("("+ argument +" given)");
				System.exit(1);
			}
		} // end primo if
		else if (type == Type.Size) {
			if(argument.length()==0) {
				argumentAsNumber=0;
				return;
			}
			char suffix = argument.charAt(argument.length() -1);
			long factor = 1;
			if ( suffix < '0' || suffix > '9' )
	    	{
	    		switch ( suffix )
	    		{
	    			case 'k':
	    			case 'K':
	    				factor = 1000;
	    				break;
	    			case 'm':
	    			case 'M':
	    				factor = 1000000;
	    				break;
	    			case 'g':
	    			case 'G':
	    				factor = 1000000000;
	    				break;
	    			case 't':
	    			case 'T':
	    				factor = 1000000000000L;
	    				break;
	    			default:
						System.err.println("ERROR: Unrecognized unit (\"" + suffix + "\") in argument to -" + identifier + ". If specified, unit must be one of [kKmMgGtT].");
						System.exit(1);
	    		}
	    		
	    		argument = argument.substring(0,argument.length()-1);
	    	}
			boolean fail = false;
			long b = (long) argumentAsNumber;
			try {
				argumentAsNumber = Float.parseFloat(argument);
			}catch(Exception e) {
				fail = true;
			}
			if ((argumentAsNumber <= 0) ||(b != argumentAsNumber)){
				fail = true;
			}
			if (fail) {
				System.err.println("ERROR: Argument to  must be a whole number, optionally followed by one of [kKmMgGtT])");
				System.exit(1);
			}
			argumentAsNumber *= factor;
		}
    }
}
    
	protected ArrayList<String> categoryNames = new ArrayList<String>();
	protected Map<String,Option> options = new HashMap<String,Option>();
	protected Map<String,Option> optionsAvaible = new HashMap<String,Option>();
	protected ArrayList<String> arguments = new ArrayList<String>();
	private Map<String,String> optionsNamesByIdentifier = new HashMap<String,String>();
	private Map<String,ArrayList<String>> optionsNamesByCategory = new HashMap<String,ArrayList<String>>();
	private ArrayList<String> categories = new ArrayList<String>();
	private Map<String,String> categoryDisplayNames = new HashMap<String,String>();
	String name;
	String summary;
	String description;
	String argumentString;
	
	public Command(){
            addAvailableOption("help", new Option(Type.Boolean, "h", "", "Help", "",0,0));
            addAvailableOption("kmer", new Option(Type.Integer, "k", "Sketch", "K-mer size. Hashes will be based on strings of this many nucleotides. Canonical nucleotides are used by default (see Alphabet options below).", "21", 1, 32));
            addAvailableOption("windowed", new Option(Type.Boolean, "W", "Sketch", "Windowed", "",0,0));
            addAvailableOption("window", new Option(Type.Integer, "L", "Window", "Window length. Hashes that are minima in any window of this size will be stored.", "10000",0,0));
            //addAvailableOption("error", new Option(Type.Number, "e", "Sketch", "Error bound. The (maximum) number of min-hashes in each sketch will be one divided by this number squared.", "0.05"));
            addAvailableOption("sketchSize", new Option(Type.Integer, "s", "Sketch", "Sketch size. Each sketch will have at most this many non-redundant min-hashes.", "1000",0,0));
            addAvailableOption("verbose", new Option(Type.Boolean, "v", "Output", "Verbose", "",0,0));
            addAvailableOption("silent", new Option(Type.Boolean, "s", "Output", "Silent", "",0,0));
            addAvailableOption("individual", new Option(Type.Boolean, "i", "Sketch", "Sketch individual sequences, rather than whole files, e.g. for multi-fastas of single-chromosome genomes or pair-wise gene comparisons.", "",0,0));
            addAvailableOption("warning", new Option(Type.Number, "w", "Sketch", "Probability threshold for warning about low k-mer size.", "0.01", 0, 1));
            addAvailableOption("reads", new Option(Type.Boolean, "r", "Sketch", "Input is a read set. See Reads options below. Incompatible with -i.", "",0,0));
            addAvailableOption("seed", new Option(Type.Integer, "S", "Sketch", "Seed to provide to the hash function.", "42", 0, 0xFFFFFFFF));
            addAvailableOption("memory", new Option(Type.Size, "b", "Reads", "Use a Bloom filter of this size (raw bytes or with K/M/G/T) to filter out unique k-mers. This is useful if exact filtering with -m uses too much memory. However, some unique k-mers may pass erroneously, and copies cannot be counted beyond 2. Implies -r.","",0,0));
            addAvailableOption("minCov", new Option(Type.Integer, "m", "Reads", "Minimum copies of each k-mer required to pass noise filter for reads. Implies -r.", "1",0,0));
            addAvailableOption("targetCov", new Option(Type.Number, "c", "Reads", "Target coverage. Sketching will conclude if this coverage is reached before the end of the input file (estimated by average k-mer multiplicity). Implies -r.","",0,0));
            addAvailableOption("genome", new Option(Type.Size, "g", "Reads", "Genome size (raw bases or with K/M/G/T). If specified, will be used for p-value calculation instead of an estimated size from k-mer content. Implies -r.","",0,0));
            addAvailableOption("noncanonical", new Option(Type.Boolean, "n", "Alphabet", "Preserve strand (by default, strand is ignored by using canonical DNA k-mers, which are alphabetical minima of forward-reverse pairs). Implied if an alphabet is specified with -a or -z.", "",0,0));
            addAvailableOption("protein", new Option(Type.Boolean, "a", "Alphabet", "Use amino acid alphabet (A-Z, except BJOUXZ). Implies -n, -k 9.", "",0,0));
            addAvailableOption("alphabet", new Option(Type.Stringa, "z", "Alphabet", "Alphabet to base hashes on (case ignored by default; see -Z). K-mers with other characters will be ignored. Implies -n.", "",0,0));
            addAvailableOption("case", new Option(Type.Boolean, "Z", "Alphabet", "Preserve case in k-mers and alphabet (case is ignored by default). Sequence letters whose case is not in the current alphabet will be skipped when sketching.", "",0,0));
            addAvailableOption("threads", new Option(Type.Integer, "p", "", "Parallelism. This many threads will be spawned for processing.", "1",0,0));
            addAvailableOption("pacbio", new Option(Type.Boolean, "pacbio", "", "Use default settings for PacBio sequences.", "",0,0));
            addAvailableOption("illumina", new Option(Type.Boolean, "illumina", "", "Use default settings for Illumina sequences.", "",0,0));
            addAvailableOption("nanopore", new Option(Type.Boolean, "nanopore", "", "Use default settings for Oxford Nanopore sequences.", "",0,0));
            addAvailableOption("factor", new Option(Type.Number, "f", "Window", "Compression factor", "100",0,0));

            addCategory("", "");
            addCategory("Input", "Input");
            addCategory("Output", "Output");
            addCategory("Sketch", "Sketching");
            addCategory("Window", "Sketching (windowed)");
            addCategory("Reads", "Sketching (reads)");
            addCategory("Alphabet", "Sketching (alphabet)");
			}
	
	/*public void setArgument(String argumentnew) {
		this.argument = argumentnew;
		long a = (long) argumentAsNumber;
		if(type == Type.Number || type == Type.Integer) { // primo if
			if (argument.length() == 0) {
				argumentAsNumber = 0;
				return;
			}
			boolean failed = false;
			try {
				argumentAsNumber = Float.parseFloat(argument);
				if(argumentMin != argumentMax && (argumentAsNumber < argumentMin || argumentAsNumber > argumentMax)) {
					failed = true;
				}
				else if(type == Type.Integer && (a != argumentAsNumber))
					failed = true;
			}catch(Exception e) {
				
				failed= true;
				e.printStackTrace();
			}
			if (failed) {
				System.err.println("ERROR: Argument must be an integer");
				if(argumentMin != argumentMax) {
					System.err.println("between "+ argumentMin + "and " + argumentMax);
				}
				System.err.println("("+ argument +" given)");
				System.exit(1);
			}
		} // end primo if
		else if (type == Type.Size) {
			if(argument.length()==0) {
				argumentAsNumber=0;
				return;
			}
			char suffix = argument.charAt(argument.length() -1);
			long factor = 1;
			if ( suffix < '0' || suffix > '9' )
	    	{
	    		switch ( suffix )
	    		{
	    			case 'k':
	    			case 'K':
	    				factor = 1000;
	    				break;
	    			case 'm':
	    			case 'M':
	    				factor = 1000000;
	    				break;
	    			case 'g':
	    			case 'G':
	    				factor = 1000000000;
	    				break;
	    			case 't':
	    			case 'T':
	    				factor = 1000000000000L;
	    				break;
	    			default:
						System.err.println("ERROR: Unrecognized unit (\"" + suffix + "\") in argument to -" + identifier + ". If specified, unit must be one of [kKmMgGtT].");
						System.exit(1);
	    		}
	    		
	    		argument = argument.substring(0,argument.length()-1);
	    	}
			boolean fail = false;
			long b = (long) argumentAsNumber;
			try {
				argumentAsNumber = Float.parseFloat(argument);
			}catch(Exception e) {
				fail = true;
			}
			if ((argumentAsNumber <= 0) ||(b != argumentAsNumber)){
				fail = true;
			}
			if (fail) {
				System.err.println("ERROR: Argument to  must be a whole number, optionally followed by one of [kKmMgGtT])");
				System.exit(1);
			}
			argumentAsNumber *= factor;
		}
			
		} // fine setArgument
		*/
		public void addOption(String name, Option option) {
			options.put(name, option);
			categoryNames.add(name); // array list di nomi categorie
			optionsNamesByCategory.put(option.category,categoryNames);
			optionsNamesByIdentifier.put(option.identifier,name);
		}
		public void print() {
			//String[] columns ;
			ArrayList<String> columns = new ArrayList<String>();
			columns.add(0,"mash " + name + " [options] " + argumentString);
			//columns[0]= "mash " + name + " [options] " + argumentString;
			System.out.println(columns); 
			System.out.println("Description");
			columns.remove(0);
			columns.add(0,description);
			//columns[0] = "";
		//	columns[0]= description;
			System.out.println(columns);
			if (options.size() == 0) return;
			System.out.println("Options");
			columns.clear();
			columns.add(0,"Option");
			columns.add(1,"Description (range) [default]");
			ArrayList< Pair<Integer,String>> dividers = new ArrayList<Pair<Integer,String>>();
			for(String a : categories) {
				if(a != null && (optionsNamesByCategory.get(a).size()!=0)) 
				//	Pair<Integer,String> pair = new Pair<Integer,String>(columns.get(0).length(), categoryDisplayNames.get(a)+ "...");
				
					dividers.add(new Pair<Integer,String>(columns.get(0).length(), categoryDisplayNames.get(a)+ "..."));
				
				for(Iterator<Entry<String, ArrayList<String>>> j = optionsNamesByCategory.entrySet().iterator(); j!= optionsNamesByCategory.get(a); j.next() ) {
					Option option = options.get(j);
					String optionString = "-" + option.identifier;
					String range;
					if (option.type != Type.Boolean) {
						String type;
						switch(option.type) {
						case Boolean:
							break;
						case Number:
							type = "num";
							break;
						case Integer:
							type = "int";
							break;
						case Size:
							type = "size";
							break;
						case File:
							type = "path";
							break;
						case Stringa:
							type = "text";
							break;
						}//fine switch
						optionString+= "<"+ option.type+">";
					}//fine if
					columns.add(0,optionString);
					String descString= option.description;
					if(option.argumentMin!= option.argumentMax) {
						float stringMin;
						float stringMax;
						if(option.type== Type.Integer) {
							stringMin = (long) option.argumentMin;
							stringMax = (long) option.argumentMax;
							/* 
							 * Capire come implementare questo stream e assegnarci un long nonostante sia uno string
							 */
							
						}else {
							stringMin = option.argumentMin;
							stringMax = option.argumentMax;
						}
						descString+="("+Float.toString(stringMin)+ "-"+ Float.toString(stringMax)+ ")";
					}// fine if(option.argumentMin!= option.argumentMax
					if(option.argumentDefault != "") descString+= "["+ option.argumentDefault + "]";
					columns.add(1,descString);
				}// fine 2 for
			}// fine for(String a:categories)
			System.out.println(columns+" " +dividers);
		} //fine print
		
		public int run(int argc, String [] args ) {
			
				
			  for ( int i = 0; i < argc; i++ )
			    {
			        if ( args[i].charAt(0) == '-' && args[i].charAt(1) != 0 )
			        {
			            if ( Collections.frequency((Collection<?>) optionsNamesByIdentifier,  args[i]+1) == 0 )
			            {
			                System.err.println( "ERROR: Unrecognized option: " +args[i]);
			                return 1;
			            }
			            
			            Option  option = options.get(optionsNamesByIdentifier.get(args[i] + 1));
			            
			            option.active = true;
			            
			            if ( option.type != Type.Boolean )
			            {
			                i++;
			                
			                if ( i == argc )
			                {
			                   System.err.println("ERROR: -" +option.identifier +" requires an argument");
			                    return 1;
			                }
			                
			                option.setArgument(args[i]);
			            }
			        }
			        else
			        {
			            arguments.add(args[i]);
			        }
			    }
			    
			    return run(argc, args); // ???
		}
		public void useOption(String name) {
			addOption(name, optionsAvaible.get(name));
		}
		
		public void useSketchOptions() {
		    useOption("threads");
		    useOption("kmer");
		    useOption("noncanonical");
		    useOption("protein");
		    useOption("alphabet");
		    useOption("case");
		   // #ifdef COMMAND_FIND
		    useOption("windowed");
		    useOption("window");
		    useOption("factor");
	//	#endif
		    useOption("sketchSize");
		    useOption("individual");
		    useOption("seed");
		    useOption("warning");
		    useOption("reads");
		    useOption("memory");
		    useOption("minCov");
		    useOption("targetCov");
		    useOption("genome");
		}
		
		public void addAvailableOption(String name, Option option) {
			optionsAvaible.put(name, option);
		}//addAvailableOption
		
		public void addCategory(String name, String displayName) {
			if (Collections.frequency((Collection<?>)  categoryDisplayNames,name) == 0 )
			{
				categories.add(name);
			    categoryDisplayNames.put(name,displayName);
			    optionsNamesByCategory.put(name,new ArrayList<String>());
			}
		} //addCategory
		
		public void splitFile(String file, ArrayList<String> lines) throws FileNotFoundException {
			String line = null;
			BufferedReader reader = new BufferedReader(new FileReader(file));
			while(line!=null) {
			    
			     try {
					line = reader.readLine();
					lines.add(line);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					System.err.println("ERROR: Could not open" + file + "\n");
					e.printStackTrace();
				}
			}
		}//splitFile
		
	/**	public void printColumns(ArrayList<String> columns,int indent, int spacing, char  missing, int max) {
			printColumns(columns, new ArrayList<Pair<Integer,String>>() , indent, spacing, missing, max);
		}
		public void printColumns(ArrayList<String>columns,ArrayList<Pair<Integer,String>>dividers, int indent, int spacing, char missing, int max) {
			
		}
		*/
	}
	
