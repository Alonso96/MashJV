package mash;

import mash.Option.Type;

public class Option {
    enum Type{Boolean,Number,Integer,Size,File,Stringa};
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
    public Option(Type type,String identifierNew,String categoryNew,String descriptionNew, String argumentDefaultNew,float argumentMinNew,float argumentMaxNew){
    	this.type = type;
    	identifier = identifierNew;
    	category = categoryNew;
    	description = descriptionNew;
    	argumentDefault = argumentDefaultNew;
    	argumentMin = argumentMinNew;
    	argumentMax = argumentMaxNew;
    	this.active = false;
        
    
   
    }
    
    public float getArgumentAsNumber(){
        return this.argumentAsNumber;
    }
    public void setArgument(String argumentNew){
    	this.argument = argumentNew;
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

