package mash;
import mash.MurmurHash3;
public class Hash_U {
	long hash32;
	long hash64;
	
	
	public long getHash(String seq, int length, int seed, boolean use64) {
	//	char[] data = new char[16];
	   MurmurHash3 m ;
	   //implementare funzione MarmurHash3 per questo tipo di variabilig
	   MurmurHash3.LongPair result = new MurmurHash3.LongPair(); // da verificare
	   MurmurHash3.murmurhash3_x64_128(seq.getBytes(), length, seed, result);
	   Hash_U hash = new Hash_U();
	    
	    if ( use64 ) // si suppone sia sempre 64bit
	    {
	        hash.hash64 = result.val1;
	    }
	    else
	    {
	        hash.hash32 =  result.val2;
	    
	   
	}
	    return hash64;
	}
	
	/***  char data[16];
    MurmurHash3_x64_128(seq, length, seed, data);
#endif
    
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = *((hash64_t *)data);
    }
    else
    {
        hash.hash32 = *((hash32_t *)data);
    }
    
    return hash;
    ***/
	public boolean hashLessThan(Hash_U hash1, Hash_U hash2, boolean use64) {
		 if ( use64 )
		    {
		        return hash1.hash64 < hash2.hash64;
		    }
		    else
		    {
		        return hash1.hash32 < hash2.hash32;
		    }
	}

}
