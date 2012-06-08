package uk.ac.ebi.mydas.examples;

import java.net.MalformedURLException;
import java.net.URL;

public class IkmcLinks {

	protected static URL getCassetteImageLink(String alleleId) {
		return getUrl("http://www.knockoutmouse.org/targ_rep/alleles/" + alleleId
				+ "/allele-image");
	}

	protected static URL getGbFileforEsCellClone(String alleleId) {
		
		return getUrl("http://www.knockoutmouse.org/targ_rep/alleles/" + alleleId
				+ "/escell-clone-genbank-file");
	}
	
	protected static URL getMgiDetailsLink(String mgiAccession){
		String urlString="http://www.knockoutmouse.org/genedetails/"+mgiAccession;
		return getUrl(urlString);
	}
	
	private static URL getUrl(String urlString){
		URL url = null;
		try {
			url = new URL(urlString);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return url;
	}
}
