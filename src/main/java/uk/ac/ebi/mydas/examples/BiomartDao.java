package uk.ac.ebi.mydas.examples;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringWriter;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


import org.apache.http.HttpEntity;
import org.apache.http.HttpHost;
import org.apache.http.HttpResponse; import org.apache.http.NameValuePair; import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.HttpClient; import org.apache.http.client.entity.UrlEncodedFormEntity; import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost; import org.apache.http.conn.params.ConnRoutePNames;
import org.apache.http.impl.client.DefaultHttpClient; import org.apache.http.message.BasicNameValuePair;

import uk.ac.ebi.mydas.exceptions.DataSourceException;
import uk.ac.ebi.mydas.model.DasFeature;
import uk.ac.ebi.mydas.model.DasFeatureOrientation;
import uk.ac.ebi.mydas.model.DasMethod;
import uk.ac.ebi.mydas.model.DasPhase;
import uk.ac.ebi.mydas.model.DasType;



public class BiomartDao {
	private static final String xmlDocHead="<?xml version=\"1.0\" encoding=\"utf-8\"?><!DOCTYPE Query><Query virtualSchemaName=\"default\" datasetConfigVersion=\"0.6\" uniqueRows=\"1\">";
	private static final String queryEnd="</Query>";
	private final String bioMartRootUrl;
	
		
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		BiomartDao bot = new BiomartDao("http://www.i-dcc.org/biomart/martservice");
		String chromosome="4";
		int homology_arm_start_lte=136637716;//equates to segment end coord
		int homology_arm_end_gte=136276059;//equates to segment start coord
		bot.getFeaturesMaps(chromosome, homology_arm_end_gte, homology_arm_start_lte);
//		HttpResponse response=bot.postData(bot.getFeaturesQuery(chromosome, homology_arm_start_lte, homology_arm_end_gte));//"/Users/jwarren/Documents/PERL/query.xml");
//		System.out.println("-----------------------------------------------");
//		HttpResponse response2=bot.postData(bot.getGeneNamesForProjectIdQuery(92695));
		
	}
	/**
	 * Constructor takes the biomart service url e.g. http://www.i-dcc.org/biomart/martservice
	 * @param biomartServiceUrl
	 */
	public BiomartDao(String biomartServiceUrl){
		this.bioMartRootUrl = biomartServiceUrl;
		
	}
	
	private HttpResponse postData(String xmlQuery) {
		
		HttpHost proxy = new HttpHost("wwwcache.ebi.ac.uk", 8080);

	    // Create a new HttpClient and Post Header
	    HttpClient httpclient = new DefaultHttpClient();
	    HttpPost httppost = new HttpPost(bioMartRootUrl);
	    HttpResponse response=null;
	    
	    try {
	    	httpclient.getParams().setParameter(ConnRoutePNames.DEFAULT_PROXY, proxy);

	        // Add your data
	        List<NameValuePair> nameValuePairs = new ArrayList<NameValuePair>(2);
	        nameValuePairs.add(new BasicNameValuePair("query", xmlQuery));
	        httppost.setEntity(new UrlEncodedFormEntity(nameValuePairs));
	        
	        // Execute HTTP Post Request
	        System.out.println("executing request to " + httppost.getURI()+ " via " + proxy);
	        response = httpclient.execute(httppost);
    
	    } catch (ClientProtocolException e) {
	        e.printStackTrace();
	    } catch (IOException e) {
	       e.printStackTrace();
	    } finally {
            // When HttpClient instance is no longer needed,
            // shut down the connection manager to ensure
            // immediate deallocation of all system resources
            httpclient.getConnectionManager().shutdown();
        }
	    return response;
	} 
	
	
	
	private String getFeaturesQuery(String chromosome, int segStart, int segStop){
		//add fudge factor here as below a certain point no features returned
		//int fudge=5000;
		//segStart=segStart-fudge;
		//segStop=segStop+fudge;
		StringBuilder sb=new StringBuilder();
		sb.append(xmlDocHead);
		sb.append("<Dataset name=\"idcc_targ_rep\" interface=\"default\">");
		//"allele_id"
		sb.append("<Attribute name=\"pipeline\" />");
		sb.append("<Attribute name=\"mgi_accession_id\" />");
		sb.append("<Attribute name=\"ikmc_project_id\" /><Attribute name=\"design_type\" /><Attribute name=\"design_subtype\" /><Attribute name=\"cassette\" /><Attribute name=\"allele_symbol_superscript\" /><Attribute name=\"homology_arm_start\" /><Attribute name=\"homology_arm_end\" /><Attribute name=\"loxp_start\" /><Attribute name=\"loxp_end\" /><Attribute name=\"cassette_start\" /><Attribute name=\"cassette_end\" /><Attribute name=\"strand\" />");
		sb.append("<Attribute name=\"allele_id\" />");
		sb.append("<Filter name=\"homology_arm_start_lte\" value=\""+segStop+"\" />");
		sb.append("<Filter name=\"homology_arm_end_gte\" value=\""+segStart+"\" />");
		//sb.append("<Filter name=\"homology_arm_start\" value=\""+segStop+"\" />");
		//sb.append("<Filter name=\"homology_arm_end\" value=\""+segStart+"\" />");
		sb.append("<Filter name=\"chromosome\" value=\""+chromosome+"\" />");
		sb.append("</Dataset>");
		sb.append(queryEnd);
		return sb.toString();
	}
	/**
	 * query to get the gene name for a ikmcProjectId - can only call this one at a time for each id
	 * @return
	 */
	private String getGeneNamesForProjectIdQuery(String ikmcProjId){
		StringBuilder sb=new StringBuilder();
		sb.append(xmlDocHead);
		sb.append(getGeneString(ikmcProjId));
		//multiple data sets doesn't seem to work with filters or multiple filters don't work either in one dataset - we just get the last filter i.e. result for the last ikmc proj id!!!!???
		//sb.append(getGeneNames(35941));
		//sb.append(getGeneNames(84587));
		sb.append(queryEnd);
		//System.out.println("query="+sb.toString()+"queryend");
		return sb.toString();
	}
	
	private String getGeneString(String ikmcProjId){
		return "<Dataset name=\"dcc\" interface=\"default\"><Attribute name=\"status\" /><Attribute name=\"marker_symbol\" /><Filter name=\"ikmc_project_id\" value=\""+ikmcProjId+"\" /></Dataset>";
		
		
	}
	
	protected List<Map<String, String>> getFeaturesMaps(String segId, int segmentStart, int segmentStop) {
		List<Map<String,String>> listOfMaps=new ArrayList<Map<String,String>>();
		
		// send a request to the url
		System.out.println("segStart="+segmentStart+" segmentStop="+segmentStop);
		int diff=segmentStop-segmentStart;
		System.out.println("difference="+diff);
		String xmlQuery=this.getFeaturesQuery(segId, segmentStart, segmentStop);
		HttpResponse response=this.postData(xmlQuery);
        
        HttpEntity responseEntity = response.getEntity() ;
		
		if (responseEntity != null) {
			
			InputStream stream = null;
			try {
				stream = responseEntity.getContent();
			} catch (IllegalStateException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			InputStreamReader in= new InputStreamReader(stream);
			  BufferedReader bin= new BufferedReader(in);
			  String text="";
			try {
				while ((text = bin.readLine()) != null)   {
					  // Print the content on the console
					Map<String,String> map=new HashMap<String,String>();
					  System.out.println (text);
					  String []fields=text.split("\t");
					  //System.out.println("homology_arm_end_gte="+map.get("homology_arm_end_gte")+" homology_arm_start_lte="+map.get("homology_arm_start_lte"));
					  //System.out.println("field1="+fields[1]);
					  map.put("pipeline", fields[0]);
					  map.put("mgi_accession_id",fields[1]);
					  map.put("ikmc_project_id",fields[2]);
					  map.put("design_type",fields[3]);
					  map.put("design_subtype", fields[4]);
					  map.put("cassette",fields[5]);
					  map.put("allele_symbol_superscript",fields[6]);
					  //System.out.println("allele_symbol_superscript"+fields[6]);
					  map.put("homology_arm_start", fields[7]);
					  map.put("homology_arm_end",fields[8]);
					  map.put("loxp_start",fields[9]);
					  map.put("loxp_end",fields[10]);
					  map.put("cassette_start",fields[11]);
					  map.put("cassette_end", fields[12]);
					  map.put("strand",fields[13]);
					  map.put("allele_id",fields[14]);
					  //map.put("homology_arm_start_lte",fields[15]);
					  
					listOfMaps.add(map);
				}
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			  
			
		}

		return listOfMaps;
	}
	
	protected  Map<String, String> getGeneNamesMapFromIkmcProjId(String ikmcProjId) {
		
		Map<String,String> geneNames=new HashMap<String,String>();
		// send a request to the url
		String xmlQuery=this.getGeneNamesForProjectIdQuery(ikmcProjId);
				HttpResponse response=this.postData(xmlQuery);
        
        HttpEntity responseEntity = response.getEntity() ;
		
		if (responseEntity != null) {
			
			InputStream stream = null;
			try {
				stream = responseEntity.getContent();
			} catch (IllegalStateException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			InputStreamReader in= new InputStreamReader(stream);
			  BufferedReader bin= new BufferedReader(in);
			  String text="";
			try {
				while ((text = bin.readLine()) != null)   {
					  // Print the content on the console
					String []geneInfo=text.split("\t");
					//status\" /><Attribute name=\"marker_symbol
					geneNames.put("status", geneInfo[0]);
					geneNames.put("marker_symbol", geneInfo[1]);
				}
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			  
			
		}

		return geneNames;
	}

}
