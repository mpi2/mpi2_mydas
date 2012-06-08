package uk.ac.ebi.mydas.examples;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.springframework.core.NestedRuntimeException;
import org.springframework.web.client.RestClientException;
import org.springframework.web.client.RestTemplate;
import org.springframework.web.util.UriUtils;

import uk.ac.ebi.mydas.examples.solr.*;
import uk.ac.ebi.mydas.examples.solr.Result.Doc;
import uk.ac.ebi.mydas.exceptions.BadReferenceObjectException;
import uk.ac.ebi.mydas.exceptions.DataSourceException;
import uk.ac.ebi.mydas.model.DasAnnotatedSegment;
import uk.ac.ebi.mydas.model.DasComponentFeature;
import uk.ac.ebi.mydas.model.DasEntryPoint;
import uk.ac.ebi.mydas.model.DasEntryPointOrientation;
import uk.ac.ebi.mydas.model.DasFeature;
import uk.ac.ebi.mydas.model.DasFeatureOrientation;
import uk.ac.ebi.mydas.model.DasMethod;
import uk.ac.ebi.mydas.model.DasPhase;
import uk.ac.ebi.mydas.model.DasSequence;
import uk.ac.ebi.mydas.model.DasType;
import org.apache.log4j.Logger;

public class IkmcAllelesManager {
	
	static Logger logger = Logger.getLogger(IkmcAllelesManager.class);
	private final BiomartDao biomart;
	private static final String xmlDocHead = "<?xml version=\"1.0\" encoding=\"utf-8\"?><!DOCTYPE Query><Query virtualSchemaName=\"default\" datasetConfigVersion=\"0.6\" uniqueRows=\"1\">";
	private static final String queryEnd = "</Query>";
	private final String bioMartRootUrl;
	
	/**
	 * List of the types used in this data source
	 */
	private ArrayList<DasType> types;
	/**
	 * Types to be used in the data source: chromosome, gene, transcript and
	 * exon
	 */
	
	/**
	 * As this data source just have one method, it can be defined as a
	 * parameter to facilitate its use
	 */
	private final DasMethod method = new DasMethod("MGI",
			"combinatorial evidence", "ECO:0000043");

	private String database = "genotype1";
	private Collection<DasEntryPoint> entryPoints = null;
	private Map<String,Map<String,String>> ikmcGeneMap;

	public IkmcAllelesManager(String bioMartUrl) throws DataSourceException {
		this.bioMartRootUrl = bioMartUrl;
		biomart = new BiomartDao(this.bioMartRootUrl);

		// Initialize types
		types = new ArrayList<DasType>();
		//types.add(geneType);
		ikmcGeneMap=new HashMap<String,Map<String,String>>();
	}

	

	public DasAnnotatedSegment getSubmodelBySegmentId(String segmentId,
			int start, int stop) throws DataSourceException,
			BadReferenceObjectException {
		return getFeatures(segmentId, start, stop, -1);
	}

	public DasAnnotatedSegment getFeatures(String segmentId, int start,
			int stop, int maxbins) throws DataSourceException,
			BadReferenceObjectException {
			ikmcGeneMap.clear();//clear the map of biomart request results for gene symbols etc
			//as we don't want to store these for all time just per request to save mutliple biomart requests per request
		

		Collection<DasFeature> features = getFeatures(segmentId, start, stop);
		DasAnnotatedSegment segment = null;
		try {
			segment = new DasAnnotatedSegment(segmentId, new Integer(start),
					new Integer(stop), "1.0", segmentId, features);
		} catch (DataSourceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return segment;
	}

	/**
	 *  			map.put("pipeline", fields[0]);
					  map.put("mgi_accession_id",fields[1]);
					  map.put("ikmc_project_id",fields[2]);
					  map.put("design_type",fields[3]);
					  map.put("design_subtype", fields[4]);
					  map.put("cassette",fields[5]);
					  map.put("allele_symbol_superscript",fields[6]);
					  map.put("homology_arm_start", fields[7]);
					  map.put("homology_arm_end",fields[8]);
					  map.put("loxp_start",fields[9]);
					  map.put("loxp_end",fields[10]);
					  map.put("cassette_start",fields[11]);
					  map.put("cassette_end", fields[12]);
					  map.put("strand",fields[13]);
	 * @param segId
	 * @param start
	 * @param stop
	 * @return
	 */
	private Collection<DasFeature> getFeatures(String segId, int start, int stop) {
		List<Map<String, String>> listOfMaps = biomart.getFeaturesMaps(segId,
				start, stop);
		System.out.println("biomart maps size="+listOfMaps.size());
		LinkedList<DasFeature> features = new LinkedList<DasFeature>();
		int i=0;
		for (Map<String, String> map : listOfMaps) {
			//if no allel superscript assume no allele and so don't include that group in the features response
			if(!map.get("allele_symbol_superscript").equals("")){
			  //Map<String, String> geneNames = biomart.getGeneNamesMapFromIkmcProjId(ikmcProjId);
			//logger.info("gene info:"+geneNames);
//			 my @features_defaults = (
//			            'group'     => $group,
//			            'grouptype' => $self->_get_grouptype( $feature_set ),
//			            'grouplink' => $self->_get_grouplink( $feature_set ),
//			            'groupnote' => $self->_get_groupnote( $feature_set ),
//			            'note'      => $self->_get_groupnote( $feature_set ),
//			            'link'      => $self->_get_grouplink( $feature_set ),
//			            'linktxt'   => $group,
//			            'ori'       => $feature_set->{strand}
//			        );
				System.out.println("allele_symbol_superscript="+map.get("allele_symbol_superscript"));
				
			String dGroup=this.getGroup(map);
			//System.out.println("group="+dGroup);
			String dGroupType=this.getGroupType(map);
			//System.out.println("group type="+dGroupType);
			Map dGroupLink=this.getGroupLink(map, dGroup);
			String dGroupNote=this.getGroupNote(map);
			String dNote=dGroupNote;
			Collection<String> dNotes=new ArrayList<String>();
			dNotes.add(dNote);
			Map dLink=dGroupLink;
			String dId=dGroup+ " - Homology Arm start";
			DasFeatureOrientation dOrientation = DasFeatureOrientation.ORIENTATION_SENSE_STRAND;
			if(map.get("strand").equals("-")){
				dOrientation=DasFeatureOrientation.ORIENTATION_ANTISENSE_STRAND;
			}
			
			//set up an array that holds just one parent and use for all subsequent features in this set
			ArrayList<String> parents=new ArrayList<String>();
			parents.add(dGroup);
			//set up the parts array to store the name of all the sub features of the 'group'
			ArrayList<String> parts=new ArrayList<String>();
	
			//how vvi builds features- although ours will be slightly different as we are using parent parts relationships!! so top level features will be groups..
//			        
//			        # Add G5 feature
//			        push @features, {
//			            @features_defaults,
//			            'id'    => $group . ' - Homology Arm start',
//			            'label' => 'Homology Arm start',
//			            'type'  => 'Homology Arm start',
//			            'start' => $feature_set->{homology_arm_start} - 50,
//			            'end'   => $feature_set->{homology_arm_start}
//			        };
//			        
			//just want to make a test feature group once
			if(i==0){
			List<DasFeature> testFeatures=this.testFeatures(map);
			features.addAll(testFeatures);
			}
			i++;
			//end of make test feature
			DasType dType = new DasType("Homology Arm start", "", "",
					"Homology Arm start");
			DasFeature featureG5 = null;
			try {
				
				featureG5 = new DasFeature(dId, "Homology Arm start", dType, method, Integer.valueOf(map.get("homology_arm_start"))-50,
						Integer.valueOf(map.get("homology_arm_start")), /*score*/null, dOrientation,
						DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/dNotes,/*links Map<URL,String>*/ dLink,/*targets*/ null,/*parents*/ parents,
						/*parts*/null);
			} catch (DataSourceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			features.add(featureG5);
			parts.add(dId);
//			        # Add U5-U3 feature
//			        if( $feature_set->{cassette_start} <= $ensembl_end 
//			            && $feature_set->{cassette_end} >= $ensembl_start){
//			            push @features, {
//			                @features_defaults,
//			                'id'    => $group . ' - Cassette insertion',
//			                'label' => "Cassette ($feature_set->{cassette})",
//			                'type'  => "Cassette ($feature_set->{cassette})",
//			                $self->_get_position( $feature_set, 'cassette' )
//			            };
//			        }
			DasFeature featureU5U3 = null;
			if(this.getPosition(map, "cassette").get("start")<= stop && this.getPosition(map, "cassette").get("end") >= start){
			DasType cassetteType=new DasType("Cassette ("+map.get("cassette")+")","","", "");
			String u5u3FeatureId=dGroup+" - Cassette insertion";
			try {
				
				featureU5U3 = new DasFeature(u5u3FeatureId, "Cassette ("+map.get("cassette")+")", cassetteType, method, this.getPosition(map, "cassette").get("start"),
						this.getPosition(map, "cassette").get("end"), /*score*/null, dOrientation,
						DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/dNotes,/*links Map<URL,String>*/ dLink,/*targets*/ null,/*parents*/parents,
						/*parts*/null);
			} catch (DataSourceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			features.add(featureU5U3);
			parts.add(u5u3FeatureId);
			}
//			        
//			        # Add D5-D3 feature
//			        if( $feature_set->{loxp_start} && $feature_set->{loxp_end} ){
//			            if( $feature_set->{loxp_start} <= $ensembl_end 
//			             && $feature_set->{loxp_end} >= $ensembl_start){
//			                push @features, {
//			                    @features_defaults,
//			                    'id'    => $group . ' - LoxP insertion',
//			                    'label' => 'LoxP',
//			                    'type'  => 'LoxP',
//			                    $self->_get_position( $feature_set, 'loxp' )
//			                };
//			            }
//			        }
			if(!map.get("loxp_start").equals("") && !map.get("loxp_end").equals("")){
				//if(Integer.valueOf(map.get("loxp_start"))<= stop && Integer.valueOf(map.get("loxp_end")) >= start){
			DasFeature featureD5U3 = null;
			//System.out.println("creating loxp feature");
			DasType loxpType=new DasType("LoxP","LoxP","", "LoxP");
			String d5u3Id=dGroup+" - LoxP insertion";
			try {
				featureD5U3 = new DasFeature(d5u3Id, "LoxP", loxpType, method, this.getPosition(map, "loxp").get("start"),
						this.getPosition(map, "loxp").get("end"), /*score*/null, dOrientation,
						DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/dNotes,/*links Map<URL,String>*/ dLink,/*targets*/ null,/*parents*/ parents,
						/*parts*/null);
			} catch (DataSourceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			features.add(featureD5U3);
			parts.add(d5u3Id);
				//}
			}
//			        
//			        # Add G3 feature
//			        push @features, {
//			            @features_defaults,
//			            'id'    => $group . ' - Homology Arm end',
//			            'label' => 'Homology Arm end',
//			            'type'  => 'Homology Arm end',
//			            'start' => $feature_set->{homology_arm_end},
//			            'end'   => $feature_set->{homology_arm_end} + 50
//			        };
//			    }
			    
			DasType homEndType = new DasType("Homology Arm end", "", "",
					"Homology Arm end");;
			DasFeature featureG3 = null;
			try {
				featureG3 = new DasFeature(dGroup+"- Homology Arm end", "Homology Arm end", homEndType, method, Integer.valueOf(map.get("homology_arm_end")),
						Integer.valueOf(map.get("homology_arm_end"))+50, /*score*/null, dOrientation,
						DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/dNotes,/*links Map<URL,String>*/ dLink,/*targets*/ null,/*parents*/ parents,
						/*parts*/null);
			} catch (DataSourceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			features.add(featureG3);
			parts.add(dId);
			
			
			//set up the group/parent id here which uses the default for groups for this set and has to be set up last so we know what the children are for the constructor:
			DasType defaultGroupType = new DasType(dGroupType, dGroupType, "",
					dGroupType);
			DasFeature group = null;
			
			try {
				group = new DasFeature(dGroup, null, defaultGroupType, method, this.getPosition(map, "homology_arm").get("start")-50,
						this.getPosition(map, "homology_arm").get("end"), /*score*/null, dOrientation,
						DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/dNotes,/*links Map<URL,String>*/ dLink,/*targets*/ null,/*parents*/ null,
						/*parts*/parts);
			} catch (DataSourceException e) {
				e.printStackTrace();
			}
			features.push(group);
			}
		}
		//for(DasFeature feature: features){
		//	System.out.println("fId="+feature.getFeatureId()+" "+feature.toString());
		//}
		return features;

	}

	public ArrayList<DasType> getTypes() {
		System.out.println("Getting types from genotype manager");
		for (DasType type : types) {
			System.out.println(type.getCvId());
		}
		return types;
	}

	public Integer getTotalCountForType(String typeId)
			throws DataSourceException {

		return null;
	}

	public String getDatabase() {
		return database;
	}
//viveks helper methods implement each :
//	sub _get_group {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    my $gene_symbol = $self->_get_gene_symbol( $feature_set ) || "";
//	    
//	    if ( $feature_set->{allele_symbol_superscript} ) {
//	        return "$gene_symbol $feature_set->{allele_symbol_superscript}";
//	    }
//	    
//	    if ( $feature_set->{pipeline} ) {
//	        return "$feature_set->{pipeline} $gene_symbol mutant allele";
//	    }
//	    
//	    return "$gene_symbol mutant allele";
//	}
	
	private String getGroup(Map<String, String> map){
		
		String geneSymbol=this.getGeneSymbol(map);
		
		if(map.containsKey("allele_symbol_superscript")){
			return geneSymbol+map.get("allele_symbol_superscript");
		}
		
		if(!map.get("pipeline").equals("")){
			return map.get("allele_symbol_superscript")+" "+ geneSymbol+ "mutant allele";
		}
		
		return geneSymbol+" mutant allele";
	}
//
//	sub _get_grouptype {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    my $design_type     = $feature_set->{design_type} || "";
//	    my $design_subtype  = $feature_set->{design_subtype} || "";
//	    
//	    if( $design_type =~ /Knock Out/ && $design_subtype ){
//	        return "Conditional($design_subtype)";
//	    }
//	    elsif( $design_type =~ /Ins/ ){
//	        return $design_subtype ? "Insertion($design_subtype)" : "Insertion";
//	    }
//	    elsif( $design_type =~ /Deletion/ ){
//	        return $design_subtype ? "Deletion($design_subtype)" : "Deletion";
//	    }
//	    else { return 'Conditional(Frameshift)'; }
//	}
	
	private String getGroupType(Map<String, String> map){
		String groupType="Conditional(Frameshift)";
		String designType=map.get("design_type");
		String designSubType=map.get("design_subtype");
//	    if( $design_type =~ /Knock Out/ && $design_subtype ){
//        return "Conditional($design_subtype)";
//    }
		if(designType.contains("Knock Out") && !designSubType.equals("")){
			return "Conditional ("+designSubType+")";
		}
//    elsif( $design_type =~ /Ins/ ){
//        return $design_subtype ? "Insertion($design_subtype)" : "Insertion";
//    }
		else if(designType.contains("Ins")){
			if(designSubType==null || designSubType.equals(""))return "Insertion";
			return "Insertion ("+designSubType+")";
		}
//    elsif( $design_type =~ /Deletion/ ){
//        return $design_subtype ? "Deletion($design_subtype)" : "Deletion";
//    }
		else if(designType.contains("Deletion")){
			if(designSubType==null || designSubType.equals(""))return "Deletion";
			return "Deletion("+designSubType+")";
		}else{
			return groupType;
		}
//    else { return 'Conditional(Frameshift)'; }
		
	}
//
//	sub _get_grouplink {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    if ( $feature_set->{pipeline} and $feature_set->{pipeline} eq "mirKO" ){
//	        return "";
//	    }
//	    
//	    my $mgi_accession_id = $feature_set->{mgi_accession_id}
//	        or return "";
//	    
//	    return "http://www.knockoutmouse.org/genedetails/$mgi_accession_id";
//	}
	
	private Map <URL,String> getGroupLink(Map<String, String> map, String groupId){
		String pipeline=map.get("pipeline");
		Map<URL, String>linkMap=new HashMap<URL, String>();
		if(pipeline.equals("mirKO")){
			return Collections.emptyMap();
		}
		String mgiAccession=map.get("mgi_accession_id");
		linkMap.put(IkmcLinks.getMgiDetailsLink(mgiAccession), groupId + " Gene Details");
		String alleleId=map.get("allele_id");
		linkMap.put(IkmcLinks.getCassetteImageLink(alleleId),"Cassette Image" );
		linkMap.put(IkmcLinks.getGbFileforEsCellClone(alleleId), "Genbank File");
		return linkMap;
	}
//
//	sub _get_groupnote {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    return "No status found" unless ( $feature_set->{ikmc_project_id} );
//	    
//	    my $dcc_biomart_result = $self->_get_dcc_biomart_result( $feature_set );
//	    
//	    if ( $dcc_biomart_result ) {
//	        return $dcc_biomart_result->{status};
//	    }
//	    
//	    return "No status found";
//	}
	
	private String getGroupNote(Map<String, String> map){
		String ikmcId=map.get("ikmc_project_id");
		//get cached if cached otherwise make a request to biomart
		return this.getGeneInfo(map, ikmcId, "status");
		
	}
//
//	sub _get_gene_symbol {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    unless ( $feature_set->{ikmc_project_id} ){
//	        return $self->_get_gene_symbol_from_mgi( $feature_set->{mgi_accession_id} );
//	    }
//	    
//	    my $dcc_biomart_result = $self->_get_dcc_biomart_result( $feature_set );
//	    
//	    if( $dcc_biomart_result ) {
//	        return $dcc_biomart_result->{marker_symbol};
//	    }
//	    
//	    return;
//	}
	
	private String getGeneSymbol(Map<String, String> map){
		//viveks code looks for the ikmc_project_id if it's present don't search mgi using the solr search
		//if not present do a solr request using the _get_gene_symbol_from_mgi( $feature_set->{mgi_accession_id} ); method
				
		if((map.get("ikmc_project_id")==null) || map.get("ikmc_project_id").equals("")){
			System.out.println("ikmc_project_id is null or empty so should do a solr request here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			return this.getGeneSymbolFromMgi(map.get("mgi_accession_id"));
		}
		String ikmcId=map.get("ikmc_project_id");
		return getGeneInfo(map, ikmcId, "marker_symbol");
		
	}

	private String getGeneInfo(Map<String, String> map, String ikmcId, String mapKey) {
		if(!ikmcGeneMap.containsKey(ikmcId)){
		Map<String,String>geneMap= this.biomart.getGeneNamesMapFromIkmcProjId(ikmcId);
		System.out.println("requesting biomart geneNames with id="+ikmcId);
		ikmcGeneMap.put(ikmcId, geneMap);
		return geneMap.get("marker_symbol");
		}else{
			return ikmcGeneMap.get(ikmcId).get(mapKey);
		}
	}
	

//
//	sub _get_gene_symbol_from_mgi {
//	    my ( $self, $mgi_accession_id ) = @_;
//	    
//	    unless ( exists $self->{gene_symbol}{$mgi_accession_id} ) {
//	        my $ua = LWP::UserAgent->new();
//	        $ua->proxy('http', 'http://wwwcache.sanger.ac.uk:3128');
//	        
//	        my $url = "http://www.sanger.ac.uk/mouseportal/solr/select?q=$mgi_accession_id&wt=json";
//	        my $response = $ua->get( $url );
//	        die "Couldn't get gene symbol from $url" unless defined $response;
//	        
//	        if ( defined $response->content ){
//	            my $results = from_json( $response->content );
//	            $self->{gene_symbol}{$mgi_accession_id} = $results->{response}{docs}[0]{marker_symbol};
//	        }
//	        else {
//	            $self->{gene_symbol}{$mgi_accession_id} = '';
//	        }
//	    }
//	    
//	    return $self->{gene_symbol}{$mgi_accession_id};
//	}
	
	private String getGeneSymbolFromMgi(String mgiId){
		//need to connect to solr here
		return "TODO need to connect to solr to geneSymbol??";
		
	}
//
//	sub _get_dcc_biomart_result {
//	    my ( $self, $feature_set ) = @_;
//	    
//	    my $ikmc_project_id = $feature_set->{ikmc_project_id};
//	    
//	    unless ( exists $self->{dcc_biomart_result}{$ikmc_project_id} ) {        
//	        $self->{dcc_biomart_result}{$ikmc_project_id} = $self->transport('dcc')->query(
//	            "ikmc_project_id" => $ikmc_project_id
//	        );
//	    }
//	    my $biomart_results = $self->{dcc_biomart_result}{$ikmc_project_id};
//	    
//	    if ( @{$biomart_results} ) {
//	        return $biomart_results->[0];
//	    }
//	    
//	    return;
//	}
	private Map<String,String> getDccBiomartResult(Map<String, String> map){
		String ikmcProjId=map.get("ikmc_project_id");
		Map<String,String>geneMap=biomart.getGeneNamesMapFromIkmcProjId(ikmcProjId);
		return geneMap;
	}
//
//	sub _get_position {
//	    my ( $self, $feature_set, $feature_name ) = @_;
//	    
//	    if( $feature_set->{strand} eq '-' ){
//	        return (
//	            'start' => $feature_set->{$feature_name."_end"},
//	            'end'   => $feature_set->{$feature_name."_start"}
//	        );
//	    }
//	    return (
//	        'start' => $feature_set->{$feature_name."_start"},
//	        'end'   => $feature_set->{$feature_name."_end"}
//	    );
//	}
	
	private Map<String, Integer> getPosition(Map<String, String> map, String type){
		Map<String, Integer> startStop=new HashMap<String, Integer>();
		//System.out.println(type+"_start");
		if(map.get("strand").equals("-")){
			//swap start for end
			startStop.put("end", Integer.valueOf(map.get(type+"_start")));
			startStop.put("start", Integer.valueOf(map.get(type+"_end")));
			
		}else{
			startStop.put("start", Integer.valueOf(map.get(type+"_start")));
			startStop.put("end", Integer.valueOf(map.get(type+"_end")));
			
		}
		
		return startStop;
	}

	
 private List<DasFeature> testFeatures(Map<String, String> map){
	 LinkedList<DasFeature> testFeatureList=new LinkedList<DasFeature>();
	 
		DasMethod testMethod=null;
		try {
			testMethod = new DasMethod("Test Method",	"TestMethod label", "ECO:0000043");
		} catch (DataSourceException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		List<String> testNotes=new ArrayList<String>();
			testNotes.add("test note");
			Map<URL,String> testLinks=new HashMap<URL, String>();
			URL testUrl=null;
			try {
				testUrl = new URL("http://www.google.com");
			} catch (MalformedURLException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			testLinks.put(testUrl, "test link Text");

			DasFeatureOrientation dOrientation = DasFeatureOrientation.ORIENTATION_SENSE_STRAND;
			if(map.get("strand").equals("-")){
				dOrientation=DasFeatureOrientation.ORIENTATION_ANTISENSE_STRAND;
			}
			
			
		
		List<String>parents=new ArrayList<String>();
		parents.add("parentFeatureId");
		List<String>parts=new ArrayList<String>();
		
		
		DasType dType = new DasType("Homology Arm start", "", "",
				"Homology Arm start");;
		DasFeature homArmStart = null;
		try {
			
			homArmStart = new DasFeature("HomArmStartId", "Homology Arm start", dType, method, Integer.valueOf(map.get("homology_arm_start"))-50,
					Integer.valueOf(map.get("homology_arm_start")), /*score*/null, dOrientation,
					DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ parents,
					/*parts*/null);
		} catch (DataSourceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		testFeatureList.add(homArmStart);
		parts.add("HomArmStartId");
		
		//cassette set up
		DasFeature cassette = null;
		
		DasType cassetteType=new DasType("Cassette ("+map.get("cassette")+")","","", "");
		String cassetteId="CassetteId"+" - Cassette insertion";
		try {
			
			cassette = new DasFeature(cassetteId, "Cassette ("+map.get("cassette")+")", cassetteType, method, this.getPosition(map, "cassette").get("start"),
					this.getPosition(map, "cassette").get("end"), /*score*/null, dOrientation,
					DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/parents,
					/*parts*/null);
		} catch (DataSourceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		testFeatureList.add(cassette);
		parts.add(cassetteId);
		
		
		//hom arm end
		DasType homEndType = new DasType("Homology Arm end", "", "",
				"Homology Arm end");;
		DasFeature homEnd = null;
		String homEndId="homEndId"+"- Homology Arm end";
		try {
			homEnd = new DasFeature(homEndId, "Homology Arm end", homEndType, method, Integer.valueOf(map.get("homology_arm_end")),
					Integer.valueOf(map.get("homology_arm_end"))+50, /*score*/null, dOrientation,
					DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ parents,
					/*parts*/null);
		} catch (DataSourceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		testFeatureList.add(homEnd);
		parts.add(homEndId);
		
		//FRT
		DasType frtType = new DasType("FRT", "FRT", "",
				"FRT");
		DasFeature frt = null;
		try {
			
			frt = new DasFeature("frtId", "FRT", frtType, testMethod, this.getPosition(map, "homology_arm").get("start")+2800,
					this.getPosition(map, "homology_arm").get("start")+2900, /*score*/null, dOrientation,
					DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ parents,
					/*parts*/null);
		} catch (DataSourceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		testFeatureList.add(frt);
		parts.add("frtId");
		//lacZ
				DasType laczType = new DasType("lacZ", "lacZ", "",
						"lacZ");
				DasFeature lacz = null;
				try {
					
					lacz = new DasFeature("lacZId", "lacZ", laczType, testMethod, this.getPosition(map, "homology_arm").get("start")+3000,
							this.getPosition(map, "homology_arm").get("start")+4000, /*score*/null, dOrientation,
							DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ parents,
							/*parts*/null);
				} catch (DataSourceException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				testFeatureList.add(lacz);
				parts.add("lacZId");
				
				//bact /neo?
				DasType bactType = new DasType("bact", "bact", "",
						"bact");
				DasFeature bact = null;
				String bactId="bactId";
				try {
					
					bact = new DasFeature(bactId, bactId, bactType, testMethod, this.getPosition(map, "homology_arm").get("start")+4500,
							this.getPosition(map, "homology_arm").get("start")+5000, /*score*/null, dOrientation,
							DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ parents,
							/*parts*/null);
				} catch (DataSourceException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				testFeatureList.add(bact);
				parts.add(bactId);
				
				//create the parent feature
				DasType testTypeParent = new DasType("Deletion(Frameshift)", "Deletion(Frameshift)", "",
						"Deletion(Frameshift)");
				
				DasFeature testFeatureParent = null;
				try {
					
					testFeatureParent = new DasFeature("parentFeatureId", "TestFeature label Tm1a Etc", testTypeParent, testMethod, this.getPosition(map, "homology_arm").get("start")-100,
							this.getPosition(map, "homology_arm").get("end"), /*score*/null, dOrientation,
							DasPhase.PHASE_NOT_APPLICABLE, /*notes Collection<String>*/testNotes,/*links Map<URL,String>*/ testLinks,/*targets*/ null,/*parents*/ null,
							/*parts*/parts);
				} catch (DataSourceException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
		testFeatureList.push(testFeatureParent);
		return testFeatureList;
 }
}
