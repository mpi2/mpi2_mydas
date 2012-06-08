package uk.ac.ebi.mydas.examples;

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

public class MouseSolrManager {
	private final BiomartDao biomart;
	/**
	 * List of the types used in this data source
	 */
	private ArrayList<DasType> types;
	/**
	 * Types to be used in the data source: chromosome, gene, transcript and
	 * exon
	 */
	private final DasType geneType = new DasType("Gene", "", "SO:0000704",
			"Gene");;
	/**
	 * As this data source just have one method, it can be defined as a
	 * parameter to facilitate its use
	 */
	private final DasMethod method = new DasMethod("MGI",
			"combinatorial evidence", "ECO:0000043");

	

	private String database = "genotype1";
	private Collection<DasEntryPoint> entryPoints = null;
	private String solrRootUrl;
	private String solrWithBasicParams;
	private Map<String, Map<String, String>> ikmcGeneMap;

	public MouseSolrManager(String solrUrl) throws DataSourceException {
		this.solrRootUrl = solrUrl;
		this.solrWithBasicParams = solrRootUrl
				+ "/select/?version=2.2&start=0&rows=100000&indent=on&";
		biomart = new BiomartDao("http://www.i-dcc.org/biomart/martservice");
		// Initialize types
		types = new ArrayList<DasType>();
		types.add(geneType);
		ikmcGeneMap=new HashMap<String,Map<String,String>>();

	}


	public DasAnnotatedSegment getSubmodelBySegmentId(String segmentId,
			int start, int stop) throws DataSourceException,
			BadReferenceObjectException {
		return getSubmodelBySegmentId(segmentId, start, stop, -1);
	}

	public DasAnnotatedSegment getSubmodelBySegmentId(String segmentId,
			int start, int stop, int maxbins) throws DataSourceException,
			BadReferenceObjectException {

		Collection<DasFeature> features = getFeatures(start, stop, segmentId);
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

	private Collection<DasFeature> getFeatures(int start, int stop,
			String segmentId) {
		String url = solrWithBasicParams + "q=feature_start_i:[* TO "
				+ stop + "] AND feature_end_i:[" + start
				+ " TO *] AND chromosome_s:" + segmentId;

		List<Map<String, String>> listOfMaps = getFeaturesMaps(url);
		System.out.println("solr maps size=" + listOfMaps.size());
		LinkedList<DasFeature> features = new LinkedList<DasFeature>();
		//int i = 0;
		for (Map<String, String> map : listOfMaps) {
			// if no allel superscript assume no allele and so don't include
			// that group in the features response
			if (map.containsKey("allele_symbol_superscript")) {
			

				String dGroup = this.getGroup(map);
				//System.out.println("group="+dGroup);
				String dGroupType = this.getGroupType(map);
				// System.out.println("group type="+dGroupType);
				Map dGroupLink = this.getGroupLink(map, dGroup);
				String dGroupNote = this.getGroupNote(map);
				String dNote = dGroupNote;
				Collection<String> dNotes = new ArrayList<String>();
				dNotes.add(dNote);
				Map dLink = dGroupLink;
				DasFeatureOrientation dOrientation = DasFeatureOrientation.ORIENTATION_SENSE_STRAND;
				if(map.containsKey("strand")){
				if (map.get("strand").equals("-")) {
					dOrientation = DasFeatureOrientation.ORIENTATION_ANTISENSE_STRAND;
				}
				}
				// set up an array that holds just one parent and use for all
				// subsequent features in this set
				ArrayList<String> parents = new ArrayList<String>();
				parents.add(dGroup);
				// set up the parts array to store the name of all the sub
				// features of the 'group'
				ArrayList<String> parts = new ArrayList<String>();

				
				
				DasFeature featureG5 = null;
				String startString="Homology Arm start";
				String endString="Homology Arm end";
				if(map.get("strand").equals("-")){
					startString="Homology Arm end";
					endString="Homology Arm start";
				}
				String startId=dGroup +" - "+ startString;
				String endId=dGroup +" - "+  endString;
				
				DasType dType = new DasType(startString, "", "",
						startString);
				try {

					featureG5 = new DasFeature(startId, startString,
							dType, method, Integer.valueOf(map
									.get("feature_start")),
							Integer.valueOf(map.get("feature_start"))+50, /* score */
							null, dOrientation, DasPhase.PHASE_NOT_APPLICABLE, /*
																				 * notes
																				 * Collection
																				 * <
																				 * String
																				 * >
																				 */
							dNotes,/* links Map<URL,String> */dLink,/* targets */
							null,/* parents */parents,
							/* parts */null);
				} catch (DataSourceException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				features.add(featureG5);
				parts.add(startId);
				
				DasType homEndType = new DasType(endString, "", "",
						endString);
				
				DasFeature featureG3 = null;
				
								
				try {
					featureG3 = new DasFeature(endId,
							endString, homEndType, method,
							Integer.valueOf(map.get("feature_end"))-50,
							Integer.valueOf(map.get("feature_end")) , /* score */
							null, dOrientation, DasPhase.PHASE_NOT_APPLICABLE, /*
																				 * notes
																				 * Collection
																				 * <
																				 * String
																				 * >
																				 */
							dNotes,/* links Map<URL,String> */dLink,/* targets */
							null,/* parents */parents,
							/* parts */null);
				} catch (DataSourceException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				features.add(featureG3);
				parts.add(endId);

				
				DasFeature featureU5U3 = null;
				if (this.getPosition(map, "cassette").get("start") <= stop
						&& this.getPosition(map, "cassette").get("end") >= start) {
					DasType cassetteType = new DasType("Cassette ("
							+ map.get("cassette") + ")", "", "", "");
					String u5u3FeatureId = dGroup + " - Cassette insertion";
					try {

						featureU5U3 = new DasFeature(u5u3FeatureId,
								"Cassette (" + map.get("cassette") + ")",
								cassetteType, method, this.getPosition(map,
										"cassette").get("start"), this
										.getPosition(map, "cassette")
										.get("end"), /* score */null,
								dOrientation, DasPhase.PHASE_NOT_APPLICABLE, /*
																			 * notes
																			 * Collection
																			 * <
																			 * String
																			 * >
																			 */
								dNotes,/* links Map<URL,String> */dLink,/* targets */
								null,/* parents */parents,
								/* parts */null);
					} catch (DataSourceException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					features.add(featureU5U3);
					parts.add(u5u3FeatureId);
				}
				
				if (map.containsKey("loxp_start")&& map.containsKey("loxp_end")) {
					// if(Integer.valueOf(map.get("loxp_start"))<= stop &&
					// Integer.valueOf(map.get("loxp_end")) >= start){
					DasFeature featureD5U3 = null;
					// System.out.println("creating loxp feature");
					DasType loxpType = new DasType("LoxP", "LoxP", "", "LoxP");
					String d5u3Id = dGroup + " - LoxP insertion";
					try {
						featureD5U3 = new DasFeature(d5u3Id, "LoxP", loxpType,
								method, this.getPosition(map, "loxp").get(
										"start"), this.getPosition(map, "loxp")
										.get("end"), /* score */null,
								dOrientation, DasPhase.PHASE_NOT_APPLICABLE, /*
																			 * notes
																			 * Collection
																			 * <
																			 * String
																			 * >
																			 */
								dNotes,/* links Map<URL,String> */dLink,/* targets */
								null,/* parents */parents,
								/* parts */null);
					} catch (DataSourceException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					features.add(featureD5U3);
					parts.add(d5u3Id);
					// }
				}
				

				
				// set up the group/parent id here which uses the default for
				// groups for this set and has to be set up last so we know what
				// the children are for the constructor:
				DasType defaultGroupType = new DasType(dGroupType, dGroupType,
						"", dGroupType);
				DasFeature group = null;

				try {
					group = new DasFeature(dGroup, null, defaultGroupType,
							method, Integer.parseInt(map.get("feature_start")) - 50, 
									Integer.parseInt(map.get("feature_end")), /* score */
							null, dOrientation, DasPhase.PHASE_NOT_APPLICABLE, /*
																				 * notes
																				 * Collection
																				 * <
																				 * String
																				 * >
																				 */
							dNotes,/* links Map<URL,String> */dLink,/* targets */
							null,/* parents */null,
							/* parts */parts);
				} catch (DataSourceException e) {
					e.printStackTrace();
				}
				features.push(group);
			}
		}
		
		return features;

	}

	private List<Map<String, String>> getFeaturesMaps(String url) {

		List<Map<String, String>> listOfMaps = new ArrayList<Map<String, String>>();
		//System.out.println("maps size=" + listOfMaps.size());
		// send a request to the url
		RestTemplate template = new RestTemplate();
		System.out.println("url to solr after encoding=" + url);
		Response response = null;
		try {
			response = template.getForObject(url, Response.class);
		} catch (RestClientException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		// parse the solr xml so we can generate DAS features
		Result result = response.getResult();

		int number = result.getNumFound();
		//System.out.println(number);
		List<Doc> docs = result.getDoc();
		for (Doc doc : docs) {
			Map<String, String> map = new HashMap<String, String>();

			List<Object> fields = doc.getStrOrIntOrFloat();
			for (Object field : fields) {
				

				if (field instanceof Str) {
					Str str = (Str) field;
					String name = str.getName();
					String id = str.getValue();
					// map.put("pipeline", fields[0]);
					// map.put("mgi_accession_id",fields[1]);
					// map.put("ikmc_project_id",fields[2]);
					// map.put("design_type",fields[3]);
					// map.put("design_subtype", fields[4]);
					// map.put("cassette",fields[5]);
					// map.put("allele_symbol_superscript",fields[6]);
					// map.put("strand",fields[13]);
					// map.put("allele_id",fields[14]);
					if (name.equals("pipeline")) {
						map.put("pipeline_s", id);
					}
					if (name.equals("mgi_accession_id_s")) {

						map.put("mgi_accession_id", id);
					}
					if (name.equals("ikmc_project_id_s")) {
						map.put("ikmc_project_id", id);
					}
					if (name.equals("design_type_s")) {
						map.put("design_type", id);
					}
					if (name.equals("cassette_s")) {
						map.put("cassette", id);
					}
					if (name.equals("allele_symbol_superscript_s")) {
						map.put("allele_symbol_superscript", id);
					}
					if (name.equals("strand_s")) {
						map.put("strand", id);
					}
					if (name.equals("allele_id_s")) {
						map.put("allele_id", id);
					}

				} else if (field instanceof Int) {

					Int positionObject = (Int) field;
					String intName = positionObject.getName();
					// map.put("homology_arm_start", fields[7]);
					// map.put("homology_arm_end",fields[8]);
					// map.put("loxp_start",fields[9]);
					// map.put("loxp_end",fields[10]);
					// map.put("cassette_start",fields[11]);
					// map.put("cassette_end", fields[12]);
					int genericInteger = positionObject.getValue();
					if (intName.equals("feature_start_i")) {

						map.put("feature_start",
								String.valueOf(genericInteger));
					}
					if (intName.equals("feature_end_i")) {

						map.put("feature_end",
								String.valueOf(genericInteger));
					}
					if (intName.equals("loxp_start_i")) {

						map.put("loxp_start", String.valueOf(genericInteger));
					}
					if (intName.equals("loxp_end_i")) {

						map.put("loxp_end", String.valueOf(genericInteger));
					}
					if (intName.equals("cassette_start_i")) {

						map.put("cassette_start",
								String.valueOf(genericInteger));
					}
					if (intName.equals("cassette_end_i")) {

						map.put("cassette_end", String.valueOf(genericInteger));
					}

				}
				
			}
			listOfMaps.add(map);
		}

		return listOfMaps;
	}

	public ArrayList<DasType> getTypes() {
		//System.out.println("Getting types from genotype manager");
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

	// viveks helper methods implement each :
	// sub _get_group {
	// my ( $self, $feature_set ) = @_;
	//
	// my $gene_symbol = $self->_get_gene_symbol( $feature_set ) || "";
	//
	// if ( $feature_set->{allele_symbol_superscript} ) {
	// return "$gene_symbol $feature_set->{allele_symbol_superscript}";
	// }
	//
	// if ( $feature_set->{pipeline} ) {
	// return "$feature_set->{pipeline} $gene_symbol mutant allele";
	// }
	//
	// return "$gene_symbol mutant allele";
	// }

	private String getGroup(Map<String, String> map) {

		String geneSymbol = this.getGeneSymbol(map);

		if (map.containsKey("allele_symbol_superscript")) {
			return geneSymbol +"<"+ map.get("allele_symbol_superscript")+">";
		}

		if (!map.get("pipeline").equals("")) {
			return map.get("allele_symbol_superscript") + " " + geneSymbol
					+ "mutant allele";
		}

		return geneSymbol + " mutant allele";
	}

	//
	// sub _get_grouptype {
	// my ( $self, $feature_set ) = @_;
	//
	// my $design_type = $feature_set->{design_type} || "";
	// my $design_subtype = $feature_set->{design_subtype} || "";
	//
	// if( $design_type =~ /Knock Out/ && $design_subtype ){
	// return "Conditional($design_subtype)";
	// }
	// elsif( $design_type =~ /Ins/ ){
	// return $design_subtype ? "Insertion($design_subtype)" : "Insertion";
	// }
	// elsif( $design_type =~ /Deletion/ ){
	// return $design_subtype ? "Deletion($design_subtype)" : "Deletion";
	// }
	// else { return 'Conditional(Frameshift)'; }
	// }

	private String getGroupType(Map<String, String> map) {
		String groupType = "Conditional(Frameshift)";
		String designType = map.get("design_type");
		String designSubType = map.get("design_subtype");
		// if( $design_type =~ /Knock Out/ && $design_subtype ){
		// return "Conditional($design_subtype)";
		// }
		if (designType != null && designSubType != null) {
			if (designType.contains("Knock Out") && !designSubType.equals("")) {
				return "Conditional (" + designSubType + ")";
			}

			// elsif( $design_type =~ /Ins/ ){
			// return $design_subtype ? "Insertion($design_subtype)" :
			// "Insertion";
			// }
			else if (designType.contains("Ins")) {
				if (designSubType == null || designSubType.equals(""))
					return "Insertion";
				return "Insertion (" + designSubType + ")";
			}
			// elsif( $design_type =~ /Deletion/ ){
			// return $design_subtype ? "Deletion($design_subtype)" :
			// "Deletion";
			// }
			else if (designType.contains("Deletion")) {
				if (designSubType == null || designSubType.equals(""))
					return "Deletion";
				return "Deletion(" + designSubType + ")";
			}
		}
		return groupType;

		// else { return 'Conditional(Frameshift)'; }

	}

	//
	// sub _get_grouplink {
	// my ( $self, $feature_set ) = @_;
	//
	// if ( $feature_set->{pipeline} and $feature_set->{pipeline} eq "mirKO" ){
	// return "";
	// }
	//
	// my $mgi_accession_id = $feature_set->{mgi_accession_id}
	// or return "";
	//
	// return "http://www.knockoutmouse.org/genedetails/$mgi_accession_id";
	// }

	private Map<URL, String> getGroupLink(Map<String, String> map,
			String groupId) {
		
		Map<URL, String> linkMap = new HashMap<URL, String>();
		if(map.containsKey("pipeline")){
		String pipeline = map.get("pipeline");
		if (pipeline.equals("mirKO")) {
			return Collections.emptyMap();
		}
		}
		String mgiAccession = map.get("mgi_accession_id");
		linkMap.put(IkmcLinks.getMgiDetailsLink(mgiAccession), groupId
				+ " Gene Details");
		String alleleId = map.get("allele_id");
		linkMap.put(IkmcLinks.getCassetteImageLink(alleleId), "Cassette Image");
		linkMap.put(IkmcLinks.getGbFileforEsCellClone(alleleId), "Genbank File");
		return linkMap;
	}

	//
	// sub _get_groupnote {
	// my ( $self, $feature_set ) = @_;
	//
	// return "No status found" unless ( $feature_set->{ikmc_project_id} );
	//
	// my $dcc_biomart_result = $self->_get_dcc_biomart_result( $feature_set );
	//
	// if ( $dcc_biomart_result ) {
	// return $dcc_biomart_result->{status};
	// }
	//
	// return "No status found";
	// }

	private String getGroupNote(Map<String, String> map) {
		String ikmcId = map.get("ikmc_project_id");
		// get cached if cached otherwise make a request to biomart
		return this.getGeneInfo(map, ikmcId, "status");

	}

	//
	// sub _get_gene_symbol {
	// my ( $self, $feature_set ) = @_;
	//
	// unless ( $feature_set->{ikmc_project_id} ){
	// return $self->_get_gene_symbol_from_mgi( $feature_set->{mgi_accession_id}
	// );
	// }
	//
	// my $dcc_biomart_result = $self->_get_dcc_biomart_result( $feature_set );
	//
	// if( $dcc_biomart_result ) {
	// return $dcc_biomart_result->{marker_symbol};
	// }
	//
	// return;
	// }

	private String getGeneSymbol(Map<String, String> map) {
		// viveks code looks for the ikmc_project_id if it's present don't
		// search mgi using the solr search
		// if not present do a solr request using the _get_gene_symbol_from_mgi(
		// $feature_set->{mgi_accession_id} ); method

		if ((map.get("ikmc_project_id") == null)
				|| map.get("ikmc_project_id").equals("")) {
			System.out
					.println("ikmc_project_id is null or empty so should do a solr request here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			return this.getGeneSymbolFromMgi(map.get("mgi_accession_id"));
		}
		String ikmcId = map.get("ikmc_project_id");
		return getGeneInfo(map, ikmcId, "marker_symbol");

	}

	private String getGeneInfo(Map<String, String> map, String ikmcId,
			String mapKey) {
		if (!ikmcGeneMap.containsKey(ikmcId)) {
			Map<String, String> geneMap = this.biomart
					.getGeneNamesMapFromIkmcProjId(ikmcId);
			System.out
					.println("requesting biomart geneNames with id=" + ikmcId);
			ikmcGeneMap.put(ikmcId, geneMap);
			return geneMap.get("marker_symbol");
		} else {
			return ikmcGeneMap.get(ikmcId).get(mapKey);
		}
	}

	//
	// sub _get_gene_symbol_from_mgi {
	// my ( $self, $mgi_accession_id ) = @_;
	//
	// unless ( exists $self->{gene_symbol}{$mgi_accession_id} ) {
	// my $ua = LWP::UserAgent->new();
	// $ua->proxy('http', 'http://wwwcache.sanger.ac.uk:3128');
	//
	// my $url =
	// "http://www.sanger.ac.uk/mouseportal/solr/select?q=$mgi_accession_id&wt=json";
	// my $response = $ua->get( $url );
	// die "Couldn't get gene symbol from $url" unless defined $response;
	//
	// if ( defined $response->content ){
	// my $results = from_json( $response->content );
	// $self->{gene_symbol}{$mgi_accession_id} =
	// $results->{response}{docs}[0]{marker_symbol};
	// }
	// else {
	// $self->{gene_symbol}{$mgi_accession_id} = '';
	// }
	// }
	//
	// return $self->{gene_symbol}{$mgi_accession_id};
	// }

	private String getGeneSymbolFromMgi(String mgiId) {
		// need to connect to solr here
		return "TODO need to connect to solr to geneSymbol??";

	}

	//
	// sub _get_dcc_biomart_result {
	// my ( $self, $feature_set ) = @_;
	//
	// my $ikmc_project_id = $feature_set->{ikmc_project_id};
	//
	// unless ( exists $self->{dcc_biomart_result}{$ikmc_project_id} ) {
	// $self->{dcc_biomart_result}{$ikmc_project_id} =
	// $self->transport('dcc')->query(
	// "ikmc_project_id" => $ikmc_project_id
	// );
	// }
	// my $biomart_results = $self->{dcc_biomart_result}{$ikmc_project_id};
	//
	// if ( @{$biomart_results} ) {
	// return $biomart_results->[0];
	// }
	//
	// return;
	// }
	private Map<String, String> getDccBiomartResult(Map<String, String> map) {
		String ikmcProjId = map.get("ikmc_project_id");
		Map<String, String> geneMap = biomart
				.getGeneNamesMapFromIkmcProjId(ikmcProjId);
		return geneMap;
	}

	//
	// sub _get_position {
	// my ( $self, $feature_set, $feature_name ) = @_;
	//
	// if( $feature_set->{strand} eq '-' ){
	// return (
	// 'start' => $feature_set->{$feature_name."_end"},
	// 'end' => $feature_set->{$feature_name."_start"}
	// );
	// }
	// return (
	// 'start' => $feature_set->{$feature_name."_start"},
	// 'end' => $feature_set->{$feature_name."_end"}
	// );
	// }

	private Map<String, Integer> getPosition(Map<String, String> map,
			String type) {
		Map<String, Integer> startStop = new HashMap<String, Integer>();
		// System.out.println(type+"_start");
		if (map.get("strand").equals("-")) {
			// swap start for end
			startStop.put("end", Integer.valueOf(map.get(type + "_start")));
			startStop.put("start", Integer.valueOf(map.get(type + "_end")));

		} else {
			startStop.put("start", Integer.valueOf(map.get(type + "_start")));
			startStop.put("end", Integer.valueOf(map.get(type + "_end")));

		}

		return startStop;
	}

}
