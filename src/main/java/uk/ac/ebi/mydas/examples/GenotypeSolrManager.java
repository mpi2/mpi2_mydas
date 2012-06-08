package uk.ac.ebi.mydas.examples;

import java.io.UnsupportedEncodingException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

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

public class GenotypeSolrManager {
	/**
	 * List of the types used in this data source
	 */
	private ArrayList<DasType> types;
	/**
	 * Types to be used in the data source: chromosome, gene, transcript and
	 * exon
	 */
	private DasType geneType;
	/**
	 * As this data source just have one method, it can be defined as a
	 * parameter to facilitate its use
	 */
	private DasMethod method;

	private Connection connection;

	private String database = "genotype1";
	private Collection<DasEntryPoint> entryPoints = null;
	private String solrRootUrl;
	private String solrWithBasicParams;

	public GenotypeSolrManager(String solrUrl) throws DataSourceException {
		this.solrRootUrl = solrUrl;
		this.solrWithBasicParams=solrRootUrl+"/select/?version=2.2&start=0&rows=10000&indent=on&";
		method = new DasMethod("not_recorded", "not_recorded", "ECO:0000037");

		// Initialize types
		types = new ArrayList<DasType>();

		String typeId = "AA";
		System.out.println("typeid=" + typeId);
		types.add(new DasType(typeId, "", "SO:0000694", ""));

	}

	public void close() {
		try {
			connection.close();
		} catch (Exception e) { /* ignore close errors */
		}

	}

	public DasAnnotatedSegment getSubmodelBySegmentId(String segmentId,
			int start, int stop) throws DataSourceException,
			BadReferenceObjectException {
		return getSubmodelBySegmentId(segmentId, start, stop, -1);
	}

	public DasAnnotatedSegment getSubmodelBySegmentId(String segmentId,
			int start, int stop, int maxbins) throws DataSourceException,
			BadReferenceObjectException {

		// http://localhost:8983/solr/select/?q=position_i%3A%5B592942+TO+595687%5D&version=2.2&start=0&rows=10000&indent=on
		// %3A%5B592942+TO+595687%5D
		String url = solrWithBasicParams + "q=position_i:[" + start + " TO " + stop
				+ "] AND chr_s:"+segmentId;// need to add segment id
		// String url = solrRootUrl + "q=position_i:[" + start + " TO " + stop
		// + "]";// need to add segment id
		// http://localhost:8983/solr/select/?version=2.2&start=0&rows=10000&indent=on&q=position_i:[592942%20TO%20595687]
		System.out.println(url);
		Collection<DasFeature> features = getGenotypeFeatures(url);
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

	private Collection<DasFeature> getGenotypeFeatures(String url) {

		Collection<DasFeature> features = new ArrayList();
		// send a request to the url
		RestTemplate template = new RestTemplate();
		//logger.debug("url to solr=" + url);
		// try {
		// url=UriUtils.encodeUri(url, "UTF-8");
		// } catch (UnsupportedEncodingException e2) {
		// // TODO Auto-generated catch block
		// e2.printStackTrace();
		// }
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
		System.out.println(number);
		List<Doc> docs = result.getDoc();
		for (Doc doc : docs) {
			String typeId = "";
			String id = "";
			int position = 0;

			List<Object> fields = doc.getStrOrIntOrFloat();
			for (Object field : fields) {
				if (field instanceof Str) {
					Str str = (Str) field;
					String name = str.getName();
					if (name.equals("chr_s")) {

					}
					if (name.equals("genotype_s")) {
						typeId = str.getValue();
					}
					if (name.equals("id")) {
						id = str.getValue();
					}

				} else if (field instanceof Int) {
					Int positionObject = (Int) field;

					position = positionObject.getValue();

				}
			}

			// rs.getInt("chromosome");
			DasType type = new DasType(typeId, "", "SO:0000694", "");
			DasMethod method = null;
			try {
				method = new DasMethod("23AndMe", "microarray", "");
			} catch (DataSourceException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			DasFeature feature = null;
			try {
				feature = new DasFeature(id, id, type, method, position,
						position, new Double(1),
						DasFeatureOrientation.ORIENTATION_NOT_APPLICABLE,
						DasPhase.PHASE_NOT_APPLICABLE, null, null, null, null,
						null);
			} catch (DataSourceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			features.add(feature);

		}

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

		
		return 10;
	}

	public String getDatabase() {
		return database;
	}

}
