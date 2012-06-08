package uk.ac.ebi.mydas.examples;

import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;

import javax.servlet.ServletContext;

import org.apache.log4j.Logger;

import uk.ac.ebi.mydas.configuration.DataSourceConfiguration;
import uk.ac.ebi.mydas.configuration.PropertyType;
import uk.ac.ebi.mydas.controller.CacheManager;
import uk.ac.ebi.mydas.datasource.RangeHandlingAnnotationDataSource;
import uk.ac.ebi.mydas.exceptions.BadReferenceObjectException;
import uk.ac.ebi.mydas.exceptions.CoordinateErrorException;
import uk.ac.ebi.mydas.exceptions.DataSourceException;
import uk.ac.ebi.mydas.exceptions.UnimplementedFeatureException;
import uk.ac.ebi.mydas.extendedmodel.DasUnknownFeatureSegment;
import uk.ac.ebi.mydas.model.*;

/**
 * Ikmc Alleles data source based on old sanger ikmc_products data source
 * [ikmc_products]
description    = IKMC products
adaptor        = ikmc_products
state          = on
doc_href       = http://www.i-dcc.org/targ_rep/
coordinates    = NCBIM_37,Chromosome,Mus musculus -> 6:30040571,30229278
capabilities   = features -> 1.0
transport      = martservice
martservice    = http://www.i-dcc.org/biomart/martservice/
dataset        = idcc_targ_rep
attributes     = pipeline,mgi_accession_id,ikmc_project_id,design_type,design_subtype,cassette,allele_symbol_superscript,homology_arm_start,homology_arm_end,loxp_start,loxp_end,cassette_start,cassette_end,strand
timeout        = 10
dcc.martservice = http://www.i-dcc.org/biomart/martservice/
dcc.dataset     = dcc
dcc.attributes  = status,marker_symbol
stylesheet       = <<EOT
<DASSTYLE>
  <STYLESHEET version="1.0">
  <CATEGORY id="default">
    <TYPE id="default">
      <GLYPH>
          <BOX>
            <FGCOLOR>black</FGCOLOR>
            <FONT>sanserif</FONT>
            <BGCOLOR>black</BGCOLOR>
          </BOX>
        </GLYPH>
      </TYPE>
    </CATEGORY>
  </STYLESHEET>
</DASSTYLE>
EOT

 * 
 * LPGL
 *@author JWarren
 *
 */
public class IkmcAllelesDataSource implements RangeHandlingAnnotationDataSource{
	CacheManager cacheManager = null;
	ServletContext svCon;
	Map<String, PropertyType> globalParameters;
	DataSourceConfiguration config;
	IkmcAllelesManager ikmcAllelesManager;
	static Logger logger = Logger.getLogger(IkmcAllelesDataSource.class);

	public void init(ServletContext servletContext,
			Map<String, PropertyType> globalParameters,
			DataSourceConfiguration dataSourceConfig)
			throws DataSourceException {
		this.svCon = servletContext;
		this.globalParameters = globalParameters;
		this.config = dataSourceConfig;
		   String biomartServiceUrl="";
		 
		if (config.getDataSourceProperties().containsKey("biomartServiceUrl")){
			   biomartServiceUrl=config.getDataSourceProperties().get("biomartServiceUrl").getValue();
		   }
		
		if(biomartServiceUrl.equals("") || biomartServiceUrl==null ){
			throw new DataSourceException("a database url must be set such in the configuration for example: ");
		}
		try {
			logger.debug("=========connection params="+biomartServiceUrl);
			ikmcAllelesManager=new IkmcAllelesManager(biomartServiceUrl);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void destroy() {
		
	}

	public DasAnnotatedSegment getFeatures(String segmentId, int start, int stop, Integer maxbins) 
			throws BadReferenceObjectException, CoordinateErrorException, DataSourceException {
		if (maxbins==null)
			maxbins=-1;
		return ikmcAllelesManager.getFeatures(segmentId, start, stop,maxbins);
	}
	public DasAnnotatedSegment getFeatures(String segmentId, Integer maxbins) 
			throws BadReferenceObjectException, DataSourceException {
		if (maxbins==null)
			maxbins=-1;
		return ikmcAllelesManager.getFeatures(segmentId, -1, -1,maxbins);
	}

	public Collection<DasType> getTypes() throws DataSourceException {
		return ikmcAllelesManager.getTypes();
	}
	public URL getLinkURL(String field, String id)
			throws UnimplementedFeatureException, DataSourceException {
		throw new UnimplementedFeatureException("No implemented");
	}


    public String getEntryPointVersion() throws UnimplementedFeatureException, DataSourceException {
        return ikmcAllelesManager.getDatabase();
    }

   

    public Integer getTotalCountForType(DasType type)
			throws DataSourceException {
		return ikmcAllelesManager.getTotalCountForType(type.getId());
	}

	public void registerCacheManager(CacheManager cacheManager) {
		this.cacheManager = cacheManager;
	}

	@Override
	public DasAnnotatedSegment getFeatures(String segmentId, Integer maxbins,
			Range rows) throws BadReferenceObjectException,
			DataSourceException, UnimplementedFeatureException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Collection<DasAnnotatedSegment> getFeatures(
			Collection<String> featureIdCollection, Integer maxbins, Range rows)
			throws UnimplementedFeatureException, DataSourceException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DasAnnotatedSegment getFeatures(String segmentId, int start,
			int stop, Integer maxbins, Range rows)
			throws BadReferenceObjectException, CoordinateErrorException,
			DataSourceException, UnimplementedFeatureException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Collection<DasAnnotatedSegment> getFeatures(
			Collection<String> featureIdCollection, Integer maxbins)
			throws UnimplementedFeatureException, DataSourceException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Collection<DasEntryPoint> getEntryPoints(Integer start, Integer stop)
			throws UnimplementedFeatureException, DataSourceException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getTotalEntryPoints() throws UnimplementedFeatureException,
			DataSourceException {
		// TODO Auto-generated method stub
		return 0;
	}
}