//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vhudson-jaxb-ri-2.2.1.1-4 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2011.11.24 at 01:00:07 PM GMT 
//


package uk.ac.ebi.mydas.examples.solr;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for anonymous complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType>
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;attribute name="waitFlush" type="{http://www.w3.org/2001/XMLSchema}boolean" />
 *       &lt;attribute name="waitSearcher" type="{http://www.w3.org/2001/XMLSchema}boolean" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "")
@XmlRootElement(name = "commit")
public class Commit {

    @XmlAttribute(name = "waitFlush")
    protected Boolean waitFlush;
    @XmlAttribute(name = "waitSearcher")
    protected Boolean waitSearcher;

    /**
     * Gets the value of the waitFlush property.
     * 
     * @return
     *     possible object is
     *     {@link Boolean }
     *     
     */
    public Boolean isWaitFlush() {
        return waitFlush;
    }

    /**
     * Sets the value of the waitFlush property.
     * 
     * @param value
     *     allowed object is
     *     {@link Boolean }
     *     
     */
    public void setWaitFlush(Boolean value) {
        this.waitFlush = value;
    }

    /**
     * Gets the value of the waitSearcher property.
     * 
     * @return
     *     possible object is
     *     {@link Boolean }
     *     
     */
    public Boolean isWaitSearcher() {
        return waitSearcher;
    }

    /**
     * Sets the value of the waitSearcher property.
     * 
     * @param value
     *     allowed object is
     *     {@link Boolean }
     *     
     */
    public void setWaitSearcher(Boolean value) {
        this.waitSearcher = value;
    }

}
