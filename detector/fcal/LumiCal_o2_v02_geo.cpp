#include <DD4hep/DetFactoryHelper.h>
#include "DD4hep/DetType.h"
#include <XML/Layering.h>
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"

#include <string>

using dd4hep::Assembly;
using dd4hep::BitField64;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::Box;
using dd4hep::Cone;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::DetType;
using dd4hep::IntersectionSolid;
using dd4hep::Layering;
using dd4hep::Layer;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::PolyhedraRegular;
using dd4hep::Position;
using dd4hep::Readout;
using dd4hep::Ref_t;
using dd4hep::Rotation3D;
using dd4hep::RotationY;
using dd4hep::RotationZ;
using dd4hep::RotationZYX;
using dd4hep::SensitiveDetector;
using dd4hep::SubtractionSolid;
using dd4hep::Transform3D;
using dd4hep::Translation3D;
using dd4hep::Trapezoid;
using dd4hep::Tube;
using dd4hep::ConeSegment;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::UnionSolid;
using dd4hep::IntersectionSolid;
using dd4hep::Segmentation;
using dd4hep::EllipticalTube;

using dd4hep::rec::LayeredCalorimeterData;

static Ref_t create_detector(Detector& theDetector,
                                               xml_h element,
                                               SensitiveDetector /*sens*/) {

    //Materials
    Material air = theDetector.air();
    
    //Access to the XML File
    xml_det_t     xmlLumiCal    = element;
    const std::string detName   = xmlLumiCal.nameStr();
    
    DetElement sdet ( detName, xmlLumiCal.id() );
    
    // --- create an envelope volume and position it into the world ---------------------
    
    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector, element , sdet ) ;
    
    sdet.setTypeFlag( DetType::CALORIMETER |  DetType::ENDCAP  | DetType::ELECTROMAGNETIC |  DetType::FORWARD ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;
    
    //-----------------------------------------------------------------------------------

    //Parameters we have to know about
    dd4hep::xml::Component xmlParameter = xmlLumiCal.child(_Unicode(parameter));
    const double fullCrossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));

    dd4hep::xml::Dimension dimensions =  xmlLumiCal.dimensions();
    
    //LumiCal Dimensions
    const double lcalInnerZ = dimensions.inner_z();
    const double lcalOuterZ = dimensions.outer_z() ;
    const double rInnerStart = dimensions.rmin1();
    const double rOuterStart = dimensions.rmax1();
    const double rInnerEnd = dimensions.rmin2();
    const double rOuterEnd = dimensions.rmax2();

    //const double zHalf = dimensions.zhalf();
    const double lcalThickness = Layering(xmlLumiCal).totalThickness();
    //const double lcalThickness = 2*zHalf;
    const double lcalCentreZ = lcalOuterZ-lcalThickness*0.5;
 
    
    //double LumiCal_cell_size      = theDetector.constant<double>("LumiCal_cell_size");

    //** DD4hep/TGeo seems to need rad (as opposed to the manual)
    const double phi1 = 0 ;
    const double phi2 = 360.0*dd4hep::degree;
    const double thetaInn = atan(( rInnerEnd - rInnerStart ) / lcalThickness);
    const double thetaOut = atan(( rOuterEnd - rOuterStart ) / lcalThickness) ;


    // counter for the current layer to be placed
    int thisLayerId = 0;
   
    
    //Envelope to place the layers in
    //Moving from cylindrical to conical geometry
    //DD4hep::Geometry::Tube envelopeTube (lcalInnerR, lcalOuterR, lcalThickness*0.5 );
    ConeSegment envelopeCone (lcalThickness*0.5, rInnerStart, rOuterStart, rInnerEnd, rOuterEnd, phi1, phi2) ;
    Volume envelopeVol(detName+"_module",envelopeCone,air);
    envelopeVol.setVisAttributes(theDetector,xmlLumiCal.visStr());
 



    ////////////////////////////////////////////////////////////////////////////////
    // Create all the layers
    ////////////////////////////////////////////////////////////////////////////////

    //Loop over all the layer (repeat=NN) sections
    //This is the starting point to place all layers, we need this when we have more than one layer block
    double referencePosition = -lcalThickness*0.5;
    for(dd4hep::xml::Collection_t coll(xmlLumiCal,_U(layer)); coll; ++coll)  {
        dd4hep::xml::Component xmlLayer(coll); //we know this thing is a layer
        
        
        //This just calculates the total size of a single layer
        //Why no convenience function for this?
        double layerThickness = 0;
        for(dd4hep::xml::Collection_t l(xmlLayer,_U(slice)); l; ++l)
            layerThickness += xml_comp_t(l).thickness();
        
        std::cout << "Total Length "    << lcalThickness/dd4hep::cm  << " cm" << std::endl;
        std::cout << "Layer Thickness " << layerThickness/dd4hep::cm << " cm" << std::endl;
        
	// Initialisation for rInn1 and rOut1 for the first conical layer
	double rInn1 = rInnerStart+ 0.1*dd4hep::cm ;
	double rOut1 = rOuterStart-0.1*dd4hep::cm ;
	double rInn2 = 0 ;
	double rOut2 = 0 ; 

        //Loop for repeat=NN
        for(int i=0, repeat=xmlLayer.repeat(); i<repeat; ++i)  {
            
            std::string layer_name = detName + dd4hep::xml::_toString(thisLayerId,"_layer%d");
	    //Need to create layers from conical segments
            //DD4hep::Geometry::Tube layer_base(lcalInnerR,lcalOuterR,layerThickness*0.5);
	    //Definition of inner and outer end radii for each layer
	    rInn2 = rInn1 + tan(thetaInn)*layerThickness ;
	    rOut2 = rOut1 + tan(thetaOut)*layerThickness ;

	    std::cout << " Starting radii for layer " << i << " rinn " << rInn1 << " rout " << rOut1 << std::endl ;

	    ConeSegment layer_base(layerThickness*0.5,rInn1,rOut1,rInn2,rOut2,phi1,phi2);

	    std::cout << " angle of the layer " << tan((rOut2-rOut1)/layerThickness) << std::endl ;
            
            Volume layer_vol(layer_name,layer_base,air);
            
	    int sliceID=0;
            double inThisLayerPosition = -layerThickness*0.5;
            
            double nRadiationLengths=0.;
            double nInteractionLengths=0.;
            double thickness_sum=0;
            
            //DD4hep::DDRec::LayeredCalorimeterData::Layer caloLayer ;

	    double rSliceInn2, rSliceOut2 ;
	    double rSliceInn1 = rInn1;
	    double rSliceOut1 = rOut1;
            
            for(dd4hep::xml::Collection_t collSlice(xmlLayer,_U(slice)); collSlice; ++collSlice)  {
                dd4hep::xml::Component compSlice = collSlice;
                const double      slice_thickness = compSlice.thickness();
                const std::string sliceName = layer_name + dd4hep::xml::_toString(sliceID,"slice%d");
                Material   slice_material  = theDetector.material(compSlice.materialStr());
                

                //DD4hep::Geometry::Tube sliceBase(lcalInnerR,lcalOuterR,slice_thickness/2);
		rSliceInn2 = rSliceInn1 + tan(thetaInn)*slice_thickness ;
		rSliceOut2 = rSliceOut1 + tan(thetaOut)*slice_thickness ;

		ConeSegment sliceBase(slice_thickness/2.,rSliceInn1,rSliceOut1,rSliceInn2,rSliceOut2,phi1,phi2);
                
                Volume slice_vol (sliceName,sliceBase,slice_material);
                
                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;

                nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
                nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
                thickness_sum += slice_thickness/2;
                
                slice_vol.setAttributes(theDetector,compSlice.regionStr(),compSlice.limitsStr(),compSlice.visStr());
                layer_vol.placeVolume(slice_vol,Position(0,0,inThisLayerPosition+slice_thickness*0.5));
                    
                inThisLayerPosition += slice_thickness;
                ++sliceID;
		rSliceInn1 = rSliceInn2 ;
		rSliceOut1 = rSliceOut2 ;
            }//For all slices in this layer
	    
            //Why are we doing this for each layer, this just needs to be done once and then placed multiple times
            //Do we need unique IDs for each piece?
            layer_vol.setVisAttributes(theDetector,xmlLayer.visStr());
            
            Position layer_pos(0,0,referencePosition+0.5*layerThickness);
            referencePosition += layerThickness;
	    //if (i==1 || i==20 || i==39){
            PlacedVolume pv = envelopeVol.placeVolume(layer_vol,layer_pos);
            pv.addPhysVolID("layer",thisLayerId);
	    //}
	    rInn1 = rInn2 ;
	    rOut1 = rOut2 ;

            ++thisLayerId;
            
        }//for all layers
        
    }// for all layer collections




    const Position bcForwardPos (std::tan(0.5*fullCrossingAngle)*lcalCentreZ,0.0, lcalCentreZ);
    const Position bcBackwardPos(std::tan(0.5*fullCrossingAngle)*lcalCentreZ,0.0,-lcalCentreZ);
    const Rotation3D bcForwardRot ( RotationY(fullCrossingAngle*0.5 ) );
    const Rotation3D bcBackwardRot( RotationZYX ( (M_PI), (M_PI-fullCrossingAngle*0.5), (0.0)));

    PlacedVolume pv = envelope.placeVolume(envelopeVol, Transform3D( bcForwardRot, bcForwardPos ) );
    pv.addPhysVolID("barrel", 1);
    
    PlacedVolume pv2 = envelope.placeVolume(envelopeVol, Transform3D( bcBackwardRot, bcBackwardPos ) );
    pv2.addPhysVolID("barrel", 2);
    

    return sdet;
}
                                               
DECLARE_DETELEMENT(LumiCal_o2_v02,create_detector)
