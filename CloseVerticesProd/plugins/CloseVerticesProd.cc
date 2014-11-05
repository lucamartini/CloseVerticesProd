// -*- C++ -*-
//
// Package:    CloseVerticesProd/CloseVerticesProd
// Class:      CloseVerticesProd
// 
/**\class CloseVerticesProd CloseVerticesProd.cc CloseVerticesProd/CloseVerticesProd/plugins/CloseVerticesProd.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luca Martini
//         Created:  Wed, 05 Nov 2014 09:54:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include <Vertex.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "TVector3.h"
#include "TMath.h"
#include "vector"
#include <utility> 
#include <algorithm> 

//
// class declaration
//

using namespace edm;
using namespace reco;
using namespace std; 

const double mumass = .105658;

class CloseVerticesProd : public edm::EDProducer {
   public:
      explicit CloseVerticesProd(const edm::ParameterSet&);
      ~CloseVerticesProd();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
       InputTag PrimaryVertexCollectionTag_;
       EDGetTokenT<VertexCollection> PrimaryVertexCollectionToken_;
       
       InputTag DimuonVertexCollectionTag_;
       EDGetTokenT<VertexCollection> DimuonVertexCollectionToken_;

       const double MinCos_;
       unsigned int MaxPrimaryVerticesPerDimuon_;

};

struct sort_pred {
    bool operator()(const std::pair<Double_t, Vertex> &left, const std::pair<Double_t, Vertex> &right) {
        return left.first > right.first;
    }
};

//
// constructors and destructor
//
CloseVerticesProd::CloseVerticesProd(const edm::ParameterSet& iConfig) :
  PrimaryVertexCollectionTag_(iConfig.getParameter<edm::InputTag>("PrimaryVertexCollection")),
  PrimaryVertexCollectionToken_(consumes<VertexCollection>(PrimaryVertexCollectionTag_)),
  DimuonVertexCollectionTag_(iConfig.getParameter<edm::InputTag>("DimuonVertexCollection")),
  DimuonVertexCollectionToken_(consumes<VertexCollection>(DimuonVertexCollectionTag_)),
  MinCos_(iConfig.getUntrackedParameter<double>("MinCos", 0.99)),
  MaxPrimaryVerticesPerDimuon_(iConfig.getUntrackedParameter<unsigned int>("MaxPrimaryVerticesPerDimuon", 3))
{
    produces<VertexCollection>("DimuonPVs");  
}


CloseVerticesProd::~CloseVerticesProd()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called to produce the data  ------------
void
CloseVerticesProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get dimuon vertices
  VertexCollection DimuonVertexCollection;
  Handle<VertexCollection> DimuonVertexCollectionHandle;
  bool foundDimuonVertexCollection = iEvent.getByToken(DimuonVertexCollectionToken_, DimuonVertexCollectionHandle);
  if (foundDimuonVertexCollection) DimuonVertexCollection = *DimuonVertexCollectionHandle;
  else throw cms::Exception("CollectionNotFound") << DimuonVertexCollectionTag_.label() << " " << DimuonVertexCollectionTag_.instance() << " " << DimuonVertexCollectionTag_.process() << " does not exist" << endl;

  // get primary vertices
  VertexCollection PrimaryVertexCollection;
  Handle<VertexCollection> PrimaryVertexCollectionHandle;
  bool foundPrimaryVertexCollection = iEvent.getByToken(PrimaryVertexCollectionToken_, PrimaryVertexCollectionHandle);
  if (foundPrimaryVertexCollection) PrimaryVertexCollection = *PrimaryVertexCollectionHandle;
  else throw cms::Exception("CollectionNotFound") << PrimaryVertexCollectionTag_.label() << " " << PrimaryVertexCollectionTag_.instance() << " " << PrimaryVertexCollectionTag_.process() << " does not exist" << endl;
// cout << "dimuon vx = " << DimuonVertexCollection.size() << "   PVs = " << PrimaryVertexCollection.size() << endl;
  // make output
  auto_ptr<VertexCollection> vertexOutCollection(new VertexCollection());
  // loop on dimuons
  VertexCollection::iterator DimuonIt;
  VertexCollection::iterator eDimuonIt;
  for (DimuonIt = DimuonVertexCollection.begin(), eDimuonIt = DimuonVertexCollection.end(); DimuonIt != eDimuonIt; ++DimuonIt) {
    Vertex DimuonVertex = *DimuonIt;
    // check if the vertex actually consists of exactly two muon, throw exception if not
    //if (DimuonVertex.tracksSize() != 2) throw cms::Exception("BadLogic") << "the dimuon vertex must have exactly two muons by definition. Instead, it has " << DimuonVertex.tracksSize() << " tracks" << endl;

    TVector3 dimuonp3(DimuonVertex.p4(mumass).Px(), DimuonVertex.p4(mumass).Py(), DimuonVertex.p4(mumass).Pz());

    vector<pair<Double_t, Vertex> > v_CosAngle_vertex;

    // loop on PVs
    VertexCollection::iterator PrimaryIt;
    VertexCollection::iterator ePrimaryIt;
    for (PrimaryIt = PrimaryVertexCollection.begin(), ePrimaryIt = PrimaryVertexCollection.end(); PrimaryIt != ePrimaryIt; ++PrimaryIt) {
      Vertex PrimaryVertex = *PrimaryIt;

      TVector3 primaryp3(PrimaryVertex.p4().Px(), PrimaryVertex.p4().Py(), PrimaryVertex.p4().Pz());

      Double_t CosAngle = TMath::Cos(dimuonp3.Angle(primaryp3));
      if (CosAngle < MinCos_) continue;

      // make pairs of angles and PVs
      pair<Double_t, Vertex> CosAngle_vertex(CosAngle, PrimaryVertex);
      v_CosAngle_vertex.push_back(CosAngle_vertex);
    }

    // sort from angle = 0, up
    std::sort(v_CosAngle_vertex.begin(), v_CosAngle_vertex.end(), sort_pred());
// cout << endl;
// for (unsigned int i = 0; i < v_CosAngle_vertex.size(); i++) cout << v_CosAngle_vertex.at(i).first << endl;
// cout << endl;

    // put only the closest PVs
    for (unsigned int i = 0; i < MaxPrimaryVerticesPerDimuon_ && i < v_CosAngle_vertex.size(); i++) {
      vertexOutCollection->push_back(v_CosAngle_vertex.at(i).second);
    }    

  }

  iEvent.put(vertexOutCollection, "DimuonPVs");
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
CloseVerticesProd::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CloseVerticesProd::endJob() {
}

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CloseVerticesProd::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CloseVerticesProd);
