#include "Systematics.h"


Systematics::Systematics(){}

Systematics::Systematics(std::unordered_map<std::string, PartStats> const &distats){

}
Systematics::~Systematics(){}

void Systematics::init(){

}

void Systematics::shiftJet(Particle& jet, TLorentzVector recoJet, std::string& syst_name, int syst){

  jet.addP4Syst(recoJet, syst);
  return;

}

void Systematics::shiftParticle(Particle& jet, TLorentzVector recoJet, double const& corrJetPt, double& corrJetMass, std::string& syst_name, int syst){

  TLorentzVector shiftedRecoJet;
  // Set the new components of the 4-momentum
  shiftedRecoJet.SetPtEtaPhiM(corrJetPt, recoJet.Eta(), recoJet.Phi(), corrJetMass);

  // std::cout << "I get here in shift particle" << std::endl;

  jet.addP4Syst(shiftedRecoJet, syst);
  return;

}

void Systematics::shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst){
  if (genLep == TLorentzVector(0,0,0,0)) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }
  double ratio = ((genLep.Pt()*scale) + (recoLep.Pt() - genLep.Pt())*resolution)/recoLep.Pt();
  //cout<<"ratio  "<<ratio<<"  "<<scale<<"  "<<resolution    <<std::endl;
   //add the shifted part up
   dPx+=recoLep.Px()*(ratio-1);
   dPy+=recoLep.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recoLep*=ratio;
   lepton.addP4Syst(recoLep, syst);
   return;
}


void Systematics::loadScaleRes(const PartStats& smear, const PartStats& syst, std::string syst_name) {
  scale = 1;
  resolution = 1;
  if(smear.bfind("SmearTheParticle")) {
    scale = smear.dmap.at("PtScaleOffset");
    resolution = smear.dmap.at("PtResolutionOffset");
  }
  if(syst_name.find("_Res_")) {
    resolution = syst_name.find("_Up") ? 1 + syst.dmap.at("res") : 1 - syst.dmap.at("res");
    scale=1;
  } else if(syst_name.find("_Scale_")) {
    scale = syst_name.find("_Up") ? 1+syst.dmap.at("scale") : 1- syst.dmap.at("scale");
    resolution=1;
  }
}

void Systematics::shiftLepton(Lepton& lepton,
                              TLorentzVector recoLep,
                              TLorentzVector genLep,
                              double& dPx,
                              double& dPy,
                              int syst,
                              std::string syst_name) {
  if (genLep == TLorentzVector(0.0, 0.0, 0.0, 0.0)) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }

  double scale_shift = 0.0;
  bool is_electron_in_transition_region = false;
  switch (lepton.type) {
    case PType::Electron: {
      const double eta = std::abs(recoLep.Eta());
      if (eta < 1.44) {
        // barrel
        scale_shift = 0.01;

      } else if ((eta > 1.57) and (eta < 2.50)) {
        // endcap
        scale_shift = 0.025;

      } else {
        // transition region or forward region
        scale_shift = 0.00;
        is_electron_in_transition_region = true;
      }
      break;
    }

    case PType::Muon: {
      scale_shift = 0.05;
      break;
    }

    case PType::Tau: {
      std::cerr << "it is deprecated to use Systematics::shiftLepton for taus in favour of Analyzer::smearTaus" << std::endl;
      break;
    }

    default: {
      std::cerr << "Systematics::shiftLepton got an unknown Lepton::type: " << static_cast<int>(lepton.type) << std::endl;
      assert(false);
    }
      
  }

  if (is_electron_in_transition_region) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }

  double scale = 1.0;
  double resolution = 1.0;
  if (syst_name.find("_Scale_")) {
    scale = syst_name.find("_Up") ? 1 + scale_shift : 1 - scale_shift;
  }

  double ratio = ((genLep.Pt() * scale) + (recoLep.Pt() - genLep.Pt()) * resolution) / recoLep.Pt();

   //add the shifted part up
   dPx += recoLep.Px() * (ratio - 1);
   dPy += recoLep.Py() * (ratio - 1);

   //WARNING change the particle content for the particle
   recoLep *= ratio;
   lepton.addP4Syst(recoLep, syst);

   return;
}
