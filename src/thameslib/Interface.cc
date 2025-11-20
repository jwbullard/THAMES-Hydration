/**
@file Interface.cc
@brief Definition of methods for the Interface class.

*/
#include "Interface.h"

using std::vector;

bool cmp(const Site *s1, const Site *s2) { return s1->getWmc() < s2->getWmc(); }

bool affinitySort(const Isite s1, const Isite s2) {
  return s1.getAffinityInt() > s2.getAffinityInt();
}

Interface::Interface() {
  microPhaseId_ = 0;
  growthSites_.clear();
  dissolutionSites_.clear();
}

Interface::Interface(const bool verbose) {
  microPhaseId_ = 0;
  growthSites_.clear();
  dissolutionSites_.clear();

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif
}

Interface::Interface(ChemicalSystem *csys, vector<Site *> gv,
                     vector<Site *> dv, unsigned int pid, const bool verbose) {
  int j;
  int i;
  int aftyInt;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif

  microPhaseId_ = pid;
  chemSys_ = csys;

  affinityInt_.clear();
  affinityInt_ = chemSys_->getAffinityInt(microPhaseId_);
  int affSize = affinityInt_.size();
  cout << endl
       << "Interface::Interface - affinity for microPhaseId_ = "
       << microPhaseId_ << "   affinityInt_.size() = " << affSize
       << " : " << endl;
  for (int i = 0; i < affSize; i++) {
    cout << "   i = " << i << "    affinityInt_[" << i << "] = "
         << affinityInt_[i] << endl;
  }

  dissolutionSites_.clear();
  growthSites_.clear();

  ///
  /// create growth interface for microPhaseId_
  ///

  int gvsize = gv.size();
  for (j = 0; j < gvsize; j++) {
    aftyInt = 0;
    for (i = 0; i < NN_NNN; i++) {
      aftyInt += affinityInt_[gv[j]->nb(i)->getMicroPhaseId()];
    }

    ///
    /// Add to the list of Isites.  An Isite is an object consisting
    /// of a pointer to a site and an affinity value
    ///

    growthSites_.push_back(Isite(gv[j]->getId(), aftyInt));
  }

  ///
  /// create dissolution interface for microPhaseId_
  ///

  int dvsize = dv.size();
  for (j = 0; j < dvsize; j++) {
    aftyInt = 0;
    for (i = 0; i < NN_NNN; i++) {
      aftyInt += affinityInt_[dv[j]->nb(i)->getMicroPhaseId()];
    }
    dissolutionSites_.push_back(Isite(dv[j]->getId(), aftyInt));
  }

} // End of constructors

Interface::~Interface() {
  growthSites_.clear();
  dissolutionSites_.clear();
}

void Interface::addGrowthSite(Site *loc) {
  // vector<Isite>::iterator p, q, start, end;
  // start = growthSites_.begin();
  // end = growthSites_.end();

  // loc->setInGrowInterfacePos(microPhaseId_, growthSitesSize_);
  int aftyInt = 0;
  for (int i = 0; i < NN_NNN; i++) {
    aftyInt += affinityInt_[loc->nb(i)->getMicroPhaseId()];
  }
  Isite tisite(loc->getId(), aftyInt);
  // q = lower_bound(start, end, tisite, affinitySort);
  // growthSites_.insert(q, tisite);
  growthSites_.push_back(tisite);
}

void Interface::addDissolutionSite(Site *loc) {
  Isite tisite(loc->getId(), 0);
  dissolutionSites_.push_back(tisite);
}

void Interface::removeGrowthSite(int pos0, int pos1) {
  // if (pos0 != pos1)
  growthSites_[pos0] = growthSites_[pos1];
  growthSites_.pop_back();
}

void Interface::removeDissolutionSite(int pos0, int pos1) {
  // if (pos0 != pos1)
  // try {
  dissolutionSites_[pos0] = dissolutionSites_[pos1];
  dissolutionSites_.pop_back();
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Interface::removeDissolutionSite pos0/pos1 = " << pos0 <<" / " << pos1 << " => exit" << endl;
  //  exit(1);
  //}
}
