/*
@file Site.cc
@brief Method definitions for the Site class.

*/

#include "Site.h"

using std::cout;
using std::endl;
using std::vector;

Site::Site() {
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  nbSA_.clear();
  stressFreeVolume_ = trueVolume_ = 0.0;
  damage_ = false;
  expstrain_ = 0.0;
}

Site::Site(int xp, int yp, int zp, int xs, int ys, int zs, int neigh,
           ChemicalSystem *csys, const bool verbose) {
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  nbSA_.clear();
  stressFreeVolume_ = trueVolume_ = 1.0;
  damage_ = false;
  x_ = xp;
  y_ = yp;
  z_ = zp;
  visit_ = 0;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif

  growth_.clear();

  id_ = x_ + (xs * y_) + ((xs * ys) * z_);

  if (neigh > 0)
    nb_.resize(neigh, 0);

  chemSys_ = csys;

  expstrain_ = 0.0;

  inGrowInterfacePos_.clear();
  inGrowInterfacePos_.resize(chemSys_->getNumMicroPhases(), -1);
  inDissInterfacePos_ = -1;

  inGrowthVectorPos_.clear();
  inGrowthVectorPos_.resize(chemSys_->getNumMicroPhases(), -1);
  inDissolutionVectorPos_ = -1;
}

// void Site::calcWmc(void) {
//   wmc_ = chemSys_->getMicroPhasePorosity(getMicroPhaseId());
//   for (int i = 0; i < NN_NNN; i++) {
//     wmc_ += chemSys_->getMicroPhasePorosity(nb_[i]->getMicroPhaseId());
//   }
// }

vector<int> Site::getXYZ() {
  vector<int> v(3, 0);
  v[0] = x_;
  v[1] = y_;
  v[2] = z_;
  return v;
}

void Site::setGrowthSite(int pid) {
  vector<int>::iterator start = growth_.begin();
  vector<int>::iterator end = growth_.end();
  vector<int>::iterator p = find(start, end, pid);
  if (p == growth_.end())
    growth_.push_back(pid);
}

void Site::removeGrowthSite(int pid) {
  bool found = false;
  int size = growth_.size();
  int i = -1;

  for (i = 0; i < size; i++) {
    if (growth_[i] == pid) {
      growth_[i] = growth_[size - 1];
      growth_.pop_back();
      found = true;
      break;
    }
  }
  if (found == false) {
    std::clog << endl << " stop - void removeGrowthSite(int pid) " << endl;
    std::clog.flush();
    std::clog << endl
              << "i size pid " << i << " " << size << " " << pid << endl;
    exit(1);
  }
}

void Site::setStressFreeVolume(double vol) {
  if (vol < 0) {
    std::clog
        << "in the setStrfreevolume function...volume should not be negative."
        << endl;
    std::cerr
        << "in the setStrfreevolume function...volume should not be negative."
        << endl;
    exit(1);
  } else {
    stressFreeVolume_ = vol;
  }
}

void Site::setTrueVolume(double vol) {
  if (vol < 0) {
    std::clog
        << "in the setTrueVolume function...volume should not be negative."
        << endl;
    std::cerr
        << "in the setTrueVolume function...volume should not be negative."
        << endl;
    exit(1);
  } else {
    trueVolume_ = vol;
  }
  return;
}
