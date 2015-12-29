#include "info_pfc.h"

using namespace std;
using namespace fastjet;


MOD::InfoPFC::InfoPFC(int pdgId, string tag) : _pdgId(pdgId), _tag(tag) {}

const int MOD::InfoPFC::pdgId() const {
  return _pdgId;
}

const string MOD::InfoPFC::tag() const {
  return _tag;
}

