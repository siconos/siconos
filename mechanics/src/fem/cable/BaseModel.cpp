#include "BaseModel.h"

BaseModel::BaseModel() {}

BaseModel::~BaseModel() {}

void BaseModel::from_json(const json &j, const std::string &a_header)
{
  if (j.contains(a_header)) {
    from_json(j[a_header]);
  }
}
