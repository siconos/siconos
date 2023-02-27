#pragma once

#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;
using ojson = nlohmann::ordered_json;

class BaseModel {
public:
  BaseModel();
  virtual ~BaseModel();

  void from_json(const json &j, const std::string &a_header);

protected:
  virtual void from_json(const json &j) = 0;
};
