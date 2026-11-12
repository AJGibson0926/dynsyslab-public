#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <string>

namespace dsl {

std::string lower(std::string value);
std::string token(const std::string& value);
bool has(const std::string& text, const std::string& part);
std::string property(const CatalogItem& item, const char* name);
std::string compactSystem(const std::string& system);
std::string titleWords(std::string value);
std::string compactPhysics(const CatalogItem& item);

}  // namespace dsl
