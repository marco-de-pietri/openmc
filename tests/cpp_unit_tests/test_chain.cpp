#include <catch2/catch_test_macros.hpp>
#include <pugixml.hpp>

#include "openmc/chain.h"

TEST_CASE("Breakup paths have no conventional MT product")
{
  pugi::xml_document doc;
  REQUIRE(doc.load_string(R"xml(<nuclide name="O16">
    <reaction type="(n,na)" target="C12" breakup="true"/>
    <reaction type="(n,gamma)" target="O17"/>
  </nuclide>)xml"));
  openmc::ChainNuclide parent(doc.document_element());
  const auto& products = parent.reaction_products();
  REQUIRE(products.count(22) == 0);
  REQUIRE(products.at(102).size() == 1);
  REQUIRE(products.at(102).front().name == "O17");
}

TEST_CASE("Li7 breakup does not register an MT33 product")
{
  pugi::xml_document doc;
  REQUIRE(doc.load_string(R"xml(<nuclide name="Li7">
    <reaction type="(n,nt)" target="He4" breakup="true"/>
  </nuclide>)xml"));
  openmc::ChainNuclide parent(doc.document_element());
  REQUIRE(parent.reaction_products().empty());
}
