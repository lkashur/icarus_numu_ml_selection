/**
 * @file csv_maker.h
 * @brief Header file defining a dummy SpillMultiVar for dumping particle info.
 * @author justin.mueller@colostate.edu
*/
#ifndef CSV_MAKER_H
#define CSV_MAKER_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "cuts.h"
#include "variables.h"

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

std::ofstream _selected("selection.log");
std::ofstream _signal("signal.log");
std::ofstream gamma_output("gamma_output.log");

#define GUARD(VAL) std::isinf(VAL) ? -9999 : VAL
#define OUT(STREAM,TAG) STREAM << std::fixed << TAG << ","
#define CSV(VAL) VAL << ","

/**
 * Writes information about signal events to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry. 
 */

const SpillMultiVar kSignal([](const caf::SRSpillProxy* sr)
{
  // Loop over true interactions
  for(auto const & ti : sr->dlp_true)
    {
      // 1mu1pi0
      if(cuts::signal_1mu1pi0(ti) && cuts::fiducial_containment_cut(ti))
	{
	  string eff_cat;
	  // Matched to selected reco interaction
	  if(cuts::matched(ti) && cuts::all_1muNph_cut(sr->dlp[ti.match[0]]))
	    {
	      eff_cat = "MATCH_SELECTED";
	    }
	  // Matched to other reco interaction
	  else if(cuts::matched(ti) && !cuts::all_1muNph_cut(sr->dlp[ti.match[0]]))
	    {
	      eff_cat = "MATCH_NOT_SELECTED";
	    }
	  // Not matched
	  else if(!cuts::matched(ti))
	    {
	      eff_cat = "NO_MATCH";
	    }

	  _signal << CSV(vars::image_id(ti))
		  << CSV(vars::id(ti))
		  << CSV(vars::neutrino_energy(ti))
		  << CSV(cuts::topology(ti))
		  << CSV(vars::pi0_leading_photon_energy(ti))
		  << CSV(vars::pi0_subleading_photon_energy(ti))
		  << CSV(vars::pi0_mass(ti))
		  << CSV(eff_cat)
		  << endl;

	} // end 1mu1pi0 signal loop
      
    } // end true interaction loop

  return std::vector<double>{1};

}
);


/**
 * Writes information about selected events to output file.
 * @param sr is an SRSpillProxy that attaches to the StandardRecord of the current spill.
 * @return a vector with a single dummy entry.
 */

const SpillMultiVar kSelected([](const caf::SRSpillProxy* sr)
{
  // Loop over reco interactions
  for(auto const & ri : sr->dlp)
    {
      // Passes selection
      if(cuts::all_1muNph_cut(ri))
	{
	  string pur_cat;
	  int category_topo;
	  // Matched
	  if(cuts::matched(ri))
	    {
	      const auto & ti = sr->dlp_true[ri.match[0]];
	      category_topo = vars::category_topology_1mu1pi0(ti);
	      if(vars::category_1mu1pi0(ti) == 0)
		{
		  pur_cat = "MATCH_SIGNAL";
		}
	      else
		{
		  pur_cat = "MATCH_NOT_SIGNAL";
		}
	    }
	  // Not matched
	  else
	    {
	      pur_cat = "NO_MATCH";
	      category_topo = -1;
	    }
	  
	  _selected << CSV(vars::image_id(ri))
		    << CSV(vars::id(ri))
		    << CSV(vars::visible_energy(ri))
		    << CSV(cuts::topology(ri))
		    << CSV(category_topo)
		    << CSV(vars::pi0_leading_photon_energy(ri))
		    << CSV(vars::pi0_subleading_photon_energy(ri))
		    << CSV(vars::pi0_mass(ri))
		    << CSV(pur_cat)
		    << endl;
	  
	} // end selected interaction loop
    } // end reco interaction loop

  return std::vector<double>{1};
}
);


#endif
