/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/
#ifndef CUTS_H
#define CUTS_H

#include <functional>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>

#include <iostream>

namespace cuts
{
    /**
     * Apply a cut on whether a match exists.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to select on.
     * @return true if the object is matched.
    */
    template<class T>
        bool matched(const T & obj) { return obj.match.size() > 0; }

    /**
     * Apply a cut on the quality of the reconstruction using the overlap
     * fraction as a discriminant.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to selection on.
     * @return true if the object is well reconstructed.
    */
    template<class T>
        bool wellreco(const T & obj) { return matched(obj) && obj.match[0] > 0.9; }

    /**
     * Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
    */
    template<class T>
        bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.fmatched == 1; }

    /**
     * Check if the particle meets final state signal requirements.
     * Particles must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param particle to check.
     * @return true if the particle is a final state signal particle.
    */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
                if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
			       energy = p.energy_init; // p.energy_deposit

                if((p.pid == 2 && energy > 143.425) || (p.pid < 2 && energy > 40))
                    passes = true;
            }
            return passes;
        }

    /**
     * Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & interaction)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : interaction.particles)
            {
                if(final_state_signal(p))
                    ++counts[p.pid];
            }
            return counts;
        }

    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
        std::string topology(const T & interaction)
        {
            std::vector<uint32_t> counts(count_primaries(interaction));
            std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
        }

    /**
     * Apply no cut (all interactions/particles passed).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to select on.
     * @return true (always).
     */
    template<class T>
        bool no_cut(const T & obj) { return true; }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & interaction)
        {
            return interaction.is_fiducial && !(interaction.vertex[0] > 210.215 && interaction.vertex[1] > 60 && (interaction.vertex[2] > 290 && interaction.vertex[2] < 390));
        }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
        bool containment_cut(const T & interaction) { return interaction.is_contained; }

    /**
     * Apply a 1mu1p topological cut. The interaction must have a topology
     * matching 1mu1p as defined by the conditions in the topology() function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
     */
    template<class T>
        bool topological_1mu1p_cut(const T & interaction) { return topology<T>(interaction) == "0ph0e1mu0pi1p"; }
    
    /**
     * Apply a 1muNp topological cut. The interaction must have a topology
     * matching 1muNp as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muNp topology.
     */
    template<class T>
        bool topological_1muNp_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] >= 1;
        }

    /**
     * Apply a 1muX topological cut. The interaction must have a topology
     * matching 1muX as defined by the conditions in the count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muX topology.
     */
    template<class T>
        bool topological_1muX_cut(const T & interaction)
        {
            std::vector<uint32_t> c(count_primaries(interaction));
            return c[2] == 1;
        }

    /**
     * Apply a 1muNph topological cut (N > 1).  The interaction must have a topology
     * matching 1muNph as defined by the conditions in the count_primaries() function.
     * @tparm T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1muNph topology.
     */
    template<class T>
      bool topological_1muNph_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
	return c[2] == 1 && c[0] > 1;
      }

    /**
     * Define the topological 1mu1pi0 interaction classification.
     * @tparam T the type of interactino (true or reco).
     * @param interaction to select on.
     * @return true if interactin is a 1mu1pi0 neutrino interaction.
     */
    template<class T>
      bool topological_1mu1pi0_cut(const T & interaction)
      {

        int primary_muon_count = 0;
        int primary_pi0_count = 0;
	int cex_pi0_count = 0;
        unordered_map< int, vector<float> > primary_pi0_map;
	unordered_map< int, vector<float> > cex_pi0_map;
        // Loop over particles
        for(auto & p : interaction.particles)
          {
            // Primary Muon
            if(p.pid == 2 && p.is_primary && p.energy_init > 143.425)
              {
                ++primary_muon_count;
              }
            // Group primary pi0 photons
            if(p.pid == 0 & p.is_primary && p.ancestor_pdg_code == 111)
              {
                primary_pi0_map[p.ancestor_track_id].push_back(p.energy_init);
              }
          }

        // Loop over primary pi0s   
        for (auto const& pi0 : primary_pi0_map)
          {
            int num_primary_photon_daughters = 0;
            for(auto & e : pi0.second)
              {
                if(e > 40) ++num_primary_photon_daughters;
              }

            if(num_primary_photon_daughters == 2)
              {
                ++primary_pi0_count;
              }
          }
	
	// Loop over CEx pi0s
	for (auto const& pi0 : cex_pi0_map)
	  {
            int num_cex_photon_daughters = 0;
            for(auto & e : pi0.second)
              {
                ++num_cex_photon_daughters;
              }

            if(num_cex_photon_daughters == 2)
              {
                ++cex_pi0_count;
              }
          }

	if (primary_muon_count == 1 && primary_pi0_count == 1 && cex_pi0_count == 0)
          {
            return true;
          }
        else
          {
            return false;
          }
      }

    
    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
        }

    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB data.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut_data(const T & interaction)
        {
            if(!valid_flashmatch(interaction))
                return false;
            else
                return (interaction.flash_time >= -0.5) && (interaction.flash_time <= 1.4);
        }

    /**
     * Apply the CRT-PMT veto cut. The event must not have a CRT-PMT
     * match for any flash in the beam gate.
     * @param sr the spill to check for CRT-PMT matches.
     * @return true if the event passes the CRT-PMT veto.
    */
    bool crtpmt_veto(const caf::SRSpillProxy * sr)
    {
        bool crtpmt_matched(sr->ncrtpmt_matches == 0);
        for(auto const & c : sr->crtpmt_matches)
        {
            if(c.flashGateTime > 0 && c.flashGateTime < 1.6 && c.flashClassification == 0)
                crtpmt_matched = true;
        }
        return crtpmt_matched;
    }

    /**
     * Apply the CRT-PMT veto cut (data). The event must not have a CRT-PMT
     * match for any flash in the beam gate.
     * @param sr the spill to check for CRT-PMT matches.
     * @return true if the event passes the CRT-PMT veto.
    */
    bool crtpmt_veto_data(const caf::SRSpillProxy * sr)
    {
        bool crtpmt_matched(sr->ncrtpmt_matches == 0);
        for(auto const & c : sr->crtpmt_matches)
        {
            if(c.flashGateTime > -0.5 && c.flashGateTime < 1.4 && c.flashClassification == 0)
                crtpmt_matched = true;
        }
        return crtpmt_matched;
    }

    /**
     * Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1mu1p) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1mu1p_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1mu1p_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1muNp) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1muNp_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1muNp_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1muX) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_1muX_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1muX_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1mu1p), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1mu1p_cut(const T & interaction) { return topological_1mu1p_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNp), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1muNp_cut(const T & interaction) { return topological_1muNp_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }
    
    /**
     * Apply a fiducial, containment, topological (1muX), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1muX_cut(const T & interaction) { return topological_1muX_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1mu1p), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1mu1p_data_cut(const T & interaction) { return topological_1mu1p_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNp), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1muNp_data_cut(const T & interaction) { return topological_1muNp_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && containment_cut<T>(interaction); }
    
    /**
     * Apply a fiducial, containment, topological (1muX), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_1muX_data_cut(const T & interaction) { return topological_1muX_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological (1muNph) cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
      bool fiducial_containment_topological_1muNph_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_1muNph_cut<T>(\
																						    interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNph), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
      bool all_1muNph_cut(const T & interaction) { return topological_1muNph_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, topological (1muNph), and flash time cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
      bool all_1muNph_data_cut(const T & interaction) { return topological_1muNph_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut_data<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Defined the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool neutrino(const T & interaction) { return interaction.is_neutrino; }

    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
        bool cosmic(const T & interaction) { return !interaction.is_neutrino; }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_neutrino(const T & interaction) { return interaction.match.size() > 0 && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification (well-reconstructed).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a well-reconstructed neutrino interaction.
    */
    template<class T>
        bool wellreco_neutrino(const T & interaction) { return wellreco(interaction) && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_cosmic(const T & interaction) { return interaction.match.size() > 0 && cosmic(interaction); }

    /**
     * Define the true 1mu1p interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1mu1p neutrino interaction.
     */
    template<class T>
        bool signal_1mu1p(const T & interaction) { return topological_1mu1p_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muNp interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNp neutrino interaction.
     */
    template<class T>
        bool signal_1muNp(const T & interaction) { return topological_1muNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muNp interaction classification (N > 1 strictly).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNp (N > 1) neutrino interaction.
     */
    template<class T>
        bool signal_1muNp_Nnot1(const T & interaction) { return !topological_1mu1p_cut(interaction) && topological_1muNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muX interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muX neutrino interaction.
     */
    template<class T>
        bool signal_1muX(const T & interaction) { return topological_1muX_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1muX interaction classification (not 1muNp).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muX (not 1muNp) neutrino interaction.
     */
    template<class T>
        bool signal_1muX_notNp(const T & interaction) { return !topological_1muNp_cut(interaction) && topological_1muX_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true 1mu1pi0 interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1muNp neutrino interaction.
     */
    template<class T>
      bool signal_1mu1pi0(const T & interaction) { return topological_1mu1pi0_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true "other neutrino" interaction classification (1mu1p).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1mu1p(const T & interaction) { return !topological_1mu1p_cut<T>(interaction) && neutrino(interaction); }
    
    /**
     * Define the true "other neutrino" interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
      bool other_nu_1mu1pi0(const T & interaction) { return !topological_1mu1pi0_cut<T>(interaction) && neutrino(interaction); }
    
    /**
     * Define the true "other neutrino" interaction classification (1muNp).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1muNp(const T & interaction) { return !topological_1muNp_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define the true "other neutrino" interaction classification (1muX).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu_1muX(const T & interaction) { return !topological_1muX_cut<T>(interaction) && neutrino(interaction); }

    /**
     * Define true muon particle classification.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon.
    */
    template<class T>
        bool muon(const T & particle) { return particle.pid == 2 && particle.is_contained; }

    /**
     * Define true muon particle classification (matched).
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon (matched).
    */
    template<class T>
        bool matched_muon(const T & particle) { return muon(particle) && matched(particle); }

    /**
     * Define true proton particle classification.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a proton.
    */
    template<class T>
        bool proton(const T & particle) { return particle.pid == 4 && particle.is_contained; }

    /**
     * Define true proton particle classification (matched).
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a proton (matched).
    */
    template<class T>
        bool matched_proton(const T & particle) { return proton(particle) && matched(particle); }

    /**
     * Define cut for particles crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is cathode-crossing.
    */
    template<class T>
        bool cathode_crossing(const T & particle)
        {
            bool cc(!std::isnan(particle.start_point[0]) && !std::isnan(particle.end_point[0]));
            if(cc && particle.volume_id == 0)
            {
                int p0(std::copysign(1, particle.start_point[0] + 210.215));
                int p1(std::copysign(1, particle.end_point[0] + 210.215));
                cc = p0 != p1;
            }
            else if(cc && particle.volume_id == 1)
            {
                int p0(std::copysign(1, particle.start_point[0] - 210.215));
                int p1(std::copysign(1, particle.end_point[0] - 210.215));
                cc = p0 != p1;
            }
            return cc;
        }

    /**
     * Define cut for muons crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is a muon that crosses the cathode.
    */
    template<class T>
        bool cathode_crossing_muon(const T & particle) { return particle.pid == 2 && cathode_crossing(particle); }
        
    /**
     * Define cut for muons not crossing the cathode.
     * @tparam T the type of particle (true or reco).
     * @param particle to select on.
     * @return true if the particle is a muon that does not cross the cathode.
    */
    template<class T>
        bool non_cathode_crossing_muon(const T & particle) { return particle.pid == 2 && !cathode_crossing(particle); }
    
    /**
     * Define muons contained to a single TPC.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a muon contained to a single TPC.
    */
    template<class T>
        bool contained_tpc_muon(const T & particle) { return muon(particle) && !cathode_crossing(particle); }

    /**
     * Define muons contained that are well-reconstructed.
     * @tparam T the type of particle (true or reco)
     * @param particle to select on.
     * @return true if the particle is a contained muon that is
     * well-reconstructed.
    */
    template<class T>
        bool wellreco_muon(const T & particle) { return muon(particle) && wellreco(particle); }
}
#endif
