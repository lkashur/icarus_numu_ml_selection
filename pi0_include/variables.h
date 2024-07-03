/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>
#include <iostream>
#include <TVector3.h>

namespace vars
{

    /**
     * Variable for counting interactions/particles.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & obj) { return 1.0; }

    /**
     * Variable for image_id (unique identifier for the event).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the image_id of the interaction/particle.
    */
    template<class T>
        double image_id(const T & obj) { return obj.image_id; }

    /**
     * Variable for id (unique identifier for the object).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the id of the interaction/particle.
    */
    template<class T>
        double id(const T & obj) { return obj.id; }

    /**
     * Variable for the cryostat of the object.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the cryostat of the interaction/particle.
    */
    template<class T>
        double cryostat(const T & obj) { return obj.volume_id; }

    /**
     * Variable for enumerating interaction categories.  This is a basic
     * categorization using only signal, neutrino background, and cosmic
     * background as the three categories.
     * 0: 1mu1pi0 (contained and fiducial)
     * 1: 1mu1pi0 (not contained or fiducial)               
     * 2: Other nu                
     * 3: cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.   
     */
    template<class T>
      double category(const T & interaction)
      {
        // Cosmic background               
        double cat(3);

        // Signal
        if(cuts::signal_1mu1pi0(interaction))
          {
            if(cuts::fiducial_containment_cut(interaction)) cat = 0;
            else cat = 1;
          }
        // Neutrino Background                         
        else if(cuts::other_nu_1mu1pi0(interaction))
          {
            cat = 2;
          }
        return cat;
      }

    /**
     * Variable for enumerating interaction categories. This classifies the          
     * interactions based on the visible final states.
     * 0: 1mu1pi0, 1: 1muNpi0, 2: 1muCex, 3: NC, 4: Other, 5: Cosmic        
     * @tparam T the type of interaction (true or reco).             
     * @param interaction to apply the variable on.      
     * @return the enumerated category of the interaction.                                               
     */
    template<class T>
      double category_topology(const T & interaction)
      {
        // Cosmic                              
        uint16_t cat(5);
        if(interaction.is_neutrino)
          {

            int primary_muon_count = 0;
            unordered_map< int, vector<double> > primary_pi0_map;
            unordered_map< int, vector<double> > cex_pi0_map;
            int primary_pi0_count = 0;
            int cex_pi0_count = 0;

            // Loop over particles                                   
            for(auto & p : interaction.particles)
              {
                // Muon count                    
                if(p.pid == 2 && p.is_primary && p.energy_init > 143.425)
                  {
                    primary_muon_count++;
                  }
                // Group primary pi0 photons                        
                if(p.pid == 0 && p.is_primary && p.ancestor_pdg_code == 111)
                  {
                    primary_pi0_map[p.ancestor_track_id].push_back(p.energy_init);
                  }
                // Group CEx pi0 photons                           
                if(p.pid == 0 && !p.is_primary && abs(p.ancestor_pdg_code) == 211)
                  {
                    cex_pi0_map[p.ancestor_track_id].push_back(p.energy_init);
                  }
              }

            // Loop over primary pi0s                               
            for (auto const& pi0 : primary_pi0_map)
              {
                int num_primary_photon_daughters = 0;
                for(auto e : pi0.second)
                  {
                    if(e > 25) num_primary_photon_daughters++;
                  }

                if(num_primary_photon_daughters == 2)
                  {
                    primary_pi0_count++;
                  }
              }

            // Loop over CEx pi0s                                   
            for (auto const& pi0 : cex_pi0_map)
              {
                int num_cex_photon_daughters = 0;
                for(auto e : pi0.second)
                  {
                    num_cex_photon_daughters++;
                  }

                if(num_cex_photon_daughters == 2)
                  {
                    cex_pi0_count++;
                  }
              }

            // Categorization                                                         
            if(primary_muon_count == 1 && primary_pi0_count == 1 && cex_pi0_count == 0 && interaction.nu_current_type == 0 && interaction.is_contained && interaction.is_fiducial) cat =0;
            else if(primary_muon_count == 1 && primary_pi0_count > 1 && cex_pi0_count == 0 && interaction.nu_current_type == 0) cat = 1;
            else if(primary_muon_count == 1 && cex_pi0_count >= 1 && interaction.nu_current_type == 0) cat = 2;
            else if(interaction.nu_current_type == 1) cat = 3;
            else(cat) = 4;
          }
        return cat;
      }

    /**
     * Variable for enumerating interaction categories. This categorization
     * uses the interaction type (generator truth) classify the interactions
     * 0: nu_mu CC QE, 1: nu_mu CC Res, 2: nu_mu CC MEC, 3: nu_mu CC DIS, 4: nu_mu CC Coh, 5: nu_e CC, 6: NC, 7: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_interaction_mode(const T & interaction)
        {
            double cat(7);

            if(interaction.is_neutrino)
            {
                if(interaction.nu_current_type == 0)
                {
                    if(abs(interaction.nu_pdg_code) == 14)
                    {
                        if(interaction.nu_interaction_mode == 0) cat = 0;
                        else if(interaction.nu_interaction_mode == 1) cat = 1;
                        else if(interaction.nu_interaction_mode == 10) cat = 2;
                        else if(interaction.nu_interaction_mode == 2) cat = 3;
                        else if(interaction.nu_interaction_mode == 3) cat = 4;
                        else cat = 8;
                    }
                    else cat = 5;
                }
                else cat = 6;
            }

            return cat;
        }

    /**
     * Variable for counting particles in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of particles in the interaction.
     */
    template<class T>
        double count_particles(const T & interaction) { return interaction.num_particles; }

    /**
     * Variable for counting primaries in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of primaries in the interaction.
     */
    template<class T>
        double count_primaries(const T & interaction) { return interaction.num_primaries; }
    
    /**
     * Variable for energy of the neutrino primary of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the neutrino energy.
    */
    template<class T>
        double neutrino_energy(const T & interaction) { return 1000*interaction.nu_energy_init; }

    /**
     * Variable for matched interaction flash time.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the matched flash time of the interaction.
    */
    template<class T>
        double flash_time(const T & interaction)
        {
            if(!cuts::valid_flashmatch(interaction))
                return -100000.0;
            else
                return interaction.flash_time;
        }

    /**
     * Variable for particle primary categorizations.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the primary/non-primary designation of the particle.
    */
    template<class T>
        double primary(const T & particle) { return particle.is_primary ? 1 : 0; }

    /**
     * Variable for particle PID.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID of the particle.
    */
    template<class T>
        double pid(const T & particle) { return particle.pid; }

    /**
     * Variable for particle PID + primary.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID+primary information for the particle.
    */
    template<class T>
        double primary_pid(const T & particle) { return particle.pid + (particle.is_primary ? 5 : 0); }

    /**
     * Variable for particle csda_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle.
    */
    template<class T>
        double csda_ke(const T & particle) { return particle.csda_ke; }

    /**
     * Variable for particle calo_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the calo_ke of the particle.
    */
    template<class T>
        double calo_ke(const T & particle) { return particle.calo_ke; }

    /**
     * Variable for particle csda_ke (muons only).
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle (if a muon).
    */
    template<class T>
        double csda_ke_muon(const T & particle) { return (cuts::muon(particle)) ? csda_ke(particle) : -1; }

    /**
     * Variable for true particle energy deposited.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the energy deposited by the particle (truth only).
    */
    template<class T>
        double energy_deposit(const T & particle) { return particle.energy_deposit; }

    /**
     * Variable for true particle energy (total).
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the total energy of the paricle (truth only).
     */
    template<class T>
      double energy_init(const T & particle) {return particle.energy_init; }
    
    /**
     * Variable for true particle energy starting kinetic energy.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
    */
    template<class T>
        double ke_init(const T & particle)
        {
            double energy(particle.energy_init);
            switch (particle.pid)
            {
            case 1:
                energy -= ELECTRON_MASS;
                break;
            case 2:
                energy -= MUON_MASS;
                break;
            case 3:
                energy -= PION_MASS;
                break;
            case 4:
                energy -= PROTON_MASS;
                break;
            default:
                break;
            }
            return energy;
        }

    /**
     * Find indices for true pi0 decay photons.
     * @tparam T the type of interaction (true).
     * @param the interaction to operate on.
     * @return the indices of the true pi0 decay photons.
     */
    
    template <class T>
      vector<size_t> true_pi0_photon_idxs(const T & interaction)
      {
	// Output 
	vector<size_t> photon_idxs;
	
	// Temp. Storage
	unordered_map<int, vector<pair<size_t, double>> > primary_pi0_map;
	vector< pair<size_t, double> > photon_energies;

	// 1mu1pi0 signal
	if(cuts::signal_1mu1pi0(interaction))
	  {
	    // Loop over particles
	    for(size_t i(0); i < interaction.particles.size(); ++i)
	      {
		const auto & p = interaction.particles[i];

		// Primary pi0
		// Given 1mu1pi0 cut, we already know that any photons 
		// meeting the below criteria belong to a single primary pi0
		if(p.pid == 0 && p.is_primary && p.energy_init > 40 && p.ancestor_pdg_code == 111)
		  {
		    primary_pi0_map[p.ancestor_track_id].push_back({i, p.energy_init}); 
		  }
	      } // end particle loop
	    
	    // Primary pi0s
	    for (auto const& pi0 : primary_pi0_map)
	      {
		int num_primary_photon_daughters = 0;
		for(auto pair : pi0.second)
		  {
		    if(pair.second > 40) ++num_primary_photon_daughters;
		  }
		
		if(num_primary_photon_daughters == 2)
		  {
		    for(auto pair : pi0.second)
		      {
			if(pair.second > 40)
			  {
			    photon_energies.push_back(pair);
			  }
		      }
		  }
	      }

	    // Sort by photon energy (high to low)
	    sort(photon_energies.begin(), photon_energies.end(), [](auto& a, auto& b) {
		return a.second > b.second; });

	    // Fill idxs
	    photon_idxs.push_back(photon_energies[0].first);
	    photon_idxs.push_back(photon_energies[1].first); 
	  } // end 1mu1pi0 signal
	// Fill with nonsense if not a true 1mu1pi0 event
	else
	  {
	    photon_idxs.push_back(0);
            photon_idxs.push_back(0);
	  }
	return photon_idxs;
      }
    

    /**
     * Find true pi0 photon directions.
     * @tparam T the type of interaction (true).
     * @param interaction to operate on.
     * @return the true pi0 photon directions.
     */    
    template<class T>
      vector<TVector3> true_pi0_photon_dirs(const T & interaction)
      {
	// Output
	vector<TVector3> photon_dirs;
	
	vector<size_t> pi0_photon_idxs;
	pi0_photon_idxs = true_pi0_photon_idxs(interaction);
	size_t i(pi0_photon_idxs[0]);
	size_t j(pi0_photon_idxs[1]);

	// Photon 0
	const auto & p = interaction.particles[i];
	TVector3 ph0_dir (p.truth_momentum[0], p.truth_momentum[1], p.truth_momentum[2]);
	ph0_dir = ph0_dir.Unit();
	photon_dirs.push_back(ph0_dir);

	// Photon 1
	const auto & q = interaction.particles[j];
	TVector3 ph1_dir (q.truth_momentum[0], q.truth_momentum[1], q.truth_momentum[2]);
	ph1_dir = ph1_dir.Unit();
	photon_dirs.push_back(ph1_dir);

	return photon_dirs;
      }
    
    /**
     * Find indicies of pi0 decay candidate photons.
     * Assumes 1muNph cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the indices of the pi0 decay candidate photons.
     */
    template<class T>
      vector<size_t> pi0_photon_idxs_by_dir(const T & interaction)
      {
	// Output
	vector<size_t> photon_idxs;

	if(cuts::all_1muNph_cut(interaction))
	  {

	    // Temp
	    vector<pair< pair<size_t, size_t>, float> > pair_angles;
	    
	    // Loop over particles
	    for(size_t i(0); i < interaction.particles.size(); ++i)
	      {
		const auto & p = interaction.particles[i];
		if (p.pid != 0 || !p.is_primary || p.calo_ke < 40.0) continue;

		TVector3 sh0_start(p.start_point[0], p.start_point[1], p.start_point[2]);
		TVector3 sh0_dir(p.start_dir[0], p.start_dir[1], p.start_dir[2]);

		for(size_t j(0); j < interaction.particles.size(); ++j)
		  {
		    const auto & q = interaction.particles[j];
		    if(i == j) continue;
		    if(q.pid != 0 || !q.is_primary || q.calo_ke < 40.0) continue;
		    
		    // Shower start point and direction
		    TVector3 sh1_start(q.start_point[0], q.start_point[1], q.start_point[2]);
		    TVector3 sh1_dir(q.start_dir[0], q.start_dir[1], q.start_dir[2]);

		    // Find "vertex" of each photon pair
		    TVector3 c0;
		    TVector3 c1;
		    TVector3 n = sh0_dir.Cross(sh1_dir);
		    TVector3 n0 = sh0_dir.Cross(n);
		    TVector3 n1 = sh1_dir.Cross(n);
		    float s0 = (sh1_start - sh0_start).Dot(n1) / sh0_dir.Dot(n1);
		    float s1 = (sh0_start - sh1_start).Dot(n0) / sh1_dir.Dot(n0);

		    if (s0 > 0 && s1 > 0)
		      {
			c0 = sh0_start;
			c1 = sh1_start;
		      }
		    else if (s0 > 0 && s1 < 0)
		      {
			c0 = sh0_start;
			c1 = sh0_start;
		      }
		    else if (s0 < 0 && s1 > 0)
		      {
			c0 = sh1_start;
			c1 = sh1_start;
		      }
		    else
		      {
			c0 = sh0_start + s0*sh0_dir;
			c1 = sh1_start + s1*sh1_dir;
		      }

		    float d0 = (sh0_start - c0).Mag();
		    float d1 = (sh1_start - c1).Mag();

		    TVector3 vertex;
		    if (d0 == 0 || d1 == 0)
		      {
			float vertex_x = (c0[0] + c1[0]) / 2;
			float vertex_y = (c0[1] + c1[1]) / 2;
			float vertex_z = (c0[2] + c1[2]) / 2;
			vertex.SetX(vertex_x);
			vertex.SetY(vertex_y);
			vertex.SetZ(vertex_z);
		      }
		    else
		      {
			float vertex_x = ((c0[0] * d1) + (c1[0] * d0)) / (d1 + d0);
			float vertex_y = ((c0[1] * d1) + (c1[1] * d0)) / (d1 + d0);
			float vertex_z = ((c0[2] * d1) + (c1[2] * d0)) / (d1 + d0);
			vertex.SetX(vertex_x);
			vertex.SetY(vertex_y);
			vertex.SetZ(vertex_z);
		      }

		    // Find mean angular displacement between
		    // <sh_start from vertex> and <sh_dir>
		    float r0 = (sh0_start - vertex).Mag();
		    float r1 = (sh1_start - vertex).Mag();
		    float angle = 0;
		    if (r0 > 0)
		      {
			TVector3 v0 = (sh0_start - vertex).Unit();
			angle += acos( sh0_dir.Dot(v0) )/2;
		      }
		    if (r1 > 0)
		      {
			TVector3 v1 = (sh1_start - vertex).Unit();
			angle += acos( sh1_dir.Dot(v1) )/2;
		      }

		    pair_angles.push_back(make_pair(make_pair(i, j), angle));
		    sort(pair_angles.begin(), pair_angles.end(), [](const pair<pair<int, int>, float> &a, const pair<pair<int, int>, float> &b)
			 { return a.second < b.second;});

		  } // end second particle loop
	      } // end first particle loop
	    
	    // Sort by energy
	    if(calo_ke(interaction.particles[pair_angles[0].first.first]) > calo_ke(interaction.particles[pair_angles[0].first.second]))
	      {
		photon_idxs.push_back(pair_angles[0].first.first);
		photon_idxs.push_back(pair_angles[0].first.second);
	      }
	    else
	      {
		photon_idxs.push_back(pair_angles[0].first.second);
		photon_idxs.push_back(pair_angles[0].first.first);
	      }
	  }
	// If not 1muNph, return nonsense
	else
	  {
	    photon_idxs.push_back(0);
            photon_idxs.push_back(0);
	  }

	return photon_idxs;
      }

    /**
     * Find indices for pi0 decay candidate photons.
     * Assumes 1muNph cut has been made.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the indices of the pi0 decay candidate photons.
     */
    template<class T>
      vector<size_t> pi0_photon_idxs_by_energy(const T & interaction)
      {
	// Output
	vector<size_t> photon_idxs;

	if(cuts::all_1muNph_cut(interaction))
	  {
	    // Temp. storage
	    vector< pair<size_t, double> > photon_energies;

	    // Loop over particles
	    for(size_t i(0); i < interaction.particles.size(); ++i)
	      {
		const auto & p = interaction.particles[i];
		if(p.pid == 0 && p.is_primary && p.calo_ke > 40)
		  {
		    photon_energies.push_back({i, p.calo_ke});
		  }
	      }

	    // Sort by photon energy (high to low)
	    sort(photon_energies.begin(), photon_energies.end(), [](auto& a, auto& b) {
		return a.second > b.second; });
	    
	    photon_idxs.push_back(photon_energies[0].first);
            photon_idxs.push_back(photon_energies[1].first);
	    
	  }
	// If not 1muNph, return nonsense
	else
	  {
	    photon_idxs.push_back(0);
	    photon_idxs.push_back(0);
	  }
	
	return photon_idxs;
      }
    
    /**
     * Find pi0 photon directions.
     * @tparam the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the true pi0 photon directions.
     */
    
    template<class T>
      vector<TVector3> pi0_photon_dirs(const T & interaction)
      {
        // Output
        vector<TVector3> photon_dirs;

        TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);

        vector<size_t> pi0_photon_idxs;
        //pi0_photon_idxs = pi0_photon_idxs_by_energy(interaction);
	pi0_photon_idxs = pi0_photon_idxs_by_dir(interaction);
        size_t i(pi0_photon_idxs[0]);
        size_t j(pi0_photon_idxs[1]);

        // Photon 0
        const auto & p = interaction.particles[i];
        TVector3 ph0_dir;
        ph0_dir.SetX(p.start_point[0] - vertex[0]);
        ph0_dir.SetY(p.start_point[1] - vertex[1]);
        ph0_dir.SetZ(p.start_point[2] - vertex[2]);
        ph0_dir = ph0_dir.Unit();
        photon_dirs.push_back(ph0_dir);

        // Photon 1 
        const auto & q = interaction.particles[j];
        TVector3 ph1_dir;
        ph1_dir.SetX(q.start_point[0] - vertex[0]);
        ph1_dir.SetY(q.start_point[1] - vertex[1]);
        ph1_dir.SetZ(q.start_point[2] - vertex[2]);
        ph1_dir = ph1_dir.Unit();
        photon_dirs.push_back(ph1_dir);

        return photon_dirs;
      }
    

    /**
     * Find cosine of pi0 opening angle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the cosine of pi0 opening angle.
     */
    
    template<class T>
      double pi0_costheta(const T & interaction)
      {
        // Output
        double costheta;

        vector<TVector3> photon_dirs;
        // Truth
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         photon_dirs = true_pi0_photon_dirs(interaction);
                       }

        // Reco 
        else
          {
            photon_dirs = pi0_photon_dirs(interaction);
          }

        costheta = photon_dirs[0].Dot(photon_dirs[1]);

        return costheta;
      }
    

    /**
     * Find leading pi0 photon energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the leading pi0 photon energy.
     */
    
    template<class T>
      double pi0_leading_photon_energy(const T & interaction)
      {
	vector<size_t> pi0_photon_idxs;
	double energy(0);
	// Truth
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			 pi0_photon_idxs = true_pi0_photon_idxs(interaction);
			 size_t i(pi0_photon_idxs[0]);
			 energy = energy_init(interaction.particles[i]);
		       }
	// Reco
	else
	  {
	    //pi0_photon_idxs = pi0_photon_idxs_by_energy(interaction);
	    pi0_photon_idxs = pi0_photon_idxs_by_dir(interaction);
	    size_t i(pi0_photon_idxs[0]);
	    energy = calo_ke(interaction.particles[i]);
	  }
	return energy;		       
      }
    

    /**
     * Find subleading pi0 photon energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the subleading pi0 photon energy.
     */
    
    template<class T>
      double pi0_subleading_photon_energy(const T & interaction)
      {
	vector<size_t> pi0_photon_idxs;
        double energy(0);
	// Truth
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         pi0_photon_idxs = true_pi0_photon_idxs(interaction);
                         size_t i(pi0_photon_idxs[1]);
                         energy = energy_init(interaction.particles[i]);
                       }
	// Reco
	else
          {
            //pi0_photon_idxs = pi0_photon_idxs_by_energy(interaction);
	    pi0_photon_idxs = pi0_photon_idxs_by_dir(interaction);
            size_t i(pi0_photon_idxs[1]);
            energy = calo_ke(interaction.particles[i]);
	  }
	return energy;
      }

    /**
     * Variable for finding pi0 leading shower start distance to vertex.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the leading shower distance to vertex   
     */
    template<class T>
      double pi0_leading_photon_start_to_vertex(const T & interaction)
      {
	
	// Interaction vertex
	TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);
	
	// Photon info
	vector<size_t> pi0_photon_idxs;
	size_t i;
	TVector3 sh_start;
        double s_to_v(0);

	// Truth
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         pi0_photon_idxs = true_pi0_photon_idxs(interaction);
                         i = pi0_photon_idxs[0];
                       }
	// Reco
	else
	  {
	    pi0_photon_idxs = pi0_photon_idxs_by_dir(interaction);
	    i = pi0_photon_idxs[0];
	  }

	sh_start.SetX(interaction.particles[i].start_point[0]);
	sh_start.SetY(interaction.particles[i].start_point[1]);
	sh_start.SetZ(interaction.particles[i].start_point[2]);

	s_to_v = (sh_start - vertex).Mag();
	return s_to_v;
      }

    /**
     * Variable for finding pi0 subleading shower start distance to vertex.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the leading shower distance to vertex
     */
    template<class T>
      double pi0_subleading_photon_start_to_vertex(const T & interaction)
      {

        // Interaction vertex                                                                                                                    
        TVector3 vertex(interaction.vertex[0], interaction.vertex[1], interaction.vertex[2]);

        // Photon info                                                                                        
        vector<size_t> pi0_photon_idxs;
        size_t i;
        TVector3 sh_start;
        double s_to_v(0);

        // Truth                                                                                                     
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                       {
                         pi0_photon_idxs = true_pi0_photon_idxs(interaction);
                         i = pi0_photon_idxs[1];
                       }
        // Reco                                                                                                           
        else
          {
            pi0_photon_idxs = pi0_photon_idxs_by_dir(interaction);
            i = pi0_photon_idxs[1];
          }

        sh_start.SetX(interaction.particles[i].start_point[0]);
        sh_start.SetY(interaction.particles[i].start_point[1]);
        sh_start.SetZ(interaction.particles[i].start_point[2]);

        s_to_v = (sh_start - vertex).Mag();
        return s_to_v;
      }
      
    /**
     * Find pi0 invariant mass.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to operate on.
     * @return the pi0 invariant mass.
     */
    
    template<class T>
      double pi0_mass(const T & interaction)
      {
	double ph0_energy;
	double ph1_energy;
	double costheta;
	double pi0_mass;

	ph0_energy = pi0_leading_photon_energy(interaction);
	ph1_energy = pi0_subleading_photon_energy(interaction);
	costheta = pi0_costheta(interaction);
	
	pi0_mass = sqrt(2*ph0_energy*ph1_energy*(1-costheta));

	return pi0_mass;
      }
    

    /**
     * Finds the index corresponding to the leading particle of the specifed
     * particle type.
     * @tparam T the type of intearction (true or reco).
     * @param interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
    */
    template <class T>
        size_t leading_particle_index(const T & interaction, uint16_t pid)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < interaction.particles.size(); ++i)
            {
                const auto & p = interaction.particles[i];
                double energy(csda_ke(p));
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = ke_init(p);
                if(p.pid == pid && energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }

    /**
     * Variable for particle overlap (IoU) of best match.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the particle overlap of the best match (IoU).
    */
    template<class T>
        double overlap(const T & particle) { return particle.match.size() > 0 ? (double)particle.match_overlap[0] : 0.0; }

    /**
     * Variable for the cosine of the track angle within the XZ plane.
     * @tparam T the type of particle (true or reco)
     * @param particle to apply the variable on.
     * @return the cosine of the track angle within the XZ plane ()
    */
    template<class T>
        double cosine_theta_xz(const T & particle)
        {
            return particle.start_dir[2] / std::sqrt(std::pow(particle.start_dir[0], 2) + std::pow(particle.start_dir[2], 2));
        }

    /**
     * Variable for cosine theta_xz (transverse) of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine theta_xz of the leading muon.
    */
    template<class T>
        double leading_muon_cosine_theta_xz(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return cosine_theta_xz(interaction.particles[i]);
        }








    

}

#endif
