prototype module Spatial {
  proc residue {
    // Interpolate solution from SPs to FPs
    //interpolate_solution_from_sps_to_fps

    // Calculate solution at BCs
    //get_boundary_solution

  //if (diffusive flux) then {
      // Calculate discontinuous solution gradient
      //calculate_local_solution_gradient
      //calculate_discontinuous_solution_gradient

      // Calculate common solution gradient at faces
      //calculate_common_solution_gradient
      //calculate_interface_solution_gradient
      //calculate_continuous_solution_gradient

      // Calculate solution gradient at BCs
      //get_boundary_gradient
  //}

    // Calculate discontinuous flux residue component
    //calculate_local_flux_divergent
    //calculate_discontinuous_flux_divergent
      // Calculate fluxes at SPs, inviscid and viscous
      // Transform fluxes to the computational domain
      // Calculate the divergence of the discontinuous flux and add to residue
      // Interpolate discontinuous flux to FPs and project it in the face normal direction

    // Calculate continuous flux residue component
    //calculate_common_flux_divergent
    //calculate_interface_flux_divergent
    //calculate_continuous_flux_divergent
      // Calculate common fluxes
      // Calculate flux jump at interface (common flux minus local discontinuous flux)
      // Add corrections due to flux jump to residuals of the SPs

    // Return residue to RK and let it do it's shit
  }
}
