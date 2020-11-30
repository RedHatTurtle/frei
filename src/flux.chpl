prototype module Flux
{
  proc pressure( const ref u : [1..3] real(64) ) : real(64)
  {
    use Input;

    pressure = (fGamma-1.0)*(u[3] - u[2]*u[2]/u[1]/2.0);
  }

  proc euler_flux_cv( const ref u : [1..3] real(64) ) : [1..3] real(64)
  {
    euler_flux_cv[1] = u[2];
    euler_flux_cv[2] = u[2]*u[2]/u[1] + pressure(u);
    euler_flux_cv[3] = u[3] + pressure(u);
  }

  proc euler_flux_pv( const ref p : [1..3] real(64) ) : [1..3] real(64)
  {
    use Input;

    euler_flux_pv[1] = p[1]*p[2];
    euler_flux_pv[2] = p[1]*p[2]*p[2] + p[3];
    euler_flux_pv[3] = p[3]/(fGamma-1.0) + p[1]*p[2]*p[2]/2.0;
  }
}
