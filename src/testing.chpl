module Testing
{
  proc error( in reference : real, in value : real) : real
  {
    return (value - reference);
  }

  proc absolute_error( in reference : real, in value : real) : real
  {
    return abs(value - reference);
  }

  proc relative_error( in reference : real, in value : real) : real
  {
    return (value - reference) / reference;
  }
}
