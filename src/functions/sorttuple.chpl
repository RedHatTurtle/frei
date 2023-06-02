module SortTuple
{
  proc sort_tuple(in tuple : 4*int) : 4*int
  {
    if tuple[0] < tuple[2] then tuple[0] <=> tuple[2];
    if tuple[1] < tuple[3] then tuple[1] <=> tuple[3];
    if tuple[0] < tuple[1] then tuple[0] <=> tuple[1];
    if tuple[2] < tuple[3] then tuple[2] <=> tuple[3];
    if tuple[1] < tuple[2] then tuple[1] <=> tuple[2];
    return tuple;
  }
}
