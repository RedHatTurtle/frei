/*
 * Copyright 2020-2023 Hewlett Packard Enterprise Development LP
 * Copyright 2004-2019 Cray Inc.
 * Other additional copyright holders may be indicated within.
 *
 * The entirety of this work is licensed under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License.
 *
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//
// The blockDist distribution is defined with six classes:
//
//   blockDist       : distribution class
//   blockDistDom    : domain class
//   blockDistArr    : array class
//   LocblockDist    : local distribution class (per-locale instances)
//   LocblockDistDom : local domain class (per-locale instances)
//   LocblockDistArr : local array class (per-locale instances)
//
// When a distribution, domain, or array class instance is created, a
// corresponding local class instance is created on each locale that is
// mapped to by the distribution.
//

//
// TO DO List
//
// 1. refactor pid fields from distribution, domain, and array classes
//

private use DSIUtil;
private use ChapelUtil;
private use CommDiagnostics;
private use ChapelLocks;
private use ChapelDebugPrint;
private use LayoutCS;

public use SparseblockDistDist;
//
// These flags are used to output debug information and run extra
// checks when using blockDist1D.  Should these be promoted so that they can
// be used across all distributions?  This can be done by turning them
// into compiler flags or adding config parameters to the internal
// modules, perhaps called debugDists and checkDists.
//
config param debugblockDistDist = false;
config param debugblockDistDistBulkTransfer = false;

// TODO: This is no longer used, deprecate with warning because it is used
// in miniMD's release Makefile and compopts.
config const disableAliasedBulkTransfer = true;

config param disableblockDistDistBulkTransfer = false;

config param sanityCheckDistribution = false;

private config param allowDuplicateTargetLocales = false;
//
// If the testFastFollowerOptimization flag is set to true, the
// follower will write output to indicate whether the fast follower is
// used or not.  This is used in regression testing to ensure that the
// 'fast follower' optimization is working.
//
config param testFastFollowerOptimization = false;

//
// This flag is used to disable lazy initialization of the RAD cache.
//
config param disableblockDistLazyRAD = defaultDisableLazyRADOpt;

//
// blockDist Distribution Class
//
//   The fields dataParTasksPerLocale, dataParIgnoreRunningTasks, and
//   dataParMinGranularity can be changed, but changes are
//   not reflected in privatized copies of this distribution.  Perhaps
//   this is a feature, not a bug?!
//
// rank : generic rank that must match the rank of domains and arrays
//        declared over this distribution
//
// idxType: generic index type that must match the index type of
//          domains and arrays declared over this distribution
//
// boundingBox: a non-distributed domain defining a bounding box used
//              to partition the space of all indices across the array
//              of target locales; the indices inside the bounding box
//              are partitioned "evenly" across the locales and
//              indices outside the bounding box are mapped to the
//              same locale as the nearest index inside the bounding
//              box
//
// targetLocDom: a non-distributed domain over which the array of
//               target locales and the array of local distribution
//               classes are defined
//
// targetLocales: a non-distributed array containing the target
//                locales to which this distribution partitions indices
//                and data
//
// locDist: a non-distributed array of local distribution classes
//
// dataParTasksPerLocale: an integer that specifies the number of tasks to
//                        use on each locale when iterating in parallel over
//                        a blockDist-distributed domain or array
//
// dataParIgnoreRunningTasks: a boolean what dictates whether the number of
//                            task use on each locale should be limited
//                            by the available parallelism
//
// dataParMinGranularity: the minimum required number of elements per
//                        task created
//

// chpldoc TODO:
//   a good reference to
//     dataParTasksPerLocale, dataParIgnoreRunningTasks, dataParMinGranularity
//   remove the above comments to avoid duplication/maintenance?
//   talk about RAD - here or in DefaultRectangular?
//   supports RAD opt, Bulk Transfer optimization, localSubdomain
//   disableblockDistLazyRAD
//
/*
This blockDist distribution partitions indices into blocks
according to a ``boundingBox`` domain
and maps each entire block onto a locale from a ``targetLocales`` array.

The indices inside the bounding box are partitioned "evenly" across
the target locales. An index outside the bounding box is mapped to the
same locale as the nearest index inside the bounding box.

Formally, an index ``idx`` is mapped to ``targetLocales[locIdx]``,
where ``locIdx`` is computed as follows.

In the 1-dimensional case, for a blockDist distribution with:


  =================   ====================
  ``boundingBox``     ``{low..high}``
  ``targetLocales``   ``[0..N-1] locale``
  =================   ====================

we have:

  ===================    ==========================================
  if ``idx`` is ...      ``locIdx`` is ...
  ===================    ==========================================
  ``low<=idx<=high``     ``floor(  (idx-low)*N / (high-low+1)  )``
  ``idx < low``          ``0``
  ``idx > high``         ``N-1``
  ===================    ==========================================

In the multidimensional case, ``idx`` and ``locIdx`` are tuples
of indices; ``boundingBox`` and ``targetLocales`` are multi-dimensional;
the above computation is applied in each dimension.


**Example**

The following code declares a domain ``D`` distributed over
a blockDist distribution with a bounding box equal to the domain ``Space``,
and declares an array ``A`` over that domain.
The `forall` loop sets each array element
to the ID of the locale to which it is mapped.

  .. code-block:: chapel

    use blockDistDist;

    const Space = {1..8, 1..8};
    const D: domain(2) dmapped blockDist1D(boundingBox=Space) = Space;
    var A: [D] int;

    forall a in A do
      a = a.locale.id;

    writeln(A);

When run on 6 locales, the output is:

  ::

    0 0 0 0 1 1 1 1
    0 0 0 0 1 1 1 1
    0 0 0 0 1 1 1 1
    2 2 2 2 3 3 3 3
    2 2 2 2 3 3 3 3
    2 2 2 2 3 3 3 3
    4 4 4 4 5 5 5 5
    4 4 4 4 5 5 5 5


**Initializer Arguments**

The ``blockDist`` class initializer is defined as follows:

  .. code-block:: chapel

    proc blockDist1D.init(
      boundingBox: domain,
      targetLocales: [] locale  = Locales,
      dataParTasksPerLocale     = // value of  dataParTasksPerLocale      config const,
      dataParIgnoreRunningTasks = // value of  dataParIgnoreRunningTasks  config const,
      dataParMinGranularity     = // value of  dataParMinGranularity      config const,
      param rank                = boundingBox.rank,
      type  idxType             = boundingBox.idxType,
      type  sparseLayoutType    = DefaultDist)

The arguments ``boundingBox`` (a domain) and ``targetLocales`` (an array)
define the mapping of any index of ``idxType`` type to a locale
as described above.

The rank of ``targetLocales`` must match the rank of the distribution,
or be ``1``.  If the rank of ``targetLocales`` is ``1``, a greedy
heuristic is used to reshape the array of target locales so that it
matches the rank of the distribution and each dimension contains an
approximately equal number of indices.

The arguments ``dataParTasksPerLocale``, ``dataParIgnoreRunningTasks``,
and ``dataParMinGranularity`` set the knobs that are used to
control intra-locale data parallelism for blockDist-distributed domains
and arrays in the same way that the like-named config constants
control data parallelism for ranges and default-distributed domains
and arrays.

The ``rank`` and ``idxType`` arguments are inferred from the ``boundingBox``
argument unless explicitly set.  They must match the rank and index type of the
domains "dmapped" using that blockDist instance. If the ``boundingBox`` argument is
a stridable domain, the stride information will be ignored and the
``boundingBox`` will only use the lo..hi bounds.

When a ``sparse subdomain`` is created for a ``blockDist`` distributed domain, the
``sparseLayoutType`` will be the layout of these sparse domains. The default is
currently coordinate, but :class:`LayoutCS.CS` is an interesting alternative.

**Convenience Initializer Functions**

It is common for a ``blockDist`` distribution to distribute its ``boundingBox``
across all locales. In this case, a convenience function can be used to
declare variables of block-distributed domain or array type.  These functions
take a domain or list of ranges as arguments and return a block-distributed
domain or array.

  .. code-block:: chapel

    use blockDistDist;

    var blockDistDom1 = blockDist1D.createDomain({1..5, 1..5});
    var blockDistArr1 = blockDist1D.createArray({1..5, 1..5}, real);
    var blockDistDom2 = blockDist1D.createDomain(1..5, 1..5);
    var blockDistArr2 = blockDist1D.createArray(1..5, 1..5, real);

**Data-Parallel Iteration**

A `forall` loop over a blockDist-distributed domain or array
executes each iteration on the locale where that iteration's index
is mapped to.

Parallelism within each locale is guided by the values of
``dataParTasksPerLocale``, ``dataParIgnoreRunningTasks``, and
``dataParMinGranularity`` of the respective blockDist instance.
Updates to these values, if any, take effect only on the locale
where the updates are made.

**Sparse Subdomains**

When a ``sparse subdomain`` is declared as a subdomain to a blockDist-distributed
domain, the resulting sparse domain will also be blockDist-distributed. The
sparse layout used in this sparse subdomain can be controlled with the
``sparseLayoutType`` initializer argument to blockDist1D.

This example demonstrates a blockDist-distributed sparse domain and array:

  .. code-block:: chapel

    use blockDistDist;

    const Space = {1..8, 1..8};

    // Declare a dense, blockDist-distributed domain.
    const DenseDom: domain(2) dmapped blockDist1D(boundingBox=Space) = Space;

    // Declare a sparse subdomain.
    // Since DenseDom is blockDist-distributed, SparseDom will be as well.
    var SparseDom: sparse subdomain(DenseDom);

    // Add some elements to the sparse subdomain.
    // SparseDom.bulkAdd is another way to do this that allows more control.
    SparseDom += [ (1,2), (3,6), (5,4), (7,8) ];

    // Declare a sparse array.
    // This array is also blockDist-distributed.
    var A: [SparseDom] int;

    A = 1;

    writeln( "A[(1, 1)] = ", A[1,1]);
    for (ij,x) in zip(SparseDom, A) {
      writeln( "A[", ij, "] = ", x, " on locale ", x.locale);
    }

   // Results in this output when run on 4 locales:
   // A[(1, 1)] = 0
   // A[(1, 2)] = 1 on locale LOCALE0
   // A[(3, 6)] = 1 on locale LOCALE1
   // A[(5, 4)] = 1 on locale LOCALE2
   // A[(7, 8)] = 1 on locale LOCALE3


*/
class blockDist1D : BaseDist {
  param rank: int;
  type idxType = int;
  var boundingBox: domain(rank, idxType);
  var targetLocDom: domain(rank);
  var targetLocales: [targetLocDom] locale;
  var locDist: [targetLocDom] unmanaged LocblockDist(rank, idxType);
  var dataParTasksPerLocale: int;
  var dataParIgnoreRunningTasks: bool;
  var dataParMinGranularity: int;
  type sparseLayoutType = unmanaged DefaultDist;
}

//
// Local blockDist Distribution Class
//
// rank : generic rank that matches blockDist1D.rank
// idxType: generic index type that matches blockDist1D.idxType
// myChunk: a non-distributed domain that defines this locale's indices
//
class LocblockDist {
  param rank: int;
  type idxType;
  var myChunk: domain(rank, idxType);
}

//
// blockDist Domain Class
//
// rank:      generic domain rank
// idxType:   generic domain index type
// strides:   generic domain stridable parameter
// dist:      reference to distribution class
// locDoms:   a non-distributed array of local domain classes
// whole:     a non-distributed domain that defines the domain's indices
//
class blockDistDom: BaseRectangularDom {
  type sparseLayoutType;
  const dist: unmanaged blockDist1D(rank, idxType, sparseLayoutType);
  var locDoms: [dist.targetLocDom] unmanaged LocblockDistDom(rank, idxType, strides);
  var whole: domain(rank, idxType, strides);
}

//
// Local blockDist Domain Class
//
// rank: generic domain rank
// idxType: generic domain index type
// strides: generic domain stridable parameter
// myblockDist: a non-distributed domain that defines the local indices
//
class LocblockDistDom {
  param rank: int;
  type idxType;
  param strides: strideKind;
  var myblockDist: domain(rank, idxType, strides);
}

//
// blockDist Array Class
//
// eltType: generic array element type
// rank: generic array rank
// idxType: generic array index type
// strides: generic array stridable parameter
// dom: reference to domain class
// locArr: a non-distributed array of local array classes
// myLocArr: optimized reference to here's local array class (or nil)
//
class blockDistArr: BaseRectangularArr {
  type sparseLayoutType;
  var doRADOpt: bool = defaultDoRADOpt;
  var dom: unmanaged blockDistDom(rank, idxType, strides, sparseLayoutType);
  var locArr: [dom.dist.targetLocDom] unmanaged LocblockDistArr(eltType,
                                                  rank, idxType, strides);
  pragma "local field"
  var myLocArr: unmanaged LocblockDistArr(eltType, rank, idxType, strides)?;
  const SENTINEL = max(rank*int);
}

//
// Local blockDist Array Class
//
// eltType: generic array element type
// rank: generic array rank
// idxType: generic array index type
// strides: generic array stridable parameter
// locDom: reference to local domain class
// myElems: a non-distributed array of local elements
//
class LocblockDistArr {
  type eltType;
  param rank: int;
  type idxType;
  param strides: strideKind;
  const locDom: unmanaged LocblockDistDom(rank, idxType, strides);
  var locRAD: unmanaged LocRADCache(eltType, rank, idxType, strides)?; // non-nil if doRADOpt=true
  pragma "local field" pragma "unsafe"
  // may be initialized separately
  var myElems: [locDom.myblockDist] eltType;
  var locRADLock: chpl_LocalSpinlock;

  proc init(type eltType,
            param rank: int,
            type idxType,
            param strides: strideKind,
            const locDom: unmanaged LocblockDistDom(rank, idxType, strides),
            param initElts: bool) {
    this.eltType = eltType;
    this.rank = rank;
    this.idxType = idxType;
    this.strides = strides;
    this.locDom = locDom;
    this.myElems = this.locDom.myblockDist.buildArray(eltType, initElts=initElts);
  }

  // guard against dynamic dispatch resolution trying to resolve
  // write()ing out an array of sync vars and hitting the sync var
  // type's compilerError()
  override proc writeThis(f) throws {
    halt("LocblockDistArr.writeThis() is not implemented / should not be needed");
  }

  proc deinit() {
    // Elements in myElems are deinited in dsiDestroyArr if necessary.
    // Here we need to clean up the rest of the array.
    if locRAD != nil then
      delete locRAD;
  }
}


////// blockDist and LocblockDist methods ///////////////////////////////////////////

//
// blockDist initializer for clients of the blockDist distribution
//
proc blockDist1D.init(boundingBox: domain,
                  targetLocales: [] locale = Locales,
                  dataParTasksPerLocale=getDataParTasksPerLocale(),
                  dataParIgnoreRunningTasks=getDataParIgnoreRunningTasks(),
                  dataParMinGranularity=getDataParMinGranularity(),
                  param rank = boundingBox.rank,
                  type idxType = boundingBox.idxType,
                  type sparseLayoutType = unmanaged DefaultDist) {
  this.rank = rank;
  this.idxType = idxType;
  if rank != boundingBox.rank then
    compilerError("specified blockDist rank != rank of specified bounding box");
  if idxType != boundingBox.idxType then
    compilerError("specified blockDist index type != index type of specified bounding box");
  if rank != 2 && isCSType(sparseLayoutType) then
    compilerError("CS layout is only supported for 2 dimensional domains");

  if boundingBox.sizeAs(uint) == 0 then
    halt("blockDist1D() requires a non-empty boundingBox");

  this.boundingBox = boundingBox : domain(rank, idxType);

  if !allowDuplicateTargetLocales {
    var checkArr: [LocaleSpace] bool;
    for loc in targetLocales {
      if checkArr[loc.id] then
        halt("blockDistDist does not allow duplicate targetLocales");
      checkArr[loc.id] = true;
    }
  }

  const ranges = setupTargetLocRanges(rank, targetLocales);
  this.targetLocDom = {(...ranges)};
  this.targetLocales = reshape(targetLocales, this.targetLocDom);

  // Instead of 'dummyLB', we could give 'locDistTemp' a nilable element type.
  const dummyLB = new unmanaged LocblockDist(rank, idxType, dummy=true);
  var locDistTemp: [targetLocDom] unmanaged LocblockDist(rank, idxType) = dummyLB;

  // Store a const copy of the dims for RVF
  // Use a zippered coforall to get nonblocking ons but create
  // the reference into locDistTemp before launching the remote task.
  const boundingBoxDims = boundingBox.dims();
  coforall (locid, loc, locDistTempElt)
           in zip(targetLocDom, targetLocales, locDistTemp) {
    on loc {
      locDistTempElt = new unmanaged LocblockDist(rank, idxType, locid,
                                              boundingBox, boundingBoxDims,
                                              targetLocDom);
    }
  }

  delete dummyLB;
  this.locDist = locDistTemp; //make this a serial loop instead?

  // NOTE: When these knobs stop using the global defaults, we will need
  // to add checks to make sure dataParTasksPerLocale<0 and
  // dataParMinGranularity<0
  this.dataParTasksPerLocale = if dataParTasksPerLocale==0
                               then here.maxTaskPar
                               else dataParTasksPerLocale;
  this.dataParIgnoreRunningTasks = dataParIgnoreRunningTasks;
  this.dataParMinGranularity = dataParMinGranularity;

  this.sparseLayoutType = _to_unmanaged(sparseLayoutType);

  this.complete();

  if debugblockDistDist {
    writeln("Creating new blockDist distribution:");
    dsiDisplayRepresentation();
  }
}

@unstable(category="experimental", reason="'blockDist1D.redistribute()' is currently unstable due to lack of design review and is being made available as a prototype")
proc blockDist1D.redistribute(const in newBbox) {
  const newBboxDims = newBbox.dims();
  const pid = this.pid;
  coforall (locid, loc, locdist) in zip(targetLocDom, targetLocales, locDist) {
    on loc {
      const that = if _privatization then chpl_getPrivatizedCopy(this.type, pid) else this;
      that.boundingBox = newBbox;

      var inds = chpl__computeblockDist(chpl__tuplify(locid), targetLocDom, newBbox, newBboxDims);
      locdist.myChunk = {(...inds)};
    }
  }
}


proc blockDist1D.dsiAssign(other: this.type) {
  if (this.targetLocDom != other.targetLocDom ||
      || reduce (this.targetLocales != other.targetLocales)) {
    halt("blockDist distribution assignments currently require the target locale arrays to match");
  }

  this.redistribute(other.boundingBox);
}

//
// blockDist distributions are equivalent if they share the same bounding
// box and target locale set.
//
proc blockDist1D.dsiEqualDMaps(that: blockDist1D(?)) {
  return (this.rank == that.rank &&
          this.boundingBox == that.boundingBox &&
          this.targetLocales.equals(that.targetLocales));
}

//
// blockDist distributions are not equivalent to other domain maps.
//
proc blockDist1D.dsiEqualDMaps(that) param {
  return false;
}

proc blockDist1D.dsiClone() {
  return new unmanaged blockDist1D(boundingBox, targetLocales,
                   dataParTasksPerLocale, dataParIgnoreRunningTasks,
                   dataParMinGranularity,
                   rank,
                   idxType,
                   sparseLayoutType);
}

override proc blockDist1D.dsiDestroyDist() {
  coforall ld in locDist do {
    on ld do
      delete ld;
  }
}

override proc blockDist1D.dsiDisplayRepresentation() {
  writeln("boundingBox = ", boundingBox);
  writeln("targetLocDom = ", targetLocDom);
  writeln("targetLocales = ", for tl in targetLocales do tl.id);
  writeln("dataParTasksPerLocale = ", dataParTasksPerLocale);
  writeln("dataParIgnoreRunningTasks = ", dataParIgnoreRunningTasks);
  writeln("dataParMinGranularity = ", dataParMinGranularity);
  for tli in targetLocDom do
    writeln("locDist[", tli, "].myChunk = ", locDist[tli].myChunk);
}

override proc blockDist1D.dsiNewRectangularDom(param rank: int, type idxType,
                                         param strides: strideKind, inds) {
  if idxType != this.idxType then
    compilerError("blockDist domain index type does not match distribution's");
  if rank != this.rank then
    compilerError("blockDist domain rank does not match distribution's");

  const whole = createWholeDomainForInds(rank, idxType, strides, inds);

  const dummyLBD = new unmanaged LocblockDistDom(rank, idxType, strides);
  var locDomsTemp: [this.targetLocDom]
                  unmanaged LocblockDistDom(rank, idxType, strides) = dummyLBD;
  coforall (localeIdx, loc, locDomsTempElt)
           in zip(this.targetLocDom, this.targetLocales, locDomsTemp) {
    on loc {
      locDomsTempElt = new unmanaged LocblockDistDom(rank, idxType, strides,
        // todo: can/should we get a more specific return type out of getChunk?
        this.getChunk(whole, localeIdx): domain(rank, idxType, strides));
    }
  }
  delete dummyLBD;

  var dom = new unmanaged blockDistDom(rank, idxType, strides, sparseLayoutType,
                                   _to_unmanaged(this), locDomsTemp, whole);

  if debugblockDistDist {
    writeln("Creating new blockDist domain:");
    dom.dsiDisplayRepresentation();
  }
  return dom;
}

override proc blockDist1D.dsiNewSparseDom(param rank: int, type idxType,
                                    dom: domain) {
  var ret =  new unmanaged SparseblockDistDom(rank=rank, idxType=idxType,
                            sparseLayoutType=sparseLayoutType,
                            strides=dom.strides,
                            dist=_to_unmanaged(this), whole=dom._value.whole,
                            parentDom=dom);
  ret.setup();
  return ret;
}

//
// output distribution
//
proc blockDist1D.writeThis(x) throws {
  x.writeln("blockDist");
  x.writeln("-------");
  x.writeln("distributes: ", boundingBox);
  x.writeln("across locales: ", targetLocales);
  x.writeln("indexed via: ", targetLocDom);
  x.writeln("resulting in: ");
  for locid in targetLocDom do
    x.writeln("  [", locid, "] locale ", locDist(locid).locale.id,
      " owns chunk: ", locDist(locid).myChunk);
}

proc blockDist1D.dsiIndexToLocale(ind: idxType) where rank == 1 {
  return targetLocales(targetLocsIdx(ind));
}

proc blockDist1D.dsiIndexToLocale(ind: rank*idxType) {
  return targetLocales(targetLocsIdx(ind));
}

//
// compute what chunk of inds is owned by a given locale -- assumes
// it's being called on the locale in question
//
proc blockDist1D.getChunk(inds, locid) {
  // use domain slicing to get the intersection between what the
  // locale owns and the domain's index set
  //
  // TODO: Should this be able to be written as myChunk[inds] ???
  //
  // TODO: Does using David's detupling trick work here?
  //
  // Vass 2023-03: the chunk should really be computed as:
  //   const chunk = inds[locDist(locid).myChunk];
  // because we are looking for a subset of 'inds'.
  // I did not make this change because it would slightly bump comm counts in:
  //   distributions/robust/arithmetic/performance/multilocale/assignReindex
  //
  const chunk = locDist(locid).myChunk((...inds.getIndices()));
  if sanityCheckDistribution then
    if chunk.sizeAs(int) > 0 {
      if targetLocsIdx(chunk.lowBound) != locid then
        writeln("[", here.id, "] ", chunk.lowBound, " is in my chunk but maps to ",
                targetLocsIdx(chunk.lowBound));
      if targetLocsIdx(chunk.highBound) != locid then
        writeln("[", here.id, "] ", chunk.highBound, " is in my chunk but maps to ",
                targetLocsIdx(chunk.highBound));
    }
  return chunk;
}

//
// get the index into the targetLocales array for a given distributed index
//
proc blockDist1D.targetLocsIdx(ind: idxType) where rank == 1 {
  return targetLocsIdx((ind,));
}

proc blockDist1D.targetLocsIdx(ind: rank*idxType) {
  var result: rank*int;
  for param i in 0..rank-1 do
    result(i) = max(0, min(targetLocDom.dim(i).sizeAs(int)-1,
                           (((ind(i) - boundingBox.dim(i).lowBound).safeCast(int) *
                             targetLocDom.dim(i).sizeAs(int)) /
                            boundingBox.dim(i).sizeAs(int))));
  return if rank == 1 then result(0) else result;
}

// TODO: This will not trigger the bounded-coforall optimization
iter blockDist1D.activeTargetLocales(const space : domain = boundingBox) {
  const locSpace = {(...space.dims())}; // make a local domain in case 'space' is distributed
  const low  = chpl__tuplify(targetLocsIdx(locSpace.low));
  const high = chpl__tuplify(targetLocsIdx(locSpace.high));
  var dims : rank*range(low(0).type);
  for param i in 0..rank-1 {
    dims(i) = low(i)..high(i);
  }

  // In case 'locSpace' is a strided domain we need to check that the locales
  // in 'dims' actually contain indices in 'locSpace'.
  //
  // Note that we cannot use a simple stride here because it is not guaranteed
  // that each locale contains the same number of indices. For example, the
  // domain {1..10} over four locales will split like:
  //   L0: -max(int)..3
  //   L1: 4..5
  //   L2: 6..8
  //   L3: 9..max(int)
  //
  // The subset {1..10 by 4} will involve locales 0, 1, and 3.
  foreach i in {(...dims)} {
    const chunk = chpl__computeblockDist(i, targetLocDom, boundingBox, boundingBox.dims());
    // TODO: Want 'contains' for a domain. Slicing is a workaround.
    if locSpace[(...chunk)].sizeAs(int) > 0 then
      yield i;
  }
}

proc type blockDist1D.createDomain(dom: domain) {
  return dom dmapped blockDist1D(dom);
}

proc type blockDist1D.createDomain(rng: range...) {
  return createDomain({(...rng)});
}

proc type blockDist1D.createArray(dom: domain, type eltType) {
  var D = createDomain(dom);
  var A: [D] eltType;
  return A;
}

proc type blockDist1D.createArray(rng: range..., type eltType) {
  return createArray({(...rng)}, eltType);
}

proc chpl__computeblockDist(locid, targetLocBox:domain, boundingBox:domain,
                        boundingBoxDims /* boundingBox.dims() */) {
  param rank = targetLocBox.rank;
  type idxType = boundingBox.idxType;
  var inds: rank*range(idxType);
  for param i in 0..rank-1 {
    const lo = boundingBoxDims(i).lowBound;
    const hi = boundingBoxDims(i).highBound;
    const numelems = hi - lo + 1;
    const numlocs = targetLocBox.dim(i).sizeAs(int);
    const (blo, bhi) = _computeblockDist(numelems, numlocs, chpl__tuplify(locid)(i),
                                     max(idxType), min(idxType), lo);
    inds(i) = blo..bhi;
  }
  return inds;
}

proc LocblockDist.init(param rank: int,
                   type idxType,
                   locid, // the locale index from the target domain
                   boundingBox: domain,
                   boundingBoxDims /* boundingBox.dims() */,
                   targetLocDom: domain(rank)) {
  this.rank = rank;
  this.idxType = idxType;
  const inds = chpl__computeblockDist(chpl__tuplify(locid), targetLocDom, boundingBox, boundingBoxDims);
  myChunk = {(...inds)};
}

// Used to create a dummy instance.
proc LocblockDist.init(param rank, type idxType, param dummy: bool) where dummy {
  this.rank = rank;
  this.idxType = idxType;
}


////// blockDistDom and LocblockDistDom methods /////////////////////////////////////

override proc blockDistDom.dsiDisplayRepresentation() {
  writeln("whole = ", whole);
  for tli in dist.targetLocDom do
    writeln("locDoms[", tli, "].myblockDist = ", locDoms[tli].myblockDist);
}

// stopgap to avoid accessing locDoms field (and returning an array)
proc blockDistDom.getLocDom(localeIdx) do return locDoms(localeIdx);

//
// Given a tuple of scalars of type t or range(t) match the shape but
// using types rangeType and scalarType e.g. the call:
// _matchArgsShape(range(int(32)), int(32), (1:int(64), 1:int(64)..5, 1:int(64)..5))
// returns the type: (int(32), range(int(32)), range(int(32)))
//
proc _matchArgsShape(type rangeType, type scalarType, args) type {
  proc helper(param i: int) type {
    if i == args.size {
      if isCollapsedDimension(args(i)) then
        return (scalarType,);
      else
        return (rangeType,);
    } else {
      if isCollapsedDimension(args(i)) then
        return (scalarType, (... helper(i+1)));
      else
        return (rangeType, (... helper(i+1)));
    }
  }
  return helper(0);
}


iter blockDistDom.these() {
  for i in whole do
    yield i;
}

iter blockDistDom.these(param tag: iterKind) where tag == iterKind.leader {
  const maxTasks = dist.dataParTasksPerLocale;
  const ignoreRunning = dist.dataParIgnoreRunningTasks;
  const minSize = dist.dataParMinGranularity;
  const wholeLow = whole.lowBound;

  // If this is the only task running on this locale, we don't want to
  // count it when we try to determine how many tasks to use.  Here we
  // check if we are the only one running, and if so, use
  // ignoreRunning=true for this locale only.  Obviously there's a bit
  // of a race condition if some other task starts after we check, but
  // in that case there is no correct answer anyways.
  //
  // Note that this code assumes that any locale will only be in the
  // targetLocales array once.  If this is not the case, then the
  // tasks on this locale will *all* ignoreRunning, which may have
  // performance implications.
  const hereId = here.id;
  const hereIgnoreRunning = if here.runningTasks() == 1 then true
                            else ignoreRunning;
  coforall locDom in locDoms do on locDom {
    const myIgnoreRunning = if here.id == hereId then hereIgnoreRunning
      else ignoreRunning;
    // Use the internal function for untranslate to avoid having to do
    // extra work to negate the offset
    const tmpblockDist = locDom.myblockDist.chpl__unTranslate(wholeLow);
    var locOffset: rank*idxType;
    for param i in 0..tmpblockDist.rank-1 {
      const dim = tmpblockDist.dim(i);
      const aStr = if dim.hasPositiveStride() then dim.stride else -dim.stride;
      locOffset(i) = dim.low / aStr:idxType;
    }
    // Forward to defaultRectangular
    for followThis in tmpblockDist.these(iterKind.leader, maxTasks,
                                     myIgnoreRunning, minSize, locOffset) do
      yield followThis;
  }
}

//
// TODO: Abstract the addition of low into a function?
// Note relationship between this operation and the
// order/position functions -- any chance for creating similar
// support? (esp. given how frequent this seems likely to be?)
//
// TODO: Is there some clever way to invoke the leader/follower
// iterator on the local blocks in here such that the per-core
// parallelism is expressed at that level?  Seems like a nice
// natural composition and might help with my fears about how
// stencil communication will be done on a per-locale basis.
//
// TODO: rewrite the index transformations in this and similar
// leaders+followers to use un/densify, which were created for that.
// If un/densify add overhead, need to eliminate it.
//
// TODO: Can we just re-use the DefaultRectangularDom follower here?
//
iter blockDistDom.these(param tag: iterKind, followThis) where tag == iterKind.follower {
  if chpl__testParFlag then
    chpl__testParWriteln("blockDist domain follower invoked on ", followThis);

  var t: rank*range(idxType, strides = chpl_strideProduct(whole.strides,
                                              chpl_strideUnion(followThis)));
  for param i in 0..rank-1 {
    const wholeDim  = whole.dim(i);
    const followDim = followThis(i);
    var low  = wholeDim.orderToIndex(followDim.low);
    var high = wholeDim.orderToIndex(followDim.high);
    if wholeDim.hasNegativeStride() then low <=> high;
    t(i) = ( low..high by (wholeDim.stride*followDim.stride)
           ).safeCast(t(i).type);
  }
  for i in {(...t)} {
    yield i;
  }
}

//
// how to allocate a new array over this domain
//
proc blockDistDom.dsiBuildArray(type eltType, param initElts:bool) {
  const dom = this;
  const creationLocale = here.id;
  const dummyLBD = new unmanaged LocblockDistDom(rank, idxType, strides);
  const dummyLBA = new unmanaged LocblockDistArr(eltType, rank, idxType,
                                             strides, dummyLBD, false);
  var locArrTemp: [dom.dist.targetLocDom]
        unmanaged LocblockDistArr(eltType, rank, idxType, strides) = dummyLBA;
  var myLocArrTemp: unmanaged LocblockDistArr(eltType, rank, idxType, strides)?;

  // formerly in blockDistArr.setup()
  coforall (loc, locDomsElt, locArrTempElt)
           in zip(dom.dist.targetLocales, dom.locDoms, locArrTemp)
           with (ref myLocArrTemp) {
    on loc {
      const LBA = new unmanaged LocblockDistArr(eltType, rank, idxType, strides,
                                            locDomsElt,
                                            initElts=initElts);
      locArrTempElt = LBA;
      if here.id == creationLocale then
        myLocArrTemp = LBA;
    }
  }
  delete dummyLBA, dummyLBD;

  var arr = new unmanaged blockDistArr(eltType=eltType, rank=rank, idxType=idxType,
       strides=strides, sparseLayoutType=sparseLayoutType,
       dom=_to_unmanaged(dom), locArr=locArrTemp, myLocArr=myLocArrTemp);

  // formerly in blockDistArr.setup()
  if arr.doRADOpt && disableblockDistLazyRAD then arr.setupRADOpt();

  return arr;
}

// common redirects
proc blockDistDom.parSafe param {
  compilerError("this domain type does not support 'parSafe'");
}
override proc blockDistDom.dsiLow do           return whole.lowBound;
override proc blockDistDom.dsiHigh do          return whole.highBound;
override proc blockDistDom.dsiAlignedLow do    return whole.low;
override proc blockDistDom.dsiAlignedHigh do   return whole.high;
override proc blockDistDom.dsiFirst do         return whole.first;
override proc blockDistDom.dsiLast do          return whole.last;
override proc blockDistDom.dsiStride do        return whole.stride;
override proc blockDistDom.dsiAlignment do     return whole.alignment;
proc blockDistDom.dsiNumIndices do    return whole.sizeAs(uint);
proc blockDistDom.dsiDim(d) do        return whole.dim(d);
proc blockDistDom.dsiDim(param d) do  return whole.dim(d);
proc blockDistDom.dsiDims() do        return whole.dims();
proc blockDistDom.dsiGetIndices() do  return whole.getIndices();
proc blockDistDom.dsiMember(i) do     return whole.contains(i);
proc blockDistDom.doiToString() do    return whole:string;
proc blockDistDom.dsiSerialWrite(x) { x.write(whole); }
proc blockDistDom.dsiLocalSlice(param strides, ranges) do return whole((...ranges));
override proc blockDistDom.dsiIndexOrder(i) do              return whole.indexOrder(i);
override proc blockDistDom.dsiMyDist() do                   return dist;

//
// INTERFACE NOTES: Could we make dsiSetIndices() for a rectangular
// domain take a domain rather than something else?
//
proc blockDistDom.dsiSetIndices(x: domain) {
  if x.rank != rank then
    compilerError("rank mismatch in domain assignment");
  if x._value.idxType != idxType then
    compilerError("index type mismatch in domain assignment");
  whole = x;
  setup();
  if debugblockDistDist {
    writeln("Setting indices of blockDist domain:");
    dsiDisplayRepresentation();
  }
}

proc blockDistDom.dsiSetIndices(x) {
  if x.size != rank then
    compilerError("rank mismatch in domain assignment");
  if x(0).idxType != idxType then
    compilerError("index type mismatch in domain assignment");
  //
  // TODO: This seems weird:
  //
  whole.setIndices(x);
  setup();
  if debugblockDistDist {
    writeln("Setting indices of blockDist domain:");
    dsiDisplayRepresentation();
  }
}

proc blockDistDom.dsiAssignDomain(rhs: domain, lhsPrivate:bool) {
  chpl_assignDomainWithGetSetIndices(this, rhs);
}

proc blockDistDom.setup() {
  coforall (localeIdx, loc, locDomsElt)
           in zip(dist.targetLocDom, dist.targetLocales, locDoms) {
    on loc {
      locDomsElt.myblockDist = dist.getChunk(whole, localeIdx);
    }
  }
}

override proc blockDistDom.dsiDestroyDom() {
  coforall (loc, locDomsElt) in zip(dist.targetLocales, locDoms) {
    on loc {
      delete locDomsElt;
    }
  }
}

//
// Added as a performance stopgap to avoid returning a domain
//
proc LocblockDistDom.contains(i) do return myblockDist.contains(i);


////// blockDistArr and LocblockDistArr methods /////////////////////////////////////

override proc blockDistArr.dsiDisplayRepresentation() {
  for tli in dom.dist.targetLocDom {
    writeln("locArr[", tli, "].myElems = ", for e in locArr[tli].myElems do e);
    if doRADOpt && locArr[tli].locRAD != nil then
      writeln("locArr[", tli, "].locRAD = ", locArr[tli].locRAD!.RAD);
  }
}

override proc blockDistArr.dsiGetBaseDom() do return dom;

override proc blockDistArr.dsiIteratorYieldsLocalElements() param {
  return true;
}

override proc blockDistDom.dsiIteratorYieldsLocalElements() param {
  return true;
}

//
// NOTE: Each locale's myElems array must be initialized prior to
// setting up the RAD cache.
//
proc blockDistArr.setupRADOpt() {
  for localeIdx in dom.dist.targetLocDom {
    on dom.dist.targetLocales(localeIdx) {
      const myLocArr = locArr(localeIdx);
      if myLocArr.locRAD != nil {
        delete myLocArr.locRAD;
        myLocArr.locRAD = nil;
      }
      if disableblockDistLazyRAD {
        myLocArr.locRAD = new unmanaged LocRADCache(eltType, rank, idxType,
                                             strides, dom.dist.targetLocDom);
        for l in dom.dist.targetLocDom {
          if l != localeIdx {
            myLocArr.locRAD!.RAD(l) = locArr(l).myElems._value.dsiGetRAD();
          }
        }
      }
    }
  }
}

override proc blockDistArr.dsiElementInitializationComplete() {
  coforall locArrElt in locArr {
    on locArrElt {
      locArrElt.myElems.dsiElementInitializationComplete();
    }
  }
}

override proc blockDistArr.dsiElementDeinitializationComplete() {
  coforall locArrElt in locArr {
    on locArrElt {
      locArrElt.myElems.dsiElementDeinitializationComplete();
    }
  }
}

override proc blockDistArr.dsiDestroyArr(deinitElts:bool) {
  coforall locArrElt in locArr {
    on locArrElt {
      var arr = locArrElt;
      if deinitElts then
        _deinitElements(arr.myElems);
      arr.myElems.dsiElementDeinitializationComplete();
      delete arr;
    }
  }
}

inline proc blockDistArr.dsiLocalAccess(i: rank*idxType) ref {
  return if allowDuplicateTargetLocales then this.dsiAccess(i)
                                        else _to_nonnil(myLocArr).this(i);
}

//
// the global accessor for the array
//
// TODO: Do we need a global bounds check here or in targetLocsIdx?
//
// By splitting the non-local case into its own function, we can inline the
// fast/local path and get better performance.
//
inline proc blockDistArr.dsiAccess(const in idx: rank*idxType) ref {
  local {
    if const myLocArrNN = myLocArr then
      if myLocArrNN.locDom.contains(idx) then
        return myLocArrNN.this(idx);
  }
  return nonLocalAccess(idx);
}

inline proc blockDistArr.dsiBoundsCheck(i: rank*idxType) {
  return dom.dsiMember(i);
}

pragma "fn unordered safe"
proc blockDistArr.nonLocalAccess(i: rank*idxType) ref {
  if doRADOpt {
    if const myLocArr = this.myLocArr {
      var rlocIdx = dom.dist.targetLocsIdx(i);
      if !disableblockDistLazyRAD {
        if myLocArr.locRAD == nil {
          myLocArr.locRADLock.lock();
          if myLocArr.locRAD == nil {
            var tempLocRAD = new unmanaged LocRADCache(eltType, rank, idxType,
                                               strides, dom.dist.targetLocDom);
            tempLocRAD.RAD.blk = SENTINEL;
            myLocArr.locRAD = tempLocRAD;
          }
          myLocArr.locRADLock.unlock();
        }
        const locRAD = _to_nonnil(myLocArr.locRAD);
        if locRAD.RAD(rlocIdx).blk == SENTINEL {
          locRAD.lockRAD(rlocIdx);
          if locRAD.RAD(rlocIdx).blk == SENTINEL {
            locRAD.RAD(rlocIdx) =
              locArr(rlocIdx).myElems._value.dsiGetRAD();
          }
          locRAD.unlockRAD(rlocIdx);
        }
      }
      pragma "no copy" pragma "no auto destroy" var myLocRAD = myLocArr.locRAD;
      pragma "no copy" pragma "no auto destroy" var radata = _to_nonnil(myLocRAD).RAD;
      if radata(rlocIdx).shiftedData != nil {
        var dataIdx = radata(rlocIdx).getDataIndex(i);
        return radata(rlocIdx).getDataElem(dataIdx);
      }
    }
  }
  return locArr(dom.dist.targetLocsIdx(i))(i);
}

proc blockDistArr.dsiAccess(i: idxType...rank) ref do
  return dsiAccess(i);

iter blockDistArr.these() ref {
  foreach i in dom do
    yield dsiAccess(i);
}

//
// TODO: Rewrite this to reuse more of the global domain iterator
// logic?  (e.g., can we forward the forall to the global domain
// somehow?
//
iter blockDistArr.these(param tag: iterKind) where tag == iterKind.leader {
  for followThis in dom.these(tag) do
    yield followThis;
}

override proc blockDistArr.dsiStaticFastFollowCheck(type leadType) param {
  if isSubtype(leadType, blockDistArr) {
    // Comparing domains for rank/idxType/stride/sparseLayout, which allows for
    // fast followers regardless of eltType.
    //
    // TODO: Remove once 'typeExpr.field' results in a type
    var x : leadType?;
    return _to_borrowed(x!.dom.type) == _to_borrowed(this.dom.type);
  } else {
    return _to_borrowed(leadType) == _to_borrowed(this.dom.type);
  }
}

proc blockDistArr.dsiDynamicFastFollowCheck(lead: []) do
  return this.dsiDynamicFastFollowCheck(lead.domain);

proc blockDistArr.dsiDynamicFastFollowCheck(lead: domain) {
  // TODO: Should this return true for domains with the same shape?
  return lead.dist.dsiEqualDMaps(this.dom.dist) && lead._value.whole == this.dom.whole;
}

iter blockDistArr.these(param tag: iterKind, followThis, param fast: bool = false) ref where tag == iterKind.follower {
  if chpl__testParFlag {
    if fast then
      chpl__testParWriteln("blockDist array fast follower invoked on ", followThis);
    else
      chpl__testParWriteln("blockDist array non-fast follower invoked on ", followThis);
  }

  if testFastFollowerOptimization then
    writeln((if fast then "fast" else "regular") + " follower invoked for blockDist array");

  var myFollowThis: rank*range(idxType=idxType, strides=chpl_strideProduct(
                                     strides, chpl_strideUnion(followThis)));
  var lowIdx: rank*idxType;

  for param i in 0..rank-1 {
    const stride = dom.whole.dim(i).stride;
    // NOTE: Not bothering to check to see if these can fit into idxType
    var low = followThis(i).lowBound * abs(stride):idxType;
    var high = followThis(i).highBound * abs(stride):idxType;
    myFollowThis(i) = ((low..high by stride) + dom.whole.dim(i).low by followThis(i).stride).safeCast(myFollowThis(i).type);
    lowIdx(i) = myFollowThis(i).lowBound;
  }

  const myFollowThisDom = {(...myFollowThis)};
  if fast {
    //
    // TODO: The following is a buggy hack that will only work when we're
    // distributing across the entire Locales array.  I still think the
    // locArr/locDoms arrays should be associative over locale values.
    //
    var arrSection = locArr(dom.dist.targetLocsIdx(lowIdx));

    //
    // if arrSection is not local and we're using the fast follower,
    // it means that myFollowThisDom is empty; make arrSection local so
    // that we can use the local block below
    //
    if arrSection.locale.id != here.id then
      arrSection = _to_nonnil(myLocArr);

    local {
      use CTypes; // Needed to cast from c_void_ptr in the next line
      const narrowArrSection =
        __primitive("_wide_get_addr", arrSection):arrSection.type?;
      ref myElems = _to_nonnil(narrowArrSection).myElems;
      foreach i in myFollowThisDom do yield myElems[i];
    }
  } else {
    //
    // we don't necessarily own all the elements we're following
    //
    foreach i in myFollowThisDom {
      yield dsiAccess(i);
    }
  }
}

proc blockDistArr.dsiSerialRead(f) {
  chpl_serialReadWriteRectangular(f, this);
}

//
// output array
//
proc blockDistArr.dsiSerialWrite(f) {
  chpl_serialReadWriteRectangular(f, this);
}

pragma "no copy return"
proc blockDistArr.dsiLocalSlice(ranges) {
  var low: rank*idxType;
  for param i in 0..rank-1 {
    low(i) = ranges(i).low;
  }

  return locArr(dom.dist.targetLocsIdx(low)).myElems((...ranges));
}

proc _extendTuple(type t, idx: _tuple, args) {
  var tup: args.size*t;
  var j: int = 1;

  for param i in 0..args.size-1 {
    if isCollapsedDimension(args(i)) then
      tup(i) = args(i);
    else {
      tup(i) = idx(j);
      j += 1;
    }
  }
  return tup;
}

proc _extendTuple(type t, idx, args) {
  var tup: args.size*t;
  var idxTup = (idx,);
  var j: int = 1;

  for param i in 0..args.size-1 {
    if isCollapsedDimension(args(i)) then
      tup(i) = args(i);
    else {
      tup(i) = idxTup(j);
      j += 1;
    }
  }
  return tup;
}

override proc blockDistArr.dsiReallocate(bounds:rank*range(idxType,boundKind.both,strides))
{
  //
  // For the default rectangular array, this function changes the data
  // vector in the array class so that it is setup once the default
  // rectangular domain is changed.  For this distributed array class,
  // we don't need to do anything, because changing the domain will
  // change the domain in the local array class which will change the
  // data in the local array class.  This will only work if the domain
  // we are reallocating to has the same distribution, but domain
  // assignment is defined so that only the indices are transferred.
  // The distribution remains unchanged.
  //
}

override proc blockDistArr.dsiPostReallocate() {
  // Call this *after* the domain has been reallocated
  if doRADOpt then setupRADOpt();
}

proc blockDistArr.setRADOpt(val=true) {
  doRADOpt = val;
  if doRADOpt then setupRADOpt();
}

//
// the accessor for the local array -- assumes the index is local
//
// TODO: Should this be inlined?
//
inline proc LocblockDistArr.this(i) ref {
  return myElems(i);
}

override proc blockDistDom.dsiSupportsAutoLocalAccess() param {
  return true;
}

///// Privatization and serialization ///////////////////////////////////////

proc blockDist1D.init(other: blockDist1D, privateData,
                  param rank = other.rank,
                  type idxType = other.idxType,
                  type sparseLayoutType = other.sparseLayoutType) {
  this.rank = rank;
  this.idxType = idxType;
  boundingBox = {(...privateData(0))};
  targetLocDom = {(...privateData(1))};
  targetLocales = other.targetLocales;
  locDist = other.locDist;
  dataParTasksPerLocale = privateData(2);
  dataParIgnoreRunningTasks = privateData(3);
  dataParMinGranularity = privateData(4);
  this.sparseLayoutType = sparseLayoutType;
}

override proc blockDist1D.dsiSupportsPrivatization() param do return true;

proc blockDist1D.dsiGetPrivatizeData() {
  return (boundingBox.dims(), targetLocDom.dims(),
          dataParTasksPerLocale, dataParIgnoreRunningTasks, dataParMinGranularity);
}

proc blockDist1D.dsiPrivatize(privatizeData) {
  return new unmanaged blockDist1D(_to_unmanaged(this), privatizeData);
}

proc blockDist1D.dsiGetReprivatizeData() do return boundingBox.dims();

proc blockDist1D.dsiReprivatize(other, reprivatizeData) {
  boundingBox = {(...reprivatizeData)};
  targetLocDom = other.targetLocDom;
  targetLocales = other.targetLocales;
  locDist = other.locDist;
  dataParTasksPerLocale = other.dataParTasksPerLocale;
  dataParIgnoreRunningTasks = other.dataParIgnoreRunningTasks;
  dataParMinGranularity = other.dataParMinGranularity;
}

proc blockDistDom.chpl__serialize() {
  return pid;
}

// TODO: What happens when we try to deserialize on a locale that doesn't
// own a copy of the privatized class?  (can that happen?)  Could this
// be a way to lazily privatize by also making the originating locale part
// of the 'data'?
proc type blockDistDom.chpl__deserialize(data) {
  return chpl_getPrivatizedCopy(
           unmanaged blockDistDom(rank=this.rank,
                              idxType=this.idxType,
                              strides=this.strides,
                              sparseLayoutType=this.sparseLayoutType),
           data);
}

override proc blockDistDom.dsiSupportsPrivatization() param do return true;

record blockDistDomPrvData {
  var distpid;
  var dims;
  var locdoms;  //todo rvf its elements along with the rest of the record
}

proc blockDistDom.dsiGetPrivatizeData() {
  return new blockDistDomPrvData(dist.pid, whole.dims(), locDoms);
}

proc blockDistDom.dsiPrivatize(privatizeData) {
  var privdist = chpl_getPrivatizedCopy(dist.type, privatizeData.distpid);

  var locDomsTemp: [privdist.targetLocDom]
                      unmanaged LocblockDistDom(rank, idxType, strides)
    = privatizeData.locdoms;

  // in initializer we have to pass sparseLayoutType as it has no default value
  const c = new unmanaged blockDistDom(rank, idxType, strides,
                                   privdist.sparseLayoutType, privdist,
                                   locDomsTemp, {(...privatizeData.dims)});
  return c;
}

proc blockDistDom.dsiGetReprivatizeData() do return whole.dims();

proc blockDistDom.dsiReprivatize(other, reprivatizeData) {
  locDoms = other.locDoms;
  whole = {(...reprivatizeData)};
}

proc blockDistArr.chpl__serialize()
      where !(isDomainType(eltType) || isArrayType(eltType)) {
  return pid;
}

proc type blockDistArr.chpl__deserialize(data) {
  return chpl_getPrivatizedCopy(
           unmanaged blockDistArr(rank=this.rank,
                              idxType=this.idxType,
                              strides=this.strides,
                              eltType=this.eltType,
                              sparseLayoutType=this.sparseLayoutType),
           data);
}

override proc blockDistArr.dsiSupportsPrivatization() param do return true;

record blockDistArrPrvData {
  var dompid;
  var locarr;  //todo rvf its elements along with the rest of the record
}

proc blockDistArr.dsiGetPrivatizeData() {
  return new blockDistArrPrvData(dom.pid, locArr);
}

proc blockDistArr.dsiPrivatize(privatizeData) {
  var privdom = chpl_getPrivatizedCopy(dom.type, privatizeData.dompid);

  var locArrTemp: [privdom.dist.targetLocDom]
                     unmanaged LocblockDistArr(eltType, rank, idxType, strides)
    = privatizeData.locarr;

  var myLocArrTemp: unmanaged LocblockDistArr(eltType, rank, idxType, strides)?;
  for localeIdx in privdom.dist.targetLocDom do
    if locArrTemp(localeIdx).locale.id == here.id then
      myLocArrTemp = locArrTemp(localeIdx);

  const c = new unmanaged blockDistArr(eltType=eltType, rank=rank, idxType=idxType,
                      strides=strides, sparseLayoutType=sparseLayoutType,
                      dom=privdom, locArr=locArrTemp, myLocArr = myLocArrTemp);
  return c;
}


////// more /////////////////////////////////////////////////////////////////

proc blockDistArr.dsiTargetLocales() const ref {
  return dom.dist.targetLocales;
}

proc blockDistDom.dsiTargetLocales() const ref {
  return dist.targetLocales;
}

proc blockDist1D.dsiTargetLocales() const ref {
  return targetLocales;
}

proc blockDist1D.chpl__locToLocIdx(loc: locale) {
  for locIdx in targetLocDom do
    if (targetLocales[locIdx] == loc) then
      return (true, locIdx);
  return (false, targetLocDom.first);
}

// blockDist subdomains are continuous

proc blockDistArr.dsiHasSingleLocalSubdomain() param do return !allowDuplicateTargetLocales;
proc blockDistDom.dsiHasSingleLocalSubdomain() param do return !allowDuplicateTargetLocales;

// returns the current locale's subdomain

proc blockDistArr.dsiLocalSubdomain(loc: locale) {
  if (loc == here) {
    // quick solution if we have a local array
    if const myLocArrNN = myLocArr then
      return myLocArrNN.locDom.myblockDist;
    // if not, we must not own anything
    var d: domain(rank, idxType, strides);
    return d;
  } else {
    return dom.dsiLocalSubdomain(loc);
  }
}
proc blockDistDom.dsiLocalSubdomain(loc: locale) {
  const (gotit, locid) = dist.chpl__locToLocIdx(loc);
  if (gotit) {
    if loc == here {
      // If we're doing the common case of just querying our own ownership,
      // return the local, pre-computed value
      return locDoms[locid].myblockDist;
    } else {
      // Otherwise, compute it to avoid communication...
      var inds = chpl__computeblockDist(locid, dist.targetLocDom, dist.boundingBox, dist.boundingBox.dims());
      return whole[(...inds)];
    }
  } else {
    var d: domain(rank, idxType, strides);
    return d;
  }
}

proc blockDistDom.numRemoteElems(viewDom, rlo, rid) {
  // NOTE: Not bothering to check to see if rid+1, length, or rlo-1 used
  //  below can fit into idxType
  var blo, bhi:dist.idxType;
  if rid==(dist.targetLocDom.dim(rank-1).sizeAs(int) - 1) then
    bhi=viewDom.dim(rank-1).highBound;
  else {
      bhi = dist.boundingBox.dim(rank-1).lowBound +
        intCeilXDivByY((dist.boundingBox.dim(rank-1).highBound - dist.boundingBox.dim(rank-1).lowBound +1)*(rid+1):idxType,
                       dist.targetLocDom.dim(rank-1).sizeAs(idxType)) - 1:idxType;
  }

  return (bhi - (rlo - 1):idxType);
}

private proc canDoAnyToblockDist(Dest, destDom, Src, srcDom) param : bool {
  if Src.doiCanBulkTransferRankChange() == false &&
     Dest.rank != Src.rank then return false;

  use Reflection;

  // Does 'Src' support bulk transfers *to* a DefaultRectangular?
  if !canResolveMethod(Src, "doiBulkTransferToKnown", srcDom,
                       Dest.locArr[Dest.locArr.domain.first].myElems._value, destDom) {
    return false;
  }

  return !disableblockDistDistBulkTransfer;
}

// blockDist = this
proc blockDistArr.doiBulkTransferToKnown(srcDom, destClass:blockDistArr, destDom) : bool
where this.sparseLayoutType == unmanaged DefaultDist &&
      destClass.sparseLayoutType == unmanaged DefaultDist &&
      !disableblockDistDistBulkTransfer {
  return _doSimpleblockDistTransfer(destClass, destDom, this, srcDom);
}

// this = blockDist
proc blockDistArr.doiBulkTransferFromKnown(destDom, srcClass:blockDistArr, srcDom) : bool
where this.sparseLayoutType == unmanaged DefaultDist &&
      srcClass.sparseLayoutType == unmanaged DefaultDist &&
      !disableblockDistDistBulkTransfer {
  return _doSimpleblockDistTransfer(this, destDom, srcClass, srcDom);
}

proc blockDistArr.canDoOptimizedSwap(other) {
  var domsMatch = true;

  if this.dom != other.dom { // no need to check if this is true
    if domsMatch {
      for param i in 0..this.dom.rank-1 {
        if this.dom.whole.dim(i) != other.dom.whole.dim(i) {
          domsMatch = false;
        }
      }
    }
  }

  if domsMatch {
    // distributions must be equal, too
    return this.dom.dist.dsiEqualDMaps(other.dom.dist);
  }
  return false;
}

// A helper routine that will perform a pointer swap on an array
// instead of doing a deep copy of that array. Returns true
// if used the optimized swap, false otherwise
//
// TODO: stridability causes issues with RAD swap, and somehow isn't captured by
// the formal type when we check whether this resolves.
proc blockDistArr.doiOptimizedSwap(other: this.type)
  where this.strides == other.strides {

  if(canDoOptimizedSwap(other)) {
    if debugOptimizedSwap {
      writeln("blockDistArr doing optimized swap. Domains: ",
              this.dom.whole, " ", other.dom.whole, " Bounding boxes: ",
              this.dom.dist.boundingBox, " ", other.dom.dist.boundingBox);
    }
    coforall (locarr1, locarr2) in zip(this.locArr, other.locArr) {
      on locarr1 {
        locarr1.myElems <=> locarr2.myElems;
        locarr1.locRAD <=> locarr2.locRAD;
      }
    }
    return true;
  } else {
    if debugOptimizedSwap {
      writeln("blockDistArr doing unoptimized swap. Domains: ",
              this.dom.whole, " ", other.dom.whole, " Bounding boxes: ",
              this.dom.dist.boundingBox, " ", other.dom.dist.boundingBox);
    }
    return false;
  }
}


// The purpose of this overload is to provide debugging output in the event that
// debugOptimizedSwap is on and the main routine doesn't resolve (e.g., due to a
// type, stridability, or rank mismatch in the other argument). When
// debugOptimizedSwap is off, this overload will be ignored due to its where
// clause.
pragma "last resort"
proc blockDistArr.doiOptimizedSwap(other) where debugOptimizedSwap {
  writeln("blockDistArr doing unoptimized swap. Type mismatch");
  return false;
}

private proc _doSimpleblockDistTransfer(Dest, destDom, Src, srcDom) {
  if !chpl_allStridesArePositive(Dest, destDom, Src, srcDom) then return false;

  if debugblockDistDistBulkTransfer then
    writeln("In blockDist=blockDist Bulk Transfer: Dest[", destDom, "] = Src[", srcDom, "]");

  // Cache to avoid GETs
  const DestPID = Dest.pid;
  const SrcPID = Src.pid;

  coforall i in Dest.dom.dist.activeTargetLocales(destDom) {
    on Dest.dom.dist.targetLocales[i] {
      // Relies on the fact that we privatize across all locales in the
      // program, not just the targetLocales of Dest/Src.
      const dst = if _privatization then chpl_getPrivatizedCopy(Dest.type, DestPID) else Dest;
      const src = if _privatization then chpl_getPrivatizedCopy(Src.type, SrcPID) else Src;

      // Compute the local portion of the destination domain, and find the
      // corresponding indices in the source's domain.
      const localDestblockDist = dst.dom.locDoms[i].myblockDist[destDom];
      assert(localDestblockDist.sizeAs(int) > 0);
      const corSrcblockDist    = bulkCommTranslateDomain(localDestblockDist, destDom, srcDom);
      if debugblockDistDistBulkTransfer then
        writeln("  Dest[",localDestblockDist,"] = Src[", corSrcblockDist,"]");

      // The corresponding indices in the source's domain do not necessarily
      // all live on the same locale. This loop finds the chunks that live on a
      // single locale, translates those back to the destination's domain, and
      // performs the transfer.
      for srcLoc in src.dom.dist.activeTargetLocales(corSrcblockDist) {
        const localSrcChunk  = corSrcblockDist[src.dom.locDoms[srcLoc].myblockDist];
        const localDestChunk = bulkCommTranslateDomain(localSrcChunk, corSrcblockDist, localDestblockDist);
        chpl__bulkTransferArray(dst.locArr[i].myElems._value, localDestChunk,
                                src.locArr[srcLoc].myElems._value, localSrcChunk);
      }
    }
  }

  return true;
}

// Overload for any transfer *to* blockDist, if the RHS supports transfers to a
// DefaultRectangular
proc blockDistArr.doiBulkTransferFromAny(destDom, Src, srcDom) : bool
where canDoAnyToblockDist(this, destDom, Src, srcDom) {
  if !chpl_allStridesArePositive(this, destDom, Src, srcDom) then return false;

  if debugblockDistDistBulkTransfer then
    writeln("In blockDistDist.doiBulkTransferFromAny");

  coforall j in dom.dist.activeTargetLocales(destDom) {
    on dom.dist.targetLocales(j) {
      const Dest = if _privatization then chpl_getPrivatizedCopy(this.type, pid) else this;

      const inters   = Dest.dom.locDoms(j).myblockDist[destDom];
      const srcChunk = bulkCommTranslateDomain(inters, destDom, srcDom);

      if debugblockDistDistBulkTransfer then
        writeln("Dest.locArr[", j, "][", inters, "] = Src[", srcDom, "]");

      chpl__bulkTransferArray(Dest.locArr[j].myElems._value, inters, Src, srcChunk);
    }
  }

  return true;
}

// For assignments of the form: DefaultRectangular = blockDist
proc blockDistArr.doiBulkTransferToKnown(srcDom, Dest:DefaultRectangularArr, destDom) : bool
where !disableblockDistDistBulkTransfer {
  if !chpl_allStridesArePositive(this, srcDom, Dest, destDom) then return false;

  if debugblockDistDistBulkTransfer then
    writeln("In blockDistDist.doiBulkTransferToKnown(DefaultRectangular)");

  coforall j in dom.dist.activeTargetLocales(srcDom) {
    on dom.dist.targetLocales(j) {
      const Src = if _privatization then chpl_getPrivatizedCopy(this.type, pid) else this;
      const inters = Src.dom.locDoms(j).myblockDist[srcDom];

      const destChunk = bulkCommTranslateDomain(inters, srcDom, destDom);

      if debugblockDistDistBulkTransfer then
        writeln("  A[",destChunk,"] = B[",inters,"]");

      const elemActual = Src.locArr[j].myElems._value;
      chpl__bulkTransferArray(Dest, destChunk, elemActual, inters);
    }
  }

  return true;
}

// For assignments of the form: blockDist = DefaultRectangular
proc blockDistArr.doiBulkTransferFromKnown(destDom, Src:DefaultRectangularArr, srcDom) : bool
where !disableblockDistDistBulkTransfer {
  if !chpl_allStridesArePositive(this, destDom, Src, srcDom) then return false;

  if debugblockDistDistBulkTransfer then
    writeln("In blockDistArr.doiBulkTransferFromKnown(DefaultRectangular)");

  coforall j in dom.dist.activeTargetLocales(destDom) {
    on dom.dist.targetLocales(j) {
      // Grab privatized copy of 'this' to avoid extra GETs
      const Dest = if _privatization then chpl_getPrivatizedCopy(this.type, pid) else this;
      const inters = Dest.dom.locDoms(j).myblockDist[destDom];
      assert(inters.sizeAs(int) > 0);

      const srcChunk = bulkCommTranslateDomain(inters, destDom, srcDom);

      if debugblockDistDistBulkTransfer then
        writeln("  A[",inters,"] = B[",srcChunk,"]");

      const elemActual = Dest.locArr[j].myElems._value;
      chpl__bulkTransferArray(elemActual, inters, Src, srcChunk);
    }
  }

  return true;
}

override proc blockDistArr.doiCanBulkTransferRankChange() param do return true;

config param debugblockDistScan = false;

/* Boxed sync with readFE/writeEF only, but works for arbitrary types. Not
 * suitable for general use since there are races with when the value gets
 * written to, but safe for single writer, single reader case here.
 */
class BoxedSync {
  type T;
  var s: sync int; // int over bool to enable native qthread sync
  var res: T;

  proc readFE(): T {
    s.readFE();
    return res;
  }

  proc writeEF(val: T): void {
    res = val;
    s.writeEF(1);
  }

  // guard against dynamic dispatch trying to resolve write()ing the sync
  override proc writeThis(f) throws { }
}

proc blockDistArr.doiScan(op, dom) where (rank == 1) &&
                                     chpl__scanStateResTypesMatch(op) {

  // The result of this scan, which will be blockDist-distributed as well
  type resType = op.generate().type;
  var res = dom.buildArray(resType, initElts=!isPOD(resType));

  // Store one element per locale in order to track our local total
  // for a cross-locale scan as well as flags to negotiate reading and
  // writing it.
  const ref targetLocs = this.dsiTargetLocales();
  const firstLoc = targetLocs.domain.lowBound;
  var elemPerLoc: [targetLocs.domain] resType;
  var inputReady: [targetLocs.domain] sync int;
  var outputReady: [targetLocs.domain] unmanaged BoxedSync(resType)?;

  const cachedPID = pid;

  param sameStaticDist = dsiStaticFastFollowCheck(res._value.type);
  const sameDynamicDist = sameStaticDist && dsiDynamicFastFollowCheck(res);

  // Fire up tasks per participating locale
  coforall locid in targetLocs.domain {
    on targetLocs[locid] {
      const myop = op.clone(); // this will be deleted by doiScan()

      // set up some references to our LocblockDistArr descriptor, our
      // local array, local domain, and local result elements
      const thisArr = if _privatization then chpl_getPrivatizedCopy(this.type, cachedPID) else this;
      ref myLocArr = if allowDuplicateTargetLocales then thisArr.locArr[locid].myElems
                                                        else thisArr.myLocArr!.myElems;
      const ref myLocDom = myLocArr.domain;

      // Compute the local pre-scan on our local array
      const (numTasks, rngs, state, tot) =
        if !allowDuplicateTargetLocales && sameStaticDist && sameDynamicDist then
          myLocArr._value.chpl__preScan(myop, res._value.myLocArr!.myElems, myLocDom[dom])
        else
          myLocArr._value.chpl__preScan(myop, res, myLocDom[dom]);

      if debugblockDistScan then
        writeln(locid, ": ", (numTasks, rngs, state, tot));

      // Create a local ready var and store it back on the initiating
      // locale so it can notify us
      const myOutputReady = new unmanaged BoxedSync(resType);

      on elemPerLoc {
        // save our local scan total away and signal that it's ready
        outputReady[locid] = myOutputReady;
        elemPerLoc[locid] = tot;
        inputReady[locid].writeEF(1);
      }

      // the "first" locale scans the per-locale contributions as they
      // become ready
      if (locid == firstLoc) then on elemPerLoc {
        const metaop = op.clone();

        var next: resType = metaop.identity;
        for locid in targetLocs.domain {
          inputReady[locid].readFE();

          // store the scan value
          ref locVal = elemPerLoc[locid];
          locVal <=> next;

          // accumulate to prep for the next iteration
          metaop.accumulateOntoState(next, locVal);
        }

        // Iterator that yields values instead of references (to enable RVF)
        iter valIter(iterable) where isPOD(iterable.eltType) {
          for elem in iterable do yield elem;
        }

        // BigInteger values are yielded by ref to disable RVF
        iter valIter(iterable) ref where !isPOD(iterable.eltType) {
          for elem in iterable do yield elem;
        }

        // Mark that scan values are ready
        coforall (_, ready, elem) in zip(targetLocs.domain, valIter(outputReady), valIter(elemPerLoc)) do on ready {
          ready!.writeEF(elem);
        }

        delete metaop;
      }

      // block until our local value has been updated
      const myadjust = myOutputReady.readFE();
      if debugblockDistScan then
        writeln(locid, ": myadjust = ", myadjust);

      // update our state vector with our locale's adjustment value
      for s in state do
        myop.accumulateOntoState(s, myadjust);
      if debugblockDistScan then
        writeln(locid, ": state = ", state);

      // have our local array compute its post scan with the globally
      // accurate state vector
      if !allowDuplicateTargetLocales && sameStaticDist && sameDynamicDist then
        myLocArr._value.chpl__postScan(op, res._value.myLocArr!.myElems, numTasks, rngs, state);
      else
        myLocArr._value.chpl__postScan(op, res, numTasks, rngs, state);

      if debugblockDistScan then
        writeln(locid, ": ", myLocArr);

      delete myop;
      delete myOutputReady;
    }
  }
  if isPOD(resType) then res.dsiElementInitializationComplete();

  delete op;
  return res;
}

////// DSI Utils overrides ////////////////////////////////////////////////////

proc setupTargetLocRanges(param rank, specifiedLocArr) {
  var ranges: rank*range;

  if specifiedLocArr.rank == 1 {
    ranges(0) = 0..#specifiedLocArr.size;
    for param i in 1..rank-1 do
      ranges(i) = 0..#1;

  } else {
    if specifiedLocArr.rank != 1 then
      compilerError("Specified locale array has rank different then 1");
  }

  return ranges;
}
