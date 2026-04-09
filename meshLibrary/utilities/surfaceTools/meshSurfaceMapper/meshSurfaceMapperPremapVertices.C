/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "helperFunctions.H"
#include "meshSurfaceOptimizer.H"

#include <cstdlib>
#include <cstdint>
#include <cstring>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

inline bool preMapDebugEnabled()
{
    return std::getenv("CFMESH_PREMAP_DEBUG_DIGEST") != nullptr;
}

inline void hashCombineU64(uint64_t& hash, const uint64_t value)
{
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
}

inline uint64_t scalarBits(const scalar value)
{
    uint64_t bits(0);
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline uint64_t pointDigest(const point& p)
{
    uint64_t hash(1469598103934665603ULL);
    hashCombineU64(hash, scalarBits(p.x()));
    hashCombineU64(hash, scalarBits(p.y()));
    hashCombineU64(hash, scalarBits(p.z()));
    return hash;
}

void reportPreMapDigest(const char* name, const uint64_t digest)
{
    if( preMapDebugEnabled() )
        Info<< "preMapDigest " << name << " 0x" << hex << digest << dec << endl;
}

void stopAfterPreMapSubstep
(
    polyMeshGen& mesh,
    const char* substepName
)
{
    const char* requested = std::getenv("CFMESH_STOP_AFTER_PREMAP_SUBSTEP");
    if( !requested || (std::strcmp(requested, substepName) != 0) )
        return;

    Info << "Saving mesh after preMapVertices substep "
         << substepName << endl;

    mesh.write();

    std::string message("Stopping after preMapVertices substep ");
    message += substepName;
    throw message;
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::preMapVertices(const label nIterations)
{
    Info << "Smoothing mesh surface before mapping." << endl;

    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const VRWGraph& pointInFace = surfaceEngine_.pointInFaces();
    const labelList& bp = surfaceEngine_.bp();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    const triSurf& surf = meshOctree_.surface();

    List<labelledPointScalar> preMapPositions(boundaryPoints.size());
    List<DynList<scalar, 6> > faceCentreDistances(bFaces.size());

    if( preMapDebugEnabled() )
    {
        uint64_t inputPointDigest(1469598103934665603ULL);
        forAll(boundaryPoints, bpI)
        {
            hashCombineU64
            (
                inputPointDigest,
                static_cast<uint64_t>(boundaryPoints[bpI])
            );
            hashCombineU64
            (
                inputPointDigest,
                pointDigest(points[boundaryPoints[bpI]])
            );
        }
        reportPreMapDigest("inputPoints", inputPointDigest);

        uint64_t pointFacesDigest(1469598103934665603ULL);
        forAll(boundaryPoints, bpI)
        {
            hashCombineU64
            (
                pointFacesDigest,
                static_cast<uint64_t>(boundaryPoints[bpI])
            );
            forAllRow(pointFaces, bpI, pfI)
            {
                hashCombineU64
                (
                    pointFacesDigest,
                    static_cast<uint64_t>(pointFaces(bpI, pfI))
                );
                hashCombineU64
                (
                    pointFacesDigest,
                    static_cast<uint64_t>(pointInFace(bpI, pfI))
                );
            }
        }
        reportPreMapDigest("pointFaces", pointFacesDigest);

        uint64_t faceCentresDigest(1469598103934665603ULL);
        forAll(faceCentres, faceI)
        {
            hashCombineU64(faceCentresDigest, static_cast<uint64_t>(faceI));
            hashCombineU64(faceCentresDigest, pointDigest(faceCentres[faceI]));
        }
        reportPreMapDigest("faceCentres", faceCentresDigest);
    }

    # ifdef USE_OMP
    # pragma omp parallel for if(false) schedule(dynamic, 20)
    # endif
    forAll(bFaces, bfI)
    {
        const point& c = faceCentres[bfI];
        const face& bf = bFaces[bfI];

        faceCentreDistances[bfI].setSize(bf.size());

        forAll(bf, pI)
        {
            faceCentreDistances[bfI][pI] = magSqr(points[bf[pI]] - c);
        }
    }

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        //- find patches in the vicinity of a boundary face
        List<DynList<label> > boundaryPointPatches(boundaryPoints.size());
        LongList<labelPair> boundaryPointPatchPairs;

        # ifdef USE_OMP
        # pragma omp parallel if(false)
        # endif
        {
            LongList<labelPair> localPairs;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 50)
            # endif
            forAll(bFaces, bfI)
            {
                const face& bf = bFaces[bfI];
                scalar boxSize(0.0);
                forAll(bf, pI)
                {
                    boxSize =
                        Foam::max
                        (
                            boxSize,
                            mag(faceCentres[bfI] - points[bf[pI]])
                        );
                }

                const boundBox bb
                (
                    faceCentres[bfI] - vector(boxSize, boxSize, boxSize),
                    faceCentres[bfI] + vector(boxSize, boxSize, boxSize)
                );

                DynList<label> containedLeaves;
                meshOctree_.findLeavesContainedInBox(bb, containedLeaves);

                DynList<label> patches;
                forAll(containedLeaves, clI)
                {
                    DynList<label> ct;
                    meshOctree_.containedTriangles(containedLeaves[clI], ct);

                    forAll(ct, i)
                        patches.appendIfNotIn(surf[ct[i]].region());
                }

                labelList sortedPatches(patches.size());
                forAll(sortedPatches, i)
                    sortedPatches[i] = patches[i];
                sort(sortedPatches);

                scalar metric(VGREAT);
                label bestPatch(-1);
                forAll(sortedPatches, ptchI)
                {
                    const scalar m =
                        faceMetricInPatch(bfI, sortedPatches[ptchI]);

                    if( m < metric )
                    {
                        metric = m;
                        bestPatch = sortedPatches[ptchI];
                    }
                }

                forAll(bf, pI)
                    localPairs.append(labelPair(bp[bf[pI]], bestPatch));
            }

            # ifdef USE_OMP
            # pragma omp critical
            # endif
            forAll(localPairs, i)
                boundaryPointPatchPairs.append(localPairs[i]);
        }

        List<labelPair> orderedBoundaryPointPatchPairs
        (
            boundaryPointPatchPairs.size()
        );
        forAll(orderedBoundaryPointPatchPairs, i)
            orderedBoundaryPointPatchPairs[i] = boundaryPointPatchPairs[i];

        sort
        (
            orderedBoundaryPointPatchPairs,
            [](const labelPair& a, const labelPair& b)
            {
                if( a.first() < b.first() )
                    return true;
                if( a.first() > b.first() )
                    return false;

                return a.second() < b.second();
            }
        );

        forAll(orderedBoundaryPointPatchPairs, i)
        {
            const labelPair& pp = orderedBoundaryPointPatchPairs[i];
            boundaryPointPatches[pp.first()].appendIfNotIn(pp.second());
        }

        if( preMapDebugEnabled() )
        {
            uint64_t patchDigest(1469598103934665603ULL);
            forAll(boundaryPointPatches, bpI)
            {
                hashCombineU64(patchDigest, static_cast<uint64_t>(bpI));
                forAll(boundaryPointPatches[bpI], i)
                    hashCombineU64
                    (
                        patchDigest,
                        static_cast<uint64_t>(boundaryPointPatches[bpI][i] + 1)
                    );
            }

            reportPreMapDigest("patches", patchDigest);
        }

        //- use the shrinking laplace first
        # ifdef USE_OMP
        # pragma omp parallel for if(false) schedule(dynamic, 40)
        # endif
        forAll(pointFaces, bpI)
        {
            labelledPointScalar lp(bpI, vector::zero, 0.0);

            const point& p = points[boundaryPoints[bpI]];
            LongList<labelPair> faceOrder;

            forAllRow(pointFaces, bpI, pfI)
                faceOrder.append(labelPair(pointFaces(bpI, pfI), pointInFace(bpI, pfI)));

            List<labelPair> sortedFaceOrder(faceOrder.size());
            forAll(sortedFaceOrder, i)
                sortedFaceOrder[i] = faceOrder[i];

            sort
            (
                sortedFaceOrder,
                [](const labelPair& a, const labelPair& b)
                {
                    if( a.first() < b.first() )
                        return true;
                    if( a.first() > b.first() )
                        return false;

                    return a.second() < b.second();
                }
            );

            forAll(sortedFaceOrder, i)
            {
                const label bfI = sortedFaceOrder[i].first();
                const point& fc = faceCentres[bfI];
                const label pos = sortedFaceOrder[i].second();
                const scalar w
                (
                    max(magSqr(p - fc) / faceCentreDistances[bfI][pos], SMALL)
                );
                lp.coordinates() += w * fc;
                lp.scalarValue() += w;
            }

            preMapPositions[bpI] = lp;
        }

        if( preMapDebugEnabled() )
        {
            uint64_t preMapDigest(1469598103934665603ULL);
            forAll(preMapPositions, bpI)
            {
                hashCombineU64(preMapDigest, static_cast<uint64_t>(bpI));
                hashCombineU64
                (
                    preMapDigest,
                    pointDigest(preMapPositions[bpI].coordinates())
                );
                hashCombineU64
                (
                    preMapDigest,
                    scalarBits(preMapPositions[bpI].scalarValue())
                );
            }

            reportPreMapDigest("laplace", preMapDigest);
        }

        //- pointer needed in case of parallel calculation
        const VRWGraph* bpAtProcsPtr(NULL);

        if( Pstream::parRun() )
        {
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            bpAtProcsPtr = &bpAtProcs;
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();

            //- collect data to be sent to other processors
            std::map<label, LongList<labelledPointScalar> > exchangeData;
            forAll(surfaceEngine_.bpNeiProcs(), i)
                exchangeData.insert
                (
                    std::make_pair
                    (
                        surfaceEngine_.bpNeiProcs()[i],
                        LongList<labelledPointScalar>()
                    )
                );

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append
                    (
                        labelledPointScalar
                        (
                            globalPointLabel[bpI],
                            preMapPositions[bpI].coordinates(),
                            preMapPositions[bpI].scalarValue()
                        )
                    );
                }
            }

            //- exchange data with other processors
            LongList<labelledPointScalar> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- combine collected data with the available data
            forAll(receivedData, i)
            {
                const labelledPointScalar& lps = receivedData[i];

                const label bpI = globalToLocal[lps.pointLabel()];

                labelledPointScalar& lp = preMapPositions[bpI];
                lp.coordinates() += lps.coordinates();
                lp.scalarValue() += lps.scalarValue();
            }
        }

        //- create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
        LongList<parMapperHelper> parallelBndNodes;

        # ifdef USE_OMP
        # pragma omp parallel for if(false) schedule(dynamic, 50)
        # endif
        forAll(boundaryPoints, bpI)
        {
            labelledPointScalar& lps = preMapPositions[bpI];

            lps.coordinates() /= lps.scalarValue();

            const point& p = points[boundaryPoints[bpI]];

            //label patch, nearestTri;
            point pMap = p;
            scalar dSq;

            if( boundaryPointPatches[bpI].size() == 1 )
            {
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    pMap,
                    dSq,
                    nt,
                    boundaryPointPatches[bpI][0],
                    lps.coordinates()
                );
            }
            else
            {
                meshOctree_.findNearestPointToPatches
                (
                    pMap,
                    dSq,
                    lps.coordinates(),
                    boundaryPointPatches[bpI]
                );
            }

            const point newP = p + 0.5 * (pMap - p);

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);

            if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                parallelBndNodes.append
                (
                    parMapperHelper
                    (
                        newP,
                        dSq,
                        bpI,
                        -1
                    )
                );
            }
        }

        if( preMapDebugEnabled() )
        {
            uint64_t movedDigest(1469598103934665603ULL);
            const pointFieldPMG& movedPoints = surfaceEngine_.points();
            forAll(boundaryPoints, bpI)
            {
                hashCombineU64
                (
                    movedDigest,
                    static_cast<uint64_t>(boundaryPoints[bpI])
                );
                hashCombineU64
                (
                    movedDigest,
                    pointDigest(movedPoints[boundaryPoints[bpI]])
                );
            }

            reportPreMapDigest("moved", movedDigest);
        }

        //- make sure that the vertices at inter-processor boundaries
        //- are mapped onto the same location
        mapToSmallestDistance(parallelBndNodes);
        stopAfterPreMapSubstep
        (
            const_cast<polyMeshGen&>(surfaceEngine_.mesh()),
            "afterMoveBeforeUpdate"
        );

        //- update the surface geometry of the
        surfaceModifier.updateGeometry();
        stopAfterPreMapSubstep
        (
            const_cast<polyMeshGen&>(surfaceEngine_.mesh()),
            "beforeUntangle"
        );

        meshSurfaceOptimizer(surfaceEngine_, meshOctree_).untangleSurface();
        stopAfterPreMapSubstep
        (
            const_cast<polyMeshGen&>(surfaceEngine_.mesh()),
            "afterUntangle"
        );

        surfaceModifier.updateGeometry();
    }

    Info << "Finished smoothing mesh surface before mapping." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
