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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "edgeExtractor.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "triSurfaceMetaData.H"
#include "polyMeshGenGeometryModification.H"
#include "surfaceMeshGeometryModification.H"

#include <cstdlib>
#include <cstdint>
#include <cstring>

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

inline bool surfaceProjectionDebugEnabled()
{
    return std::getenv("CFMESH_SURFACE_PROJECTION_DEBUG_DIGEST") != nullptr;
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

void reportSurfaceProjectionDigest(const char* name, const uint64_t digest)
{
    if( surfaceProjectionDebugEnabled() )
    {
        Info<< "surfaceProjectionDigest "
            << name << " 0x" << hex << digest << dec << endl;
    }
}

void reportMeshDigest(polyMeshGen& mesh, const char* name)
{
    if( !surfaceProjectionDebugEnabled() )
        return;

    const pointFieldPMG& points = mesh.points();
    uint64_t meshPointDigest(1469598103934665603ULL);
    forAll(points, pointI)
    {
        hashCombineU64(meshPointDigest, static_cast<uint64_t>(pointI));
        hashCombineU64(meshPointDigest, pointDigest(points[pointI]));
    }

    reportSurfaceProjectionDigest(name, meshPointDigest);
}

void stopAfterMeshOptimisationSubstep
(
    polyMeshGen& mesh,
    const char* substepName
)
{
    const char* requested = std::getenv("CFMESH_STOP_AFTER_OPT_SUBSTEP");
    if( !requested || (std::strcmp(requested, substepName) != 0) )
        return;

    Info << "Saving mesh after meshOptimisation substep "
         << substepName << endl;

    mesh.write();

    std::string message("Stopping after meshOptimisation substep ");
    message += substepName;
    throw message;
}

void stopAfterGeneratorStep
(
    polyMeshGen& mesh,
    const char* stepName
)
{
    const char* requested = std::getenv("CFMESH_STOP_AFTER_GENERATOR_STEP");
    if( !requested || (std::strcmp(requested, stepName) != 0) )
        return;

    Info << "Saving mesh after generator step "
         << stepName << endl;

    mesh.write();

    std::string message("Stopping after generator step ");
    message += stepName;
    throw message;
}

void stopAfterSurfaceProjectionSubstep
(
    polyMeshGen& mesh,
    const char* substepName
)
{
    const char* requested =
        std::getenv("CFMESH_STOP_AFTER_SURFACE_PROJECTION_SUBSTEP");
    if( !requested || (std::strcmp(requested, substepName) != 0) )
        return;

    Info << "Saving mesh after surfaceProjection substep "
         << substepName << endl;

    mesh.write();

    std::string message("Stopping after surfaceProjection substep ");
    message += substepName;
    throw message;
}

}

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void cartesianMeshGenerator::createCartesianMesh()
{
    //- create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);
    reportMeshDigest(mesh_, "beforeCreateCartesianMesh");

    if( meshDict_.found("decomposePolyhedraIntoTetsAndPyrs") )
    {
        if( readBool(meshDict_.lookup("decomposePolyhedraIntoTetsAndPyrs")) )
            cme.decomposeSplitHexes();
    }

    cme.createMesh();
    reportMeshDigest(mesh_, "afterCreateCartesianMesh");
}

void cartesianMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    bool changed;

    reportMeshDigest(mesh_, "surfacePreparationStart");

    do
    {
        changed = false;

        checkIrregularSurfaceConnections checkConnections(mesh_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;
        reportMeshDigest(mesh_, "afterCheckIrregularSurfaceConnections");

        if( checkNonMappableCellConnections(mesh_).removeCells() )
            changed = true;
        reportMeshDigest(mesh_, "afterCheckNonMappableCellConnections");

        if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
            changed = true;
        reportMeshDigest(mesh_, "afterCheckCellConnectionsOverFaces");
    } while( changed );

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();
    reportMeshDigest(mesh_, "afterCheckBoundaryFacesSharingTwoEdges");
}

void cartesianMeshGenerator::mapMeshToSurface()
{
    if( surfaceProjectionDebugEnabled() )
    {
        const pointFieldPMG& points = mesh_.points();
        uint64_t meshPointDigest(1469598103934665603ULL);
        forAll(points, pointI)
        {
            hashCombineU64(meshPointDigest, static_cast<uint64_t>(pointI));
            hashCombineU64(meshPointDigest, pointDigest(points[pointI]));
        }
        reportSurfaceProjectionDigest("meshPointsBeforeMse", meshPointDigest);
    }

    //- calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    if( surfaceProjectionDebugEnabled() )
    {
        const pointFieldPMG& points = mesh_.points();
        uint64_t meshPointDigest(1469598103934665603ULL);
        forAll(points, pointI)
        {
            hashCombineU64(meshPointDigest, static_cast<uint64_t>(pointI));
            hashCombineU64(meshPointDigest, pointDigest(points[pointI]));
        }
        reportSurfaceProjectionDigest("meshPointsAfterMse", meshPointDigest);

        const labelList& boundaryPoints = mse.boundaryPoints();
        uint64_t boundaryPointDigest(1469598103934665603ULL);
        forAll(boundaryPoints, bpI)
        {
            const label pointI = boundaryPoints[bpI];
            hashCombineU64(boundaryPointDigest, static_cast<uint64_t>(bpI));
            hashCombineU64(boundaryPointDigest, static_cast<uint64_t>(pointI));
            hashCombineU64(boundaryPointDigest, pointDigest(points[pointI]));
        }
        reportSurfaceProjectionDigest("boundaryPoints", boundaryPointDigest);

        const vectorField& faceCentres = mse.faceCentres();
        uint64_t faceCentresDigest(1469598103934665603ULL);
        forAll(faceCentres, faceI)
        {
            hashCombineU64(faceCentresDigest, static_cast<uint64_t>(faceI));
            hashCombineU64(faceCentresDigest, pointDigest(faceCentres[faceI]));
        }
        reportSurfaceProjectionDigest("faceCentres", faceCentresDigest);
    }

    //- pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);
    mapper.preMapVertices();
    stopAfterSurfaceProjectionSubstep(mesh_, "preMapVertices");

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();
    stopAfterSurfaceProjectionSubstep(mesh_, "mapVerticesOntoSurface");

    //- untangle surface faces
    meshSurfaceOptimizer(mse, *octreePtr_).untangleSurface();
    stopAfterSurfaceProjectionSubstep(mesh_, "untangleSurface");
}

void cartesianMeshGenerator::extractPatches()
{
    edgeExtractor extractor(mesh_, *octreePtr_);

    Info << "Extracting edges" << endl;
    extractor.extractEdges();

    extractor.updateMeshPatches();
}

void cartesianMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);
}

void cartesianMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();
}

void cartesianMeshGenerator::generateBoundaryLayers()
{
    //- add boundary layers
    boundaryLayers bl(mesh_);
    bl.addLayerForAllPatches();
}

void cartesianMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.refineLayers();

        labelLongList pointsInLayer;
        refLayers.pointsInBndLayer(pointsInLayer);

        meshOptimizer mOpt(mesh_);
        mOpt.lockPoints(pointsInLayer);
        mOpt.untangleBoundaryLayer();
    }
}

void cartesianMeshGenerator::optimiseFinalMesh()
{
    //- untangle the surface if needed
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

    if( true )
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);

        if( enforceConstraints )
            surfOpt.enforceConstraints();

        surfOpt.optimizeSurface();
        stopAfterMeshOptimisationSubstep(mesh_, "surfaceOptimize");
    }

    deleteDemandDrivenData(octreePtr_);

    //- final optimisation
    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();

    optimizer.optimizeMeshFV();
    stopAfterMeshOptimisationSubstep(mesh_, "optimizeMeshFV");
    optimizer.optimizeLowQualityFaces();
    stopAfterMeshOptimisationSubstep(mesh_, "optimizeLowQualityFaces");
    optimizer.optimizeBoundaryLayer(modSurfacePtr_==NULL);
    stopAfterMeshOptimisationSubstep(mesh_, "optimizeBoundaryLayer");
    optimizer.untangleMeshFV();
    stopAfterMeshOptimisationSubstep(mesh_, "untangleMeshFV");

    mesh_.clearAddressingData();

    if( modSurfacePtr_ )
    {
        polyMeshGenGeometryModification meshMod(mesh_, meshDict_);

        //- revert the mesh into the original space
        meshMod.revertGeometryModification();

        //- delete modified surface mesh
        deleteDemandDrivenData(modSurfacePtr_);
    }
}

void cartesianMeshGenerator::projectSurfaceAfterBackScaling()
{
    if( !meshDict_.found("anisotropicSources") )
        return;

    deleteDemandDrivenData(octreePtr_);
    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator
    (
        *octreePtr_,
        meshDict_
    ).createOctreeWithRefinedBoundary(20, 30);

    //- calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    //- pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    optimiseFinalMesh();
}

void cartesianMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);
}

void cartesianMeshGenerator::fillInCellAndPointLevels()
{
    mesh_.fillInCellAndPointLevels();
}

void cartesianMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}

void cartesianMeshGenerator::generateMesh()
{
    if( controller_.runCurrentStep("templateGeneration") )
    {
        createCartesianMesh();
        stopAfterGeneratorStep(mesh_, "templateGeneration");
    }

    if( controller_.runCurrentStep("surfaceTopology") )
    {
        surfacePreparation();
        stopAfterGeneratorStep(mesh_, "surfaceTopology");
    }

    if( controller_.runCurrentStep("surfaceProjection") )
    {
        mapMeshToSurface();
        stopAfterGeneratorStep(mesh_, "surfaceProjection");
    }

    if( controller_.runCurrentStep("patchAssignment") )
    {
        extractPatches();
        stopAfterGeneratorStep(mesh_, "patchAssignment");
    }

    if( controller_.runCurrentStep("edgeExtraction") )
    {
        mapEdgesAndCorners();

        optimiseMeshSurface();
        stopAfterGeneratorStep(mesh_, "edgeExtraction");
    }

    if( controller_.runCurrentStep("boundaryLayerGeneration") )
    {
        generateBoundaryLayers();
        stopAfterGeneratorStep(mesh_, "boundaryLayerGeneration");
    }

    if( controller_.runCurrentStep("meshOptimisation") )
    {
        optimiseFinalMesh();

        projectSurfaceAfterBackScaling();
        stopAfterGeneratorStep(mesh_, "meshOptimisation");
    }

    if( controller_.runCurrentStep("boundaryLayerRefinement") )
    {
        refBoundaryLayers();
        stopAfterGeneratorStep(mesh_, "boundaryLayerRefinement");
    }

    fillInCellAndPointLevels();

    renumberMesh();

    replaceBoundaries();

    controller_.workflowCompleted();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cartesianMeshGenerator::cartesianMeshGenerator(const Time& time)
:
    db_(time),
    surfacePtr_(NULL),
    modSurfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.system(),
            db_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(time),
    controller_(mesh_)
{
    try
    {
        if( true )
        {
            checkMeshDict cmd(meshDict_);
        }

        fileName surfaceFile(meshDict_.lookup("surfaceFile"));
        if( Pstream::parRun() )
            surfaceFile = ".."/surfaceFile;

        surfacePtr_ = new triSurf(db_.path()/surfaceFile);

        if( true )
        {
            //- save meta data with the mesh (surface mesh + its topology info)
            triSurfaceMetaData sMetaData(*surfacePtr_);
            const dictionary& surfMetaDict = sMetaData.metaData();

            mesh_.metaData().add("surfaceFile", surfaceFile, true);
            mesh_.metaData().add("surfaceMeta", surfMetaDict, true);
        }

        if( surfacePtr_->featureEdges().size() != 0 )
        {
            //- create surface patches based on the feature edges
            //- and update the meshDict based on the given data
            triSurfacePatchManipulator manipulator(*surfacePtr_);

            const triSurf* surfaceWithPatches =
                manipulator.surfaceWithPatches(&meshDict_);

            //- delete the old surface and assign the new one
            deleteDemandDrivenData(surfacePtr_);
            surfacePtr_ = surfaceWithPatches;
        }

        if( meshDict_.found("anisotropicSources") )
        {
            surfaceMeshGeometryModification surfMod(*surfacePtr_, meshDict_);

            modSurfacePtr_ = surfMod.modifyGeometry();

            octreePtr_ = new meshOctree(*modSurfacePtr_);
        }
        else
        {
            octreePtr_ = new meshOctree(*surfacePtr_);
        }

        meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

        generateMesh();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        WarningIn
        (
            "cartesianMeshGenerator::cartesianMeshGenerator(const Time&)"
        ) << "Meshing process terminated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshGenerator::writeMesh() const
{
    mesh_.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
