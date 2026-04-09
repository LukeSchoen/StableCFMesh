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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"
#include "labelList.H"
#include <cstdlib>
#include <cstdint>
#include <cstring>

namespace Foam
{

namespace
{

inline bool templateDebugEnabled()
{
    return std::getenv("CFMESH_TEMPLATE_DEBUG_DIGEST") != nullptr;
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

void reportTemplateMeshDigest(polyMeshGen& mesh, const char* stageName)
{
    if( !templateDebugEnabled() )
        return;

    const pointFieldPMG& points = mesh.points();
    uint64_t digest(1469598103934665603ULL);
    forAll(points, pointI)
    {
        hashCombineU64(digest, static_cast<uint64_t>(pointI));
        hashCombineU64(digest, pointDigest(points[pointI]));
    }

    Info<< "templateDigest " << stageName
        << " 0x" << hex << digest << dec << endl;
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::removeUnusedVertices()
{
    reportTemplateMeshDigest(mesh_, "beforeRemoveUnusedVertices");

    faceListPMG& faces = mesh_.faces_;
    pointFieldPMG& points = mesh_.points_;
    labelIOList& pointLevel = mesh_.pointLevel_;

    boolList usePoint(points.size(), false);
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        forAll(f, pI)
            usePoint[f[pI]] = true;
    }

    labelLongList newLabel(points.size(), -1);
    label nPoints(0);
    forAll(points, pI)
        if( usePoint[pI] )
            newLabel[pI] = nPoints++;

    //- remove unused points from the list
    forAll(newLabel, pI)
        if( (newLabel[pI] != -1) && (newLabel[pI] < pI) )
        {
            points[newLabel[pI]] = points[pI];
            if (pointLevel.size())
                pointLevel[newLabel[pI]] = pointLevel[pI];
        }
    
    points.setSize(nPoints);
    if(pointLevel.size())
        pointLevel.setSize(nPoints, -1);

    forAll(faces, faceI)
    {
        face& f = faces[faceI];

        forAll(f, pI)
            f[pI] = newLabel[f[pI]];
    }

    mesh_.updatePointSubsets(newLabel);

    mesh_.clearOut();
    this->clearOut();

    reportTemplateMeshDigest(mesh_, "afterRemoveUnusedVerticesModifier");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
