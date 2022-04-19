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
#include "VRWGraphList.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::addCells(const LongList<faceList>& cellFaces)
{
    Info << "Adding cells to the mesh" << endl;

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;
    labelIOList& cellLevel = mesh_.cellLevel_;

    VRWGraph& pointFaces = this->pointFaces();

    //- start adding cells into the mesh
    label nFaces = faces.size();
    label nCells = cells.size();

    forAll(cellFaces, cI)
    {
        const faceList& facesInCell = cellFaces[cI];

        label fI(0);
        cell c(cellFaces[cI].size());

        forAll(facesInCell, faceI)
        {
            const face& f = facesInCell[faceI];

            const label pointI = f[0];

            label fLabel(-1);
            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                if( faces[faceI] == f )
                {
                    fLabel = faceI;
                    break;
                }
            }

            if( fLabel == -1 )
            {
                faces.append(f);
                c[fI++] = nFaces;

                forAll(f, pI)
                    pointFaces.append(f[pI], nFaces);

                ++nFaces;
            }
            else
            {
                c[fI++] = fLabel;
            }
        }

        cells.append(c);
        cellLevel.append(-1); //TODO: Get cell level from neighbour?
        ++nCells;
    }

    this->clearOut();
    mesh_.clearOut();

    Info << "Finished adding cells to the mesh" << endl;
}

void polyMeshGenModifier::addCells(const VRWGraphList& cellFaces)
{
    Info << "Adding " << cellFaces.size() << " cells to the mesh" << endl;

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;
    labelIOList& cellLevel = mesh_.cellLevel_;

    VRWGraph& pointFaces = this->pointFaces();

    //- start adding cells into the mesh
    label nFaces = faces.size();
    label nCells = cells.size();

/*
    # ifdef USE_OMP
    //# pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(cellFaces, cI)
    {
        faceList facesInCell(cellFaces.sizeOfGraph(cI));
        forAll(facesInCell, fI)
        {
            facesInCell[fI].setSize(cellFaces.sizeOfRow(cI, fI));

            forAll(facesInCell[fI], pI)
                facesInCell[fI][pI] = cellFaces(cI, fI, pI);
        }

        label fI(0);
        cell c(facesInCell.size());

        forAll(facesInCell, faceI)
        {
            const face& f = facesInCell[faceI];

            const label pointI = f[0];

            label fLabel(-1);
            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                if( faces[faceI] == f )
                {
                    fLabel = faceI;
                    break;
                }
            }

            if( fLabel == -1 )
            {
                # ifdef USE_OMP
                //# pragma omp critical
                # endif
                {
                    faces.append(f);
                    c[fI++] = nFaces;

                    forAll(f, pI)
                        pointFaces.append(f[pI], nFaces);

                    ++nFaces;
                }
            }
            else
            {
                c[fI++] = fLabel;
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            cells.append(c);
            cellLevel.append(-1); //TODO: Get cellLevel from neighbouring cell?
            ++nCells;
        }
    }
*/
    // Ideas for OMP-friendly algorithm:
    // 1. Search to match existing faces
    // 2. For each new face added search in the other direction, from pointFaces
    // back to faces

    // Other idea:
    // As existing algo, but search pointFaces at the end for duplicates and eliminate

    // Other idea:
    // Do in batch. Search/add all/search/remove dups

    boolList pointLocked(pointFaces.size(), false);

    cells.setSize(nCells+cellFaces.size());
    cellLevel.setSize(nCells+cellFaces.size(), -1);

    # ifdef USE_OMP
    //# pragma omp parallel for schedule(dynamic, 50)
    # pragma omp parallel
    # pragma omp single
    # pragma omp taskloop default(shared) num_tasks(1)
    # endif
    forAll(cellFaces, cI)
    {
        faceList facesInCell(cellFaces.sizeOfGraph(cI));
        forAll(facesInCell, fI)
        {
            facesInCell[fI].setSize(cellFaces.sizeOfRow(cI, fI));

            forAll(facesInCell[fI], pI)
                facesInCell[fI][pI] = cellFaces(cI, fI, pI);
        }

        label fI(0);
        cell c(facesInCell.size());

        forAll(facesInCell, faceI)
        {
            const face& f = facesInCell[faceI];

            const label pointI = f[0];

            while(true)
            {
                bool clear = true;

                // Lock out the involved points. No other thread can add to the
                // pointFaces list involving these points until we're done
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                {
                    forAll (f, pI)
                    {
                        if (pointLocked[f[pI]]) 
                            clear = false;
                    }
                    if (clear)
                    {
                        forAll (f, pI)
                            pointLocked[f[pI]] = true;
                    }
                }
                if (clear)
                    break;

                // Points are not available. Go do some other cells in the 
                // loop before coming back and checking again
                # pragma omp taskyield
            }

            label fLabel(-1);
            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                if( faces[faceI] == f )
                {
                    fLabel = faceI;
                    break;
                }
            }

            if( fLabel == -1 )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                {
                    faces.append(f);
                    c[fI++] = nFaces;

                    forAll(f, pI)
                    {
                        pointFaces.append(f[pI], nFaces);
//                        Info << f[pI] << " " << pointFaces[f[pI]] << endl;
                    }

                    ++nFaces;
                }
            }
            else
            {
                c[fI++] = fLabel;
            }

            forAll (f, pI)
                pointLocked[f[pI]] = false;
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            //Info << nCells << " " << cI << endl;
            //cells.append(c);
            //cellLevel.append(-1); //TODO: Get cellLevel from neighbouring cell?
            //++nCells;
            cells[nCells+cI] = c;
        }
    }

    this->clearOut();
    mesh_.clearOut();

    Info << "Finished adding cells to the mesh" << endl;
}

void polyMeshGenModifier::addCell(const faceList& cellFaces)
{
    faceListPMG& faces = this->facesAccess();
    cellListPMG& cells = this->cellsAccess();
    labelIOList& cellLevel = mesh_.cellLevel_;

    label nFaces = faces.size();

    VRWGraph& pointFaces = this->pointFaces();

    cell c(cellFaces.size());
    label fI(0);

    forAll(cellFaces, faceI)
    {
        const face& f = cellFaces[faceI];

        const label pointI = f[0];

        label fLabel(-1);
        forAllRow(pointFaces, pointI, pfI)
        {
            const label faceI = pointFaces(pointI, pfI);

            if( faces[faceI] == f )
            {
                fLabel = faceI;
                break;
            }
        }

        if( fLabel == -1 )
        {
            faces.append(f);
            c[fI++] = nFaces;

            forAll(f, pI)
                pointFaces.append(f[pI], nFaces);

            ++nFaces;
        }
        else
        {
            c[fI++] = fLabel;
        }
    }

    cells.append(c);
    cellLevel.append(-1); //TODO: Get cellLevel from adjacent cell?
    mesh_.clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
