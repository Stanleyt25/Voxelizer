
//Computational Fabrication Assignment #1
// By David Levin 2014
#include <iostream>
#include <vector>
#include <chrono>
using namespace std::chrono;
#include "../include/CompFab.h"
#include "../include/Mesh.h"

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray& ray, CompFab::Triangle& triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle,
     * 0 otherwise */

     // In order for the ray to intersect the triangle it must:
     // (1) intersect the plane in which the triangle resides
     // (2) be contained within the triangle

    // n is the normal of the triangle!
    CompFab::Vec3 n = CompFab::operator%(CompFab::operator-(triangle.m_v2, triangle.m_v1), CompFab::operator-(triangle.m_v3, triangle.m_v1));
    n.CompFab::Vec3Struct::normalize();

    // d is just the direction the ray is pointing in!
    CompFab::Vec3 d = ray.m_direction;
    d.CompFab::Vec3Struct::normalize();

    // The "denominator" is the dot product of the ray's direction and the triangle normal
    float denominator = CompFab::operator*(d, n);

    // First, we check if the dot product of the ray's direction and the triangle exists (if not, they're parallel and never intersect!)
    // To account for floating point error, if the dot product is within epsilon of 0, we'll count is as zero since it's so tiny
    if (denominator < EPSILON && denominator > -EPSILON) {
        return 0;
    }
    else {
        // t is the time at which the ray intersects the plane -- we will use this to calculate the point of intersection

        float t = CompFab::operator*(CompFab::operator-(triangle.m_v1, ray.m_origin), n) / denominator;

        if (t < 0) {
            return 0;
        }

        else {
            // x is the point of intersection between the ray and plane
            CompFab::Vec3 x((ray.m_origin[0] + (t * d[0])), (ray.m_origin[1] + (t * d[1])), (ray.m_origin[2] + (t * d[2])));

            // Here, we see if the point of intersection is "contained" by the three vertices of the triangle
            if (CompFab::operator*(CompFab::operator%(CompFab::operator-(triangle.m_v2, triangle.m_v1), CompFab::operator-(x, triangle.m_v1)), n) > 0 &&
                CompFab::operator*(CompFab::operator%(CompFab::operator-(triangle.m_v3, triangle.m_v2), CompFab::operator-(x, triangle.m_v2)), n) > 0 &&
                CompFab::operator*(CompFab::operator%(CompFab::operator-(triangle.m_v1, triangle.m_v3), CompFab::operator-(x, triangle.m_v3)), n) > 0) {

                // If so, we return 1, or true!
                return 1;
            }
            else {
                return 0;
            }
        }
    }
}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    CompFab::Ray ray(voxelPos, dir);
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */

    for (auto& tri: g_triangleList) {
        if (rayTriangleIntersection(ray, tri) == 1) {
            numHits++;
        }
    }
    
    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}


int main(int argc, char **argv)
{

    unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)

    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        exit(0);
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim);
    

    
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);
    
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    auto start = high_resolution_clock::now();

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                CompFab::Vec3 voxelPos((0.5f + (double)ii) * spacing, (0.5f + (double)jj) * spacing, (0.5f + (double)kk) * spacing);
                if (numSurfaceIntersections(voxelPos, direction) % 2 != 0) {
                    g_voxelGrid->isInside(ii, jj, kk) = 1;
                }
            }
        }
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Time taken to voxelize the object is: " << duration.count() << " milliseconds" << std::endl;

    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    
    delete g_voxelGrid;
}