#include <array>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSTLReader.h>
#include "vtkSmartPointer.h"
#include "vtkIdTypeArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include <vtkXMLPolyDataWriter.h>
#include <vtkVector.h>
#include <vtkOBJReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <time.h>
#include <queue>
#include <set>
#include <fstream>
#include <algorithm>
#include "blossom5-v2.05.src/PerfectMatching.h"
#include <vtkMutableUndirectedGraph.h>
#include <vtkBoostConnectedComponents.h>

using namespace std;

const char* minor_name="Local Coord Minor Curvature Direction";
const char* euc_minor_name="Minor Curvature Direction";
const char* curvatures_name="Curvatures";
const double INV_ROOT_2 = 1.0/sqrt(2);


vtkVector3d operator+(vtkVector3d u, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=u[foo]+v[foo];
  }
  return w;
}
vtkVector3d operator-(vtkVector3d u, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=u[foo]-v[foo];
  }
  return w;
}
vtkVector3d operator*(double scalar, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=scalar*v[foo];
  }
  return w;
}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//Some topological operations
vector<int> ngb_triangles( int center_tri, vtkSmartPointer<vtkPolyData> surface){
  //vtk needs a good way to call the adjacent cells
  vector<int> ngbs;
  auto cell_pts = vtkSmartPointer<vtkIdList>::New();
  surface->GetCellPoints( center_tri , cell_pts);
  for(int bar=0; bar<3; bar ++){
    //Find neighboring cell across edge bar
    auto cell_ngb_wrap = vtkSmartPointer<vtkIdList>::New();
    surface->GetCellEdgeNeighbors( center_tri ,
      cell_pts->GetId( (bar+2)%3 ),
      cell_pts->GetId( (bar+1)%3 ),
      cell_ngb_wrap);
    if( cell_ngb_wrap->GetNumberOfIds()>0 ){
      int cell_ngb = cell_ngb_wrap->GetId(0);
      ngbs.push_back(cell_ngb);
    }
  }
  return ngbs;
}

set<int> adjacent_vertices(int seed, vtkSmartPointer<vtkPolyData> surface ){
  //Returns the adjacent ptIds of seed in the surface considered as a 2dim simplicial complex
  set<int> adj_verts;
  auto cellIdList =  vtkSmartPointer<vtkIdList>::New();	//get all connected cells
	surface->GetPointCells(seed, cellIdList);
	for(int i = 0; i < cellIdList->GetNumberOfIds(); i++){		// loop  through each cell using the seed point
	  vtkCell* cell = surface->GetCell(cellIdList->GetId(i));		// get current  cell
    auto cell_pts = cell->GetPointIds();
    for( int foo =0; foo< cell_pts->GetNumberOfIds(); foo++){
      int a_pt = cell_pts->GetId(foo);
      if( a_pt!= seed ) adj_verts.insert(a_pt);
    }
  }
  return adj_verts;
}

vector<int> cell_pair_to_vert_pair( int cell1, int cell2, vtkSmartPointer<vtkPolyData> surface){
  //Takes a pair of cells and returns the ptsIds of points they have in common
  set<int> pts1;
  vector<int> common_pts;
  auto cellIdList =  vtkSmartPointer<vtkIdList>::New();	//get all connected cells
	surface->GetCellPoints(cell1, cellIdList);
  for(int i = 0; i < cellIdList->GetNumberOfIds(); i++) pts1.insert(cellIdList->GetId(i));
  surface->GetCellPoints(cell2, cellIdList);
  for(int i = 0; i < cellIdList->GetNumberOfIds(); i++){
    int a_pt = cellIdList->GetId(i);
    if( pts1.count(a_pt) ) common_pts.push_back(a_pt);
  }
  if( common_pts.size() != 2 ) cout<< "cell_pair_to_vert_pair did not find a pair";
  sort(common_pts.begin(), common_pts.end());
  return common_pts;
}

double dot_ngbs(int tri1, int tri2, vtkSmartPointer<vtkPolyData> surface, int rot2){
  //Dot the fields of neighboring triangles. Rotate the 2nd vector by pi/2*rot2
  vtkSmartPointer<vtkDataArray> euc_field = surface->GetCellData()->GetArray(euc_minor_name);

  double field1[3];
  euc_field->GetTuple( tri1, field1);
  vtkVector3d v1(field1);

  double field2[3];
  euc_field->GetTuple( tri2, field2);
  vtkVector3d v2(field2);

  if( rot2%2 ==1 ){
    vtkSmartPointer<vtkDataArray> normals=surface->GetCellData()->GetNormals();
    double nortemp[3];
    normals->GetTuple( tri2, nortemp);
    vtkVector3d vn(nortemp);
    v2 = vn.Cross(v2);
  }
  double dotter = v1.Dot(v2);
  if( rot2 == 2 || rot2 == 3) dotter = - dotter;
  return dotter;
}

int integral_edge_weight(int vert1, int vert2, vtkSmartPointer<vtkPolyData> surface){
  //Option
  int bins = 1000;
  //
  vtkSmartPointer<vtkIdList> cellIdList =  vtkSmartPointer<vtkIdList>::New();
	surface->GetPointCells(vert1, cellIdList);
  vector<int> cells1;
  cells1.resize(cellIdList->GetNumberOfIds());
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++){		// loop  through each cell using the seed point
    cells1[i]=cellIdList->GetId(i);
  }

  surface->GetPointCells(vert2, cellIdList);
  vector<int> cells2;
  cells2.resize(cellIdList->GetNumberOfIds());
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++){		// loop  through each cell using the seed point
    cells2[i]=cellIdList->GetId(i);
  }
  vector<int> common_cells;
  sort(cells1.begin(),cells1.end());
  sort(cells2.begin(),cells2.end());
  set_intersection(cells1.begin(), cells1.end(), cells2.begin(), cells2.end(), inserter(common_cells, common_cells.begin()));
  if( common_cells.size() == 2){
    double w0 = dot_ngbs( common_cells[0], common_cells[1], surface, 0);
    double w1 = dot_ngbs( common_cells[0], common_cells[1], surface, 1);
    double w3 = dot_ngbs( common_cells[1], common_cells[0], surface, 1);
    double wght = w0*w0 - 0.5*(w1*w1 + w3*w3);
    return round( wght * bins);
  }
  else return 0;
}

vtkSmartPointer<vtkPolyData> build_small_simplicial_example()
{
  string example_name ("cylinder"); //"simple" "simple2" "cylinder" "grid"
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto triangles = vtkSmartPointer<vtkCellArray>::New();

  if (example_name.compare("simple") ==0 ){
    //Choose some points in parametric space
    double parametric_pts[6][2] = {
      {0.0,0.0},
      {0.0,1.0},
      {0.0,2.0},
      {1.0,0.0},
      {1.0,1.0},
      {2.0,0.0}
    };
    //Apply a ellipsoid parametrization to obtain points in R3
    int num_pts = sizeof(parametric_pts)/sizeof(parametric_pts[0]);
    points->SetNumberOfPoints(num_pts);
    for (int foo = 0; foo<num_pts; foo++){
      double pt[3];
      double p1 = .1 + .1 * parametric_pts[foo][0];
      double p2 = .2 + .2 * parametric_pts[foo][1];
      pt[0] = 5*cos(p1)*cos(p2);
      pt[1] = 4*sin(p1)*cos(p2);
      pt[2] = 3*sin(p2);
      points->SetPoint(foo, pt[0], pt[1], pt[2]);
    }
    int cellList[4][3] = {
      {0,1,3},
      {1,3,4},
      {1,2,4},
      {3,4,5}
    };
    for(size_t foo = 0; foo<sizeof(cellList)/sizeof(cellList[0]); foo++){
      triangles->InsertNextCell(3);
      for(size_t bar = 0; bar<sizeof(cellList[0])/sizeof(cellList[0][0]); bar++){
        triangles->InsertCellPoint(cellList[foo][bar]);
      }
    }
  }

  if (example_name.compare("simple2") == 0)
  {
    double angle0 = 0.0, angle1 = M_PI/6 , angle2 = 0.0;
    double pt_pos[6][3] = {
      {0,0,0},
      {1,0,0},
      {0,1,0},
      {0,-1,tan(angle0)},
      {-1,0,tan(angle1)},
      {1,1,tan(angle2)/sqrt(2)}
    };
    int num_pts = sizeof(pt_pos)/sizeof(pt_pos[0]);
    points->SetNumberOfPoints(num_pts);
    for (vtkIdType foo = 0; foo < num_pts; foo++) {
      double p0 = pt_pos[foo][0], p1 = pt_pos[foo][1], p2 = pt_pos[foo][2];
      points->SetPoint(foo, p0, p1, p2);
    }
    int cellList[4][3] = {
      {0,1,2},
      {0,3,1},
      {0,2,4},
      {1,5,2}
    };
    for(size_t foo = 0; foo<sizeof(cellList)/sizeof(cellList[0]); foo++){
      triangles->InsertNextCell(sizeof(cellList[0])/sizeof(cellList[0][0]));
      for(size_t bar = 0; bar<sizeof(cellList[0])/sizeof(cellList[0][0]); bar++){
        triangles->InsertCellPoint(cellList[foo][bar]);
      }
    }
  }
  if (example_name.compare("cylinder") == 0)
  {
    double cyl_wobble_amp = .3, cyl_wobble_freq=2 , cylrad=1, cylheight=1;
    //Suggest parameters .3,2,1
    int xtris = 70;
    int ytris = 70;
    int const num_pts = xtris*(ytris+1);
    vector< array<double,2> > param_pts;
    param_pts.resize(num_pts);
    int const num_tris = 2 * xtris * ytris;
    vector< array<int,3> > cell_list;
    cell_list.resize(num_tris);
    for(int foo = 0 ; foo < num_tris; foo++){
      int const sq = floor(foo/2);
      int const is_down = foo%2;
      int const pt_col = sq % xtris;
      int const pt_row = floor(sq / xtris);
      if(not is_down){
        cell_list[foo][0]= pt_col  + xtris * pt_row;
        cell_list[foo][1]=(1 + pt_col) % xtris  + xtris * (pt_row);
        cell_list[foo][2]=0 + pt_col + xtris * (pt_row+1);
      }
      else{
        cell_list[foo][0]=(1 + pt_col) % xtris + xtris * pt_row;
        cell_list[foo][1]=(1 + pt_col) % xtris + xtris * (pt_row+1);
        cell_list[foo][2]=0 + pt_col  + xtris * (pt_row+1);
      }
    }
    double const y_shift=sqrt(3)/2;
    //double const x_shift=1;
    for(int foo =0; foo < (xtris)*(ytris+1); foo++){
      int pt_col = foo%xtris;
      int pt_row = foo/xtris;
      param_pts[foo][0]=pt_col + 0.5*pt_row;
      param_pts[foo][1]=y_shift*pt_row;
    }
    //Apply parametrization to obtain points in R3
    //int num_pts = sizeof(param_pts)/sizeof(param_pts[0]);
    points->SetNumberOfPoints(num_pts);
    for (int foo = 0; foo<num_pts; foo++){
      double pt[3];
      double p1 = param_pts[foo][0];
      double p2 = param_pts[foo][1];
      double rad_now = cyl_wobble_amp * cos(2*M_PI*cyl_wobble_freq*p2/ytris ) + cylrad;
      pt[0] = rad_now * cos(2*M_PI*p1/xtris);
      pt[1] = rad_now * sin(2*M_PI*p1/xtris);
      pt[2] = 2*M_PI*cylheight*p2/ytris;
      points->SetPoint(foo, pt[0], pt[1], pt[2]);
    }
    for(int foo = 0; foo< num_tris; foo++){
      triangles->InsertNextCell(3);
      for(int bar = 0; bar<3; bar++){
        triangles->InsertCellPoint(cell_list[foo][bar]);
      }
    }
  }
  if (example_name == "grid")
  {
    int gridx=30;
    int gridy=10;
    int num_pts=(gridx+1)*(gridy+1);
    vector< array< double, 3 >> pt_pos;
    pt_pos.resize(num_pts);

    for(int foo=0; foo<num_pts; foo++){
      pt_pos[foo][0] = (foo % (gridx+1));
      pt_pos[foo][1] = floor(foo / (gridx+1) );
      pt_pos[foo][2] = 0;
    }
    points->SetNumberOfPoints(num_pts);
    for (int foo = 0; foo < num_pts; foo++) {
      double p0 = pt_pos[foo][0], p1 = pt_pos[foo][1], p2 = pt_pos[foo][2];
      points->SetPoint(foo, p0, p1, p2);
    }
    int num_triang = 2*gridx*gridy;
    vector< array<int,3> > cellList;
    cellList.resize(num_triang);

    for(int foo=0; foo<num_triang; foo++){
      int is_up = (foo+1) % 2;
      int xplace = ( foo / 2 ) % (gridx);
      int yplace = ( ( foo/ 2 ) / (gridx));
      if(is_up){
        cellList[foo][0] = xplace + (gridx+1) * yplace;
        cellList[foo][1] = 1+xplace + (gridx+1) * yplace;
        cellList[foo][2] = xplace + (gridx+1) * (1+yplace);
      }
      else{
        cellList[foo][1] = 1+xplace + (gridx+1) * (1+yplace);
        cellList[foo][0] = 1+xplace + (gridx+1) * yplace;
        cellList[foo][2] = xplace + (gridx+1) * (1+yplace);
      }
    }
    for(int foo = 0; foo<num_triang; foo++){
      triangles->InsertNextCell(3);
      for(int bar = 0; bar<3; bar++){
        triangles->InsertCellPoint(cellList[foo][bar]);
      }
    }
  }

  auto simplicial_complex = vtkSmartPointer<vtkPolyData>::New();
  simplicial_complex->SetPoints(points);
  simplicial_complex->SetPolys(triangles);
  simplicial_complex->BuildCells();
  simplicial_complex->BuildLinks();
  return simplicial_complex;
}


//Define a structure for a priority queue
struct cellPrior{
  int id;
  double confidence;

  cellPrior(int const id_, double const confidence_) : id{id_}, confidence{confidence_} {};
};
bool operator<(const cellPrior& cell1, const cellPrior& cell2){
  return fabs(cell1.confidence) < fabs(cell2.confidence);
}

struct patchPrior{
  int id;
  double confidence;
  int mother;

  patchPrior(int const id_, double const confidence_, int const mother_) : id{id_}, confidence{confidence_}, mother{mother_} {};
};
bool operator<(const patchPrior& cell1, const patchPrior& cell2){
  return fabs(cell1.confidence) < fabs(cell2.confidence);
}




//Function to align the vector field consistently
vtkSmartPointer<vtkPolyData> allign_field(vtkSmartPointer<vtkPolyData> surface){
  //Remembers the order cells were decided in the allignment
  int const num_surf_cells = surface->GetNumberOfCells();
  auto cell_order = vtkSmartPointer<vtkIntArray>::New();
  cell_order->SetNumberOfComponents(1);
  cell_order->SetNumberOfTuples( num_surf_cells );
  cell_order->SetName("Cell Order");

  //Remembers the patch cells belonged to
  auto alligned_patch = vtkSmartPointer<vtkIntArray>::New();
  alligned_patch->SetNumberOfComponents(1);
  alligned_patch->SetNumberOfTuples( num_surf_cells );
  alligned_patch->SetName("Concurrence Patch");
  vector<double> cell_to_patch;
  cell_to_patch.resize(num_surf_cells);

  vtkSmartPointer<vtkDataArray> euc_field = surface->GetCellData()->GetArray(euc_minor_name);
  vtkSmartPointer<vtkDataArray> curvatures = surface->GetCellData()->GetArray(curvatures_name);
  vtkSmartPointer<vtkDataArray> normals=surface->GetCellData()->GetNormals();


  //Partition the surface into patches where the field aligns locally within 45 degrees
  vector<int> flip_scheme;
  flip_scheme.resize(num_surf_cells);

  int cell_counter = 0;
  int patch_count = 0;
  vector<bool> cell_decided;
  cell_decided.resize(num_surf_cells);
  priority_queue< cellPrior > cell_que;
  for(int foo=0; foo<num_surf_cells; foo++){
    cell_decided[foo] = false;
    flip_scheme[foo]=0;
  }
  //Choose a cell, grow a patch until no neighboring cells align well, then grow another patch.
  for(int foo=0; foo<num_surf_cells; foo++){
    if( !cell_decided[foo] ){
      patch_count++;
      cell_que.push( cellPrior{foo,1} );
      int num_cells_in_patch = 0;
      int last_cell_seen = 0;
      while(!cell_que.empty()){
        cellPrior current = cell_que.top();
        cell_que.pop();
        cell_counter++;
        cell_order->SetTuple1( current.id, cell_counter);
        if(!cell_decided[current.id]){
          if (current.confidence < 0) flip_scheme[current.id] = 1;
          alligned_patch->SetTuple1( current.id, patch_count);
          cell_to_patch[current.id] = patch_count;
          cell_decided[current.id] = true;
          num_cells_in_patch++;
          last_cell_seen = current.id;
          //Look for neighbors
          vector<int> ngbs = ngb_triangles(current.id , surface);
          for( int cell_ngb : ngbs ){

            double current_field[3];
            euc_field->GetTuple( current.id, current_field);
            vtkVector3d current_v(current_field);

            double ngb_field_euc[3];
            euc_field->GetTuple( cell_ngb, ngb_field_euc);
            vtkVector3d ngb_v(ngb_field_euc);

            //Prioritize based on the angle ~= alignment
            //and difference cuvature magnitudes ~= liklihood direction is correct

            double conf = current_v.Dot(ngb_v);
            if( fabs(conf) > 0.9 ){
              double curvas[2];
              curvatures->GetTuple( cell_ngb, curvas );
              cellPrior get_in_line{ cell_ngb, 0};
              conf *= fabs(curvas[0]) - fabs(curvas[1]);
              if(current.confidence < 0) get_in_line.confidence = - conf;
              else get_in_line.confidence = conf;
              cell_que.push(get_in_line);
            }
          }
        }
      }
      if( num_cells_in_patch ==1 ){
        patch_count--;
        alligned_patch->SetTuple1(last_cell_seen, 0);
        cell_to_patch[last_cell_seen]=0;
      }
    }
  }
  patch_count++;
  //Too many lonely 1 triangle patches.

  //Flip the cells that were antialigned to get consistency per patch
  for(int foo=0; foo< num_surf_cells; foo++){
    if ( flip_scheme[foo] == 1 ){
      double vv[3];
      euc_field->GetTuple( foo, vv);
      euc_field->SetTuple3( foo, -vv[0], -vv[1], -vv[2]);
    }
  }
  //Extract the relative alignment between patches
  //If the entry in the patch adjacency matrix is positive,
  //the field should rotate positively (about normal) by pi/2.
  //Similarly if negative.
  //Allocate
  vector< vector<double> > patch_adjmat;
  patch_adjmat.resize(patch_count);
  for(int foo=0; foo<patch_count; foo++){
    patch_adjmat[foo].resize(patch_count);
    for(int bar=0; bar<patch_count; bar++)patch_adjmat[foo][bar]=0;
  }

  for(int foo=0; foo< num_surf_cells; foo++){
    int foopatch = cell_to_patch[foo];
    patch_adjmat[foopatch][foopatch]+=1;
    vector<int> ngbs = ngb_triangles(foo, surface);
    for( int cell_ngb : ngbs){
      int ngbpatch = cell_to_patch[cell_ngb];
      if( foopatch != ngbpatch ){

        double tmpe[3];
        euc_field->GetTuple(foo,tmpe);
        vtkVector3d foofield(tmpe);

        double tmpb[3];
        euc_field->GetTuple(cell_ngb,tmpb);
        vtkVector3d barfield(tmpb);

        double tmpn[3];
        normals->GetTuple(cell_ngb,tmpn);
        vtkVector3d barnormal(tmpn);

        double align_weight = foofield.Dot( barnormal.Cross( barfield ));
        patch_adjmat[foopatch][ngbpatch] += align_weight;
      }
    }
  }
  //Find the largest patch
  int big_patch=0;
  int high_score=0;
  for(int foo=0; foo<patch_count; foo++){
    int p = patch_adjmat[foo][foo];
    if( p > high_score){
      big_patch=foo;
      high_score=p;
    }
  }

  //Decide how to rotate patches to allign moving out from the biggest patch.
  vector<int> rot_scheme;
  rot_scheme.resize(patch_count);
  rot_scheme[big_patch] = 0;
  priority_queue< patchPrior > patch_que;
  vector<bool> patch_decided;
  patch_decided.resize(patch_count);
  for(int foo=0; foo<patch_count; foo++) patch_decided[foo]=false;
  patch_que.push( patchPrior{big_patch,0,big_patch} );
  while(!patch_que.empty()){
    patchPrior adams = patch_que.top();
    patch_que.pop();
    if(!patch_decided[adams.id]){
      int mom = adams.mother;
      rot_scheme[adams.id] = (4+rot_scheme[mom] - sgn(adams.confidence))%4 ;
      patch_decided[adams.id]=true;
      for(int foo=0; foo<patch_count; foo++){
        double confy = patch_adjmat[foo][adams.id];
        if( confy != 0 ) patch_que.push(patchPrior{foo , confy, adams.id});
      }
    }
  }
  //Rotate vectors according to the scheme.
  // for(int foo=0; foo<num_surf_cells; foo++){
  //   int patch = cell_to_patch[foo];
  //   int rotate_to = rot_scheme[patch];
  //   if ( (rotate_to%2) == 1){
  //     double tmpe[3];
  //     euc_field->GetTuple(foo,tmpe);
  //     vtkVector3d foofield(tmpe);
  //     double tmpn[3];
  //     normals->GetTuple(foo,tmpn);
  //     vtkVector3d foonorm(tmpn);
  //     foofield = foonorm.Cross(foofield);
  //     if (rotate_to == 3) foofield = -1*foofield;
  //     euc_field->SetTuple3(foo,foofield[0],foofield[1],foofield[2]);
  //   }
  //   else if ( rotate_to == 2){
  //     double tmpe[3];
  //     euc_field->GetTuple(foo,tmpe);
  //     euc_field->SetTuple3(foo, -tmpe[0], -tmpe[1], -tmpe[2]);
  //   }
  // }

  //I want to look at the rotation scheme, write it to the cells
  auto rotation_scheme = vtkSmartPointer<vtkIntArray>::New();
  rotation_scheme->SetNumberOfComponents(1);
  rotation_scheme->SetNumberOfTuples(num_surf_cells);
  rotation_scheme->SetName("Rotation Scheme");
  for(int foo=0; foo<patch_count; foo++){
    if (rot_scheme[ cell_to_patch[foo] ] > 0 ){
      rotation_scheme->SetTuple1(foo, rot_scheme[ cell_to_patch[foo] ] );
    }
  }
  surface->GetCellData()->AddArray( rotation_scheme );
  surface->GetCellData()->AddArray( cell_order );
  surface->GetCellData()->AddArray( alligned_patch );
  surface->GetCellData()->AddArray( euc_field );
  return surface;
}




vtkSmartPointer<vtkPolyData> compute_curvature_frame(vtkSmartPointer<vtkPolyData> surface){
  //Compute the minor curvature direction and principle curvatures everywhere
  int num_surf_cells = surface->GetNumberOfCells();
  vtkDataArray *normals=surface->GetCellData()->GetNormals();
  //The minor curvature field here is expressed over the standar basis in RR3
  auto euc_field = vtkSmartPointer<vtkDoubleArray>::New();
  euc_field->SetNumberOfComponents(3);
  euc_field->SetNumberOfTuples(num_surf_cells);
  euc_field->SetName(euc_minor_name);
  //Also store the principle curvatures
  auto curvatures = vtkSmartPointer<vtkDoubleArray>::New();
  curvatures->SetNumberOfComponents(2);
  curvatures->SetNumberOfTuples(num_surf_cells);
  curvatures->SetName(curvatures_name);
  //Loop through the cells to compute the shape operator locally and extract
  //the minor curvature direction
  for(vtkIdType foo=0; foo< num_surf_cells; foo++){
    //Get the points as vectors
    double temp_pts[3][3];
    vtkVector3d pts[3];
    auto cell_pts = vtkSmartPointer<vtkIdList>::New();
    surface->GetCellPoints(foo,cell_pts);
    for(int bar=0;bar<3;bar++){
      surface->GetPoint(cell_pts->GetId(bar),temp_pts[bar]);
      pts[bar]=vtkVector3d(temp_pts[bar]);
    }
    //Get the normal vector
    double temp_n[3];
    normals->GetTuple(foo, temp_n);
    vtkVector3d normal(temp_n);
    //Compute_edge_vectors
    //For each edge compute the shape operator contribution.
    //The shape operator is discretized as
    // sum_{e in edges} (angle at edge) projection matrix on tangent perp to e
    // The shape operator is symmetric and computed in local edge coordinates.
    // Note that since we only care about direction, the shape operator is not
    // locally normalized so that the curvature is incorrect.
    double edge_lengths[3];
    vtkVector3d edges[3];// tang_vects[3];
    double angles[3];
    //Compute the triangle edges
    for(int bar=0; bar<3; bar ++){
      //Edge vectors for orientation 012
      edges[bar] = pts [ (2+bar)%3 ] - pts [ (1+bar)%3 ];
      //Edge lengths
      edge_lengths[bar]=edges[bar].Norm();
      //Find neighboring cell across edge bar
    }
    //Compute a local orthonormal basis
    vtkVector3d orthobase0(edges[0]), orthobase1;
    orthobase0.Normalize();
    orthobase1 = normal.Cross(orthobase0);
    double shape00=0, shape01=0, shape11=0;
    for(int bar=0; bar<3; bar ++){
      //Find neighboring cell across edge bar
      auto cell_ngb = vtkSmartPointer<vtkIdList>::New();
      surface->GetCellEdgeNeighbors(foo,
        cell_pts->GetId( (2+bar)%3 ),
        cell_pts->GetId( (1+bar)%3 ),
        cell_ngb);
      //If there is no neighboring cell, i.e. the edge is surface boundary,
      //then we assume the surface is simply flat in that direction.
      //Otherwise add the angle scaled projection in that direction.
      if(cell_ngb->GetNumberOfIds()>0){
        normals->GetTuple(cell_ngb->GetId(0),temp_n);
        vtkVector3d ngb_normal(temp_n);
        //The sign of the angle is taken to agree with right-hand rotation
        //about the edge vector.
        double sine_of_ang = edges[bar].Dot( normal.Cross(ngb_normal) ) / edge_lengths[bar];
        if (sine_of_ang > 1) angles[bar] = M_PI/2;
        else if (sine_of_ang < -1) angles[bar] = -M_PI/2;
        else angles[bar] = asin( sine_of_ang ) ;
        angles[bar] = angles[bar]/edge_lengths[bar];
        //Compute the shape operator matrix entries shape10=shape01
        //but except the eigenvalues are associated with the opposite eigenvectors
        //as in ~H of https://graphics.stanford.edu/courses/cs468-03-fall/Papers/cohen_normalcycle.pdf
        double dot0bar = orthobase0.Dot(edges[bar]);
        double dot1bar = orthobase1.Dot(edges[bar]);
        shape00 += angles[bar] * dot0bar * dot0bar;
        shape01 += angles[bar] * dot0bar * dot1bar;
        shape11 += angles[bar] * dot1bar * dot1bar;
      }
    }
  //Compute the minor eigenvalue and the corresponding eigenvector
    double shape_trace=shape00+shape11;
    double shape_discr=shape_trace*shape_trace-4*(shape00*shape11-shape01*shape01);
    if (shape_discr<0) cout<<"Contradiction! Shape operator eigs are complex!" <<foo;
    double major_curvature, minor_curvature;
    if(shape_trace>0){
      major_curvature=(shape_trace + sqrt(shape_discr))/2.0;
      minor_curvature=(shape_trace - sqrt(shape_discr))/2.0;
    }
    else{
      major_curvature=(shape_trace - sqrt(shape_discr))/2.0;
      minor_curvature=(shape_trace + sqrt(shape_discr))/2.0;
    }
    curvatures->SetTuple2(foo, major_curvature, minor_curvature);
    double ortho_curv_vect[2];
    if (fabs(major_curvature-shape00) > fabs(major_curvature-shape11)){
      ortho_curv_vect[0] = shape01;
      ortho_curv_vect[1] = major_curvature-shape00;
    }
    else{
      ortho_curv_vect[0] = major_curvature-shape11;
      ortho_curv_vect[1] = shape01;
    }
    vtkVector3d euc_curv_vect =  ortho_curv_vect[0] * orthobase0 + ortho_curv_vect[1] * orthobase1 ;
    double magni = euc_curv_vect.Norm();
    if (magni>0) euc_curv_vect =  1/magni * euc_curv_vect ;
    //Hacky attempt at semi continuity by choosing rando hemisphere to rep RealProjective2Space
    //Pre-pre-alignment
    if(euc_curv_vect[2] < 0){
      euc_curv_vect = -1 * euc_curv_vect;
    }
    euc_field->SetTuple3( foo , euc_curv_vect[0], euc_curv_vect[1], euc_curv_vect[2]);
  }

  //Pre-alignment
  //Decide a dominant direction to roughly align the field with.

  vtkVector3d dominant_dir={0,0,0};
  for(int foo=0; foo<num_surf_cells; foo++){
    double tempf[3];
    euc_field->GetTuple( foo , tempf);
    vtkVector3d tempv(tempf);
    dominant_dir = dominant_dir + tempv;
  }
  cout << "Dominant direction " << dominant_dir[0] <<" " << dominant_dir[1] <<" " << dominant_dir[2] << endl;
  for(int foo=0; foo<num_surf_cells; foo++){
    double tempf[3];
    euc_field->GetTuple( foo , tempf);
    vtkVector3d tempv(tempf);
    if( dominant_dir.Dot(tempv) < 0){
      tempv= -1 * tempv;
      euc_field->SetTuple3( foo , tempv[0], tempv[1], tempv[2]);
    }
  }

  //Looks noisey. Try a uniform window blur on cell neighbors.
  //There should be a fast way to do this with a convolution filter....
  // bool apply_blur = false;
  // if(apply_blur){
  //   const char* blur_name="BlurredField";
  //   auto blurred_field = vtkSmartPointer<vtkDoubleArray>::New();
  //   blurred_field->SetNumberOfComponents(3);
  //   blurred_field->SetNumberOfTuples(num_surf_cells);
  //   blurred_field->SetName(blur_name);
  //   for(int foo = 0; foo<num_surf_cells; foo++){
  //     double tempf[3];
  //     euc_field->GetTuple( foo , tempf);
  //     vtkVector3d field_here(tempf);
  //     int edge_ends[3][3] = {{2,1},{0,2},{1,0}};
  //     auto cell_pts = vtkSmartPointer<vtkIdList>::New();
  //     surface->GetCellPoints(foo,cell_pts);
  //     auto cell_ngb = vtkSmartPointer<vtkIdList>::New();
  //     for(int bar=0;bar<3;bar++){
  //       surface->GetCellEdgeNeighbors(foo,
  //       cell_pts->GetId(edge_ends[bar][0]),
  //       cell_pts->GetId(edge_ends[bar][0]),
  //       cell_ngb);
  //       if(cell_ngb->GetNumberOfIds()>0){
  //         euc_field->GetTuple(cell_ngb->GetId(0), tempf);
  //         vtkVector3d field_there(tempf);
  //         field_here=field_here+field_there;
  //       }
  //     }
  //     field_here.Normalize();
  //     double * w = &field_here[0];
  //     blurred_field->SetTuple(foo, w);
  //     if (surface->GetCellData()->HasArray(blur_name))
  //       surface->GetCellData()->RemoveArray(blur_name);
  //     surface->GetCellData()->AddArray(blurred_field);
  //   }
  // }
  // If an old curvature array already exists we remove it
  if (surface->GetCellData()->HasArray(euc_minor_name)){
    surface->GetCellData()->RemoveArray(euc_minor_name);
  }
  if (surface->GetCellData()->HasArray(curvatures_name)){
    surface->GetCellData()->RemoveArray(curvatures_name);
  }
  surface->GetCellData()->AddArray(euc_field);
  surface->GetCellData()->AddArray(curvatures);

  surface = allign_field(surface);


  //Initialize memory for the minor curvautre field.
  //The field is defined per triangle of the surface.
  //The field is the direction of the shaper operator eignevalue with the lower
  //magintude curvature.
  auto minor_curv_field = vtkSmartPointer<vtkDoubleArray>::New();
  minor_curv_field->SetNumberOfComponents(2);
  minor_curv_field->SetNumberOfTuples(num_surf_cells);
  minor_curv_field->SetName(minor_name);

  //Compute the local coordinate representation over
  // local edge basis e0=p2-p1 and e1=p0-p2
  for(int foo=0; foo<num_surf_cells; foo++){
    //Get the points as vectors
    double temp_pts[3][3];
    vtkVector3d pts[3];
    auto cell_pts = vtkSmartPointer<vtkIdList>::New();
    surface->GetCellPoints( foo, cell_pts);
    for(int bar=0; bar<3; bar++){
      surface->GetPoint( cell_pts->GetId(bar), temp_pts[bar] );
      pts[bar] = vtkVector3d(temp_pts[bar]);
    }
    vtkVector3d edges[3];
    //Compute the triangle edges
    for(int bar=0; bar<3; bar ++){
      //Edge vectors for orientation 012
      edges[bar] = pts [ (2+bar)%3 ] - pts [ (1+bar)%3 ];
    }
    double e00 = edges[0].Dot(edges[0]), e01 = edges[0].Dot(edges[1]), e11 = edges[1].Dot(edges[1]);
    double edet = e00*e11-e01*e01;
    double temp[3];
    euc_field->GetTuple(foo,temp);
    vtkVector3d euc_curv_vect(temp);
    double e0dw = edges[0].Dot(euc_curv_vect)  , e1dw = edges[1].Dot(euc_curv_vect) ;
    double minor_curv_vect[2] ={
      (e11*e0dw - e01 * e1dw) /edet,
      (e00*e1dw - e01 * e0dw) /edet,
    };
    minor_curv_field->SetTuple( foo , minor_curv_vect);
  }
  if (surface->GetCellData()->HasArray(minor_name)){
    surface->GetCellData()->RemoveArray(minor_name);
  }
  surface->GetCellData()->AddArray(minor_curv_field);

  return surface;
}


vtkSmartPointer<vtkPolyData> compute_tubular_parametrization(vtkSmartPointer<vtkPolyData> surface){
  //Contstruct a function u: {mesh points}->RR with gradient ~= minor curvature field
  int num_pts = surface->GetNumberOfPoints();
  auto cell_pts = vtkSmartPointer<vtkIdList>::New();
  int num_cells = surface->GetNumberOfCells();
  //Construct the discrete gradient matrix for the surface
  //The size is slightly larger to get constraint u(0)=0
  vector< Eigen::Triplet<double> > grad_entrylist;
  grad_entrylist.reserve(4 * num_cells+1);
  vector< Eigen::Triplet<double> > curv_entrylist;
  curv_entrylist.reserve(2 * num_cells+1);
  double curv_entry[2];
  vtkSmartPointer<vtkDataArray> curv_data = surface->GetCellData()->GetArray(minor_name);
  //Porting data from vtkArray to a Eigen::Matrix
  //This should probably be done when you compute it in the first place
  for(int foo=0; foo< num_cells; foo++){
    surface->GetCellPoints(foo,cell_pts);
    int pt0=cell_pts->GetId(0), pt1=cell_pts->GetId(1), pt2=cell_pts->GetId(2);
    //cout<<foo<<"cell's pts"<<pt0<<pt1<<pt2<<endl;
    curv_data->GetTuple(foo,curv_entry);
    double veccomp0 = curv_entry[0];
    double veccomp1 = curv_entry[1];
    double vec_norm = sqrt(veccomp0 * veccomp0 + veccomp1 * veccomp1);
    //This threshold for normalization is totally arbitrary
    //and should be informed by stats on the curvature.
    double normalizing_thresh = 1e-5;
    if (vec_norm > normalizing_thresh){
      veccomp0 *= 1/vec_norm;
      veccomp1 *= 1/vec_norm;
    }
    curv_entrylist.push_back( Eigen::Triplet<double> ( 2*foo, 0, veccomp0 ) );
    curv_entrylist.push_back( Eigen::Triplet<double> ( 2*foo+1, 0, veccomp1 ) );
    //The graditent per cell has component0 = u(p2)-u(p1)
    //and component1 = u(p0)-u(p2)
    grad_entrylist.push_back( Eigen::Triplet<double> (2*foo, pt2, 1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (2*foo, pt1, -1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (2*foo+1, pt0, 1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (2*foo+1, pt2, -1) );
  }
  //Add constraint u(0)=0
  grad_entrylist.push_back( Eigen::Triplet<double> (2*num_cells, 0, 1) );
  curv_entrylist.push_back( Eigen::Triplet<double> (2*num_cells, 0, 0) );

  //Build Matrix
  Eigen::SparseMatrix<double> surf_gradient( 2*num_cells + 1, num_pts);
  surf_gradient.setFromTriplets(grad_entrylist.begin(), grad_entrylist.end());
  //Build target vector i.e. the curvature field
  Eigen::SparseMatrix<double> curv_field( 2*num_cells + 1, 1);
  curv_field.setFromTriplets(curv_entrylist.begin(), curv_entrylist.end());

  Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
  cout<<"Factoring surface gradient."<<endl;
  solver.compute(surf_gradient);
  Eigen::SparseMatrix<double> u_param(num_pts,1);
  cout<<"Constructing tubular potential function."<<endl;
  time_t time1,time2;
  time1 = time(0);
  u_param = solver.solve(curv_field);
  time2 = time(0);
  int it_took = difftime(time2,time1);
  cout<<"Integration finished in " << it_took/60 << " min and " << it_took%60 << " sec." <<endl;
    cout<<"Required " << solver.iterations() << " iterations obtaining error " << solver.error() <<endl;
  const char* tube_param_name="Tubular Parametrization";
  auto tube_param = vtkSmartPointer<vtkDoubleArray>::New();
  tube_param->SetNumberOfComponents(1);
  tube_param->SetNumberOfTuples(num_pts);
  tube_param->SetName(tube_param_name);
  for(int foo=0; foo<num_pts; foo++){
    tube_param->SetTuple1(foo, u_param.coeffRef(foo,0));
  }
  surface->GetPointData()->AddArray(tube_param);
  return surface;
}

int main(int argc, const char* argv[]){
  const char* output_name;
  auto surface= vtkSmartPointer<vtkPolyData>::New();
  bool is_input = true;
  if (is_input) {
    if (argc < 3) {
      fprintf(stderr,"Usage: %s <input> [<output.vtp>] \n",argv[0]);
      return 0;
    }

    int length = strlen(argv[1]);
    const char* ext = argv[1] + length - 3;
    bool is_vtp = false;

    if (strcmp(ext,"obj") == 0) {
      vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
      reader->SetFileName(argv[1]);
      reader->Update();
      surface = reader->GetOutput();
    }
    else if (strcmp(ext,"stl") == 0) {
      vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
      reader->SetFileName(argv[1]);
      reader->Update();
      surface = reader->GetOutput();
    }
    else if (strcmp(ext,"vtp") == 0) {
      vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
      reader->SetFileName(argv[1]);
      reader->Update();
      surface = reader->GetOutput();

      is_vtp = true;
    }
    else {
      fprintf(stderr,"Unkown file format %s\n",ext);
      return 0;
    }

    // If we have no output name but the input file is not a vtp
    if ((argc < 3) && !is_vtp) {
      fprintf(stderr,"Error: For non vtp input files an output name is required\n");
      exit(0);
    }

    if (argc > 2)
      output_name = argv[2];
    else
      output_name = argv[1];
  }
  else{
    output_name = "bongo.vtp";
    surface = build_small_simplicial_example();
  }

  // auto surface = vtkSmartPointer<vtkPolyData>::New();
  // vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  // reader->SetFileName(argv[0]);
  // reader->Update();
  // surface = reader->GetOutput();
  //  surface->BuildCells();
  auto normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

  //Compute normal vectors for all cells
  cout<<"Computing normals."<<endl;
  normalGenerator->SetInputData(surface);
  normalGenerator->ComputePointNormalsOff();
  normalGenerator->ComputeCellNormalsOn();
  normalGenerator->SetSplitting(0);
  normalGenerator->AutoOrientNormalsOn();
  normalGenerator->Update();
  surface = normalGenerator->GetOutput();
  //surface->DeleteCells();
  surface->BuildLinks();
  cout<<"Computing minor principle curvature field. "<<endl;

  surface = compute_curvature_frame(surface);
  // auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  // sphereSource->Update();
  // vtkSmartPointer<vtkPolyData> sweet_sphere = sphereSource->GetOutput();
  //
  // auto polyDataToGraphFilter = vtkSmartPointer<vtkPolyDataToGraph>::New();
  // polyDataToGraphFilter->SetInputData(sweet_sphere);
  // polyDataToGraphFilter->Update();
  //
  // // auto graph_on_a_sphere = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  // vtkSmartPointer<vtkGraph> graph_on_a_sphere = polyDataToGraphFilter->GetOutput();
  // cout<< "This graph lives on a sphere!" <<endl;
  // int vv = graph_on_a_sphere->GetNumberOfVertices();
  // int ee = graph_on_a_sphere->GetNumberOfEdges();
  // int ff = sweet_sphere->GetNumberOfCells();
  // sweet_sphere->BuildLinks();
  // cout<< vv << " vertices.\n";
  // cout<< ee << " edges.\n";
  // cout<< ff << " faces.\n";
  // cout<< "Of course " << vv << " - " << ee << " + " << ff << " = " << vv-ee+ff<<endl;
  // cout<< "So there must be " << 2-(vv-ee+ff)<< " boundary components." <<endl;

  //Build a list of the edges of the simplicial surface to have integer identifies for the blowup
  //In the original graph we direct edges (i,j) to have increasing vertex id so i<j
  //The reversed edge (j,i) has edge_id(j,i) = edge_id(i,j)+num_edges
  //So an easy test to see if two edges are the reversals of each other is if
  // e1 % num_edges == e2 % num_edges
  int num_surf_cells = surface->GetNumberOfCells();
  int num_surf_verts = surface->GetNumberOfPoints();
  int blwp_edge=0;
  vector< array<int,2> > orig_edge_list;
  orig_edge_list.reserve( 2*num_surf_cells );
  map< array<int,2> , int > orig_edge_lookup;
  int num_surf_edges=0;
  for(int foo=0; foo< num_surf_verts; foo++){
    set<int> adj_pts = adjacent_vertices(foo, surface);
    for( int bar : adj_pts){
      if( foo<bar ){
        orig_edge_list.push_back( {foo, bar} );
        orig_edge_lookup.insert( { {foo,bar} , num_surf_edges} );
        num_surf_edges++;
      }
    }
  }

  //This function allows lookups of the original edge number from the blowup
  auto find_edge_id = [&orig_edge_lookup, &num_surf_edges]( int v1, int v2){
    int x;
    if (v1 < v2) x = orig_edge_lookup[ {v1, v2} ];
    if (v1 > v2) x = orig_edge_lookup[ {v2, v1} ] +num_surf_edges;
    return x;
  };

  //Build the blowup graph by replacing every vectex v with a deg(v) sized clique
  //There is a bijection between directed edges of G and the vertices of blowup(G)
  vector< array<int,3> > blwp_edge_list;
  blwp_edge_list.reserve( 30 *surface->GetNumberOfCells() );
  for(int foo=0; foo< num_surf_verts; foo++){
    set<int> adj_pts = adjacent_vertices(foo, surface);
    for( int bar : adj_pts){
      //Add edges for the original edges in the graph
      if( foo<bar ){
        int wght = integral_edge_weight(foo, bar, surface);
        int oredge =  orig_edge_lookup[ {foo,bar} ];
        blwp_edge_list.push_back( { oredge, oredge + num_surf_edges, wght} );
      }
      //Add edges for all the cliques
      for( int qux : adj_pts){
        if (bar < qux){
          int dir_ed_1 = find_edge_id(foo,bar);
          int dir_ed_2 = find_edge_id(foo,qux);
          blwp_edge_list.push_back( { dir_ed_1, dir_ed_2, 0});
        }
      }
    }
  }

  // //I need a better way to interact with blossomV but for now I'll have to
  // //write to a txt files to communicate between the two
  // ofstream myfile;
  // myfile.open("blowup.txt");
  // myfile << 2*num_surf_edges << " " << blwp_edge_list.size()  <<endl;
  // for( auto an_edge : blwp_edge_list){
  //   myfile << an_edge[0] << " " << an_edge[1] << " " << an_edge[2] << endl;
  // }
  // myfile.close();

  //maybe we dont have to match with integers, look into that
  PerfectMatching blossom_matching( 2*num_surf_edges , blwp_edge_list.size() );
  for (auto an_edge: blwp_edge_list){
    blossom_matching.AddEdge(an_edge[0], an_edge[1], an_edge[2]);
  }
  blossom_matching.Solve();
  //blossom_matching.GetMatch()

  //Build the graph to find the cut corresponding to the blossomV matching
  auto flip_graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  flip_graph->SetNumberOfVertices(num_surf_cells);
  for(int foo=0; foo< num_surf_cells; foo++){
    vector<int> adj_tris = ngb_triangles(foo, surface);
    for(int bar : adj_tris){
      if( foo < bar ) {
        vector<int> ptpair = cell_pair_to_vert_pair( foo, bar, surface);
        int ed_id = find_edge_id( ptpair[0], ptpair[1]);
        int match_id = blossom_matching.GetMatch(ed_id);
        //the graph shouldn't have any of the edges chosen by the matching
        //add in the other edges
        if( ed_id % num_surf_edges != match_id % num_surf_edges ){
          flip_graph->AddEdge(foo,bar);
        }
      }
    }
  }

  auto connectedComponents = vtkSmartPointer<vtkBoostConnectedComponents>::New();
  connectedComponents->SetInputData( flip_graph);
  connectedComponents->Update();
  vtkGraph* outputGraph = connectedComponents->GetOutput();
  vtkIntArray* components = vtkIntArray::SafeDownCast(
    outputGraph->GetVertexData()->GetArray("component"));
  for(vtkIdType i = 0; i < components->GetNumberOfTuples(); i++)
  {
    int val = components->GetValue(i);
    // std::cout << val << std::endl;
  }

  surface->GetCellData()->AddArray(components);


  cout<<"Writing file "<<output_name<<endl;
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(output_name);
  writer->SetInputData(surface);
  writer->Write();
}
