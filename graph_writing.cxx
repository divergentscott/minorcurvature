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
#include <vtkBoostBreadthFirstSearch.h>
#include <vtkReebGraph.h>
#include <vtkGraphWriter.h>
#include <vtkGraphLayoutView.h>
#include <vtkSimple2DLayoutStrategy.h>
#include <vtkRenderWindowInteractor.h>

using namespace std;

const double M_PI = 3.14159265359;
const char* minor_name="Local Coord Minor Curvature Direction";
const char* euc_minor_name="Minor Curvature Direction";
const char* curvatures_name="Curvatures";
const char* tube_param_name="Tubular Parametrization";
const char* potential_name="Potential";
const double INV_ROOT_2 = 1.0/sqrt(2);

//RR3 vectors should have addition
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

//Define structures for a priority queue
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


//Some topological operations//
//--+++++++++++++++--//
vector<int> adjacent_triangles( int center_tri, vtkSmartPointer<vtkPolyData> surface){
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

vector<int> cell_pair_verts_shared( int cell1, int cell2, vtkSmartPointer<vtkPolyData> surface){
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
  if( common_pts.size() != 2 ) cout<< "cell_pair_verts_shared did not find a pair";
  sort(common_pts.begin(), common_pts.end());
  return common_pts;
}

vector<int> vert_pair_cells_shared( int pt1, int pt2, vtkSmartPointer<vtkPolyData> surface){
  //Takes a pair of adjacent vertices and returns the dual pair of cells
  set<int> cells1;
  vector<int> common_cells;
  auto cellIdList =  vtkSmartPointer<vtkIdList>::New();	//get all connected cells
	surface->GetPointCells(pt1, cellIdList);
  for(int i = 0; i < cellIdList->GetNumberOfIds(); i++) cells1.insert(cellIdList->GetId(i));
  surface->GetPointCells(pt2, cellIdList);
  for(int i = 0; i < cellIdList->GetNumberOfIds(); i++){
    int a_cell = cellIdList->GetId(i);
    if( cells1.count(a_cell) ) common_cells.push_back(a_cell);
  }
  sort(common_cells.begin(), common_cells.end());
  return common_cells;
}
//--+++++++++++++++--//



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

int perpin_weights(int vert1, int vert2, vtkSmartPointer<vtkPolyData> surface){
  //Computes the weights for a maxcut of parllel vs perpindicular vectors
  //Option
  int bins = 10000000;
  vector<int> common_cells =  vert_pair_cells_shared(vert1, vert2, surface);
  if( common_cells.size() == 2){
    double w0 = dot_ngbs( common_cells[0], common_cells[1], surface, 0);
    double w1 = dot_ngbs( common_cells[0], common_cells[1], surface, 1);
    double w3 = dot_ngbs( common_cells[1], common_cells[0], surface, 1);
    double wght = w0*w0 - 0.5*(w1*w1 + w3*w3);
    return round( wght * bins);
  }
  else return 0;
}

int parallel_weights(int vert1, int vert2, vtkSmartPointer<vtkPolyData> surface){
  //Computes the weights for a maxcut of forward vs backward
  //Option
  int bins = 10000000;
  vector<int> common_cells =  vert_pair_cells_shared(vert1, vert2, surface);
  if( common_cells.size() == 2){
    double wght = dot_ngbs( common_cells[0], common_cells[1], surface, 0);
    return round( wght * bins);
  }
  else return 0;
}

int perpin_weight_by_curv(int vert1, int vert2, vtkSmartPointer<vtkPolyData> surface){
  //What if we weight proportionally with the difference in curvatures
  //If two cells seem flat, we shouldn't care if they are aligned or not
  int bins = 1000000000;
  vector<int> common_cells =  vert_pair_cells_shared(vert1, vert2, surface);
  if( common_cells.size() == 2){
    auto curvels = surface->GetCellData()->GetArray(curvatures_name);
    double curvs0[2];
    curvels->GetTuple( common_cells[0], curvs0);
    double sure0 = curvs0[0]-curvs0[1];
    double curvs1[2];
    curvels->GetTuple( common_cells[1], curvs1);
    double sure1 = curvs1[0]-curvs1[1];
    double w0 = dot_ngbs( common_cells[0], common_cells[1], surface, 0);
    double w1 = dot_ngbs( common_cells[0], common_cells[1], surface, 1);
    double w3 = dot_ngbs( common_cells[1], common_cells[0], surface, 1);
    double wght = w0*w0 - 0.5*(w1*w1 + w3*w3);
    if (sure0>sure1) wght *= sure0;
    else wght *= sure1;
    return round( wght * bins);
  }
  else return 0;
}



vector<bool> maxcut_surface_cells(vtkSmartPointer<vtkPolyData> surface, int (*weight_func)(int, int, vtkSmartPointer<vtkPolyData>) ){
  //Build a list of the edges of the simplicial surface to have integer identifies for the blowup
  //In the surface 1-skeleton  we direct edges (i,j) to have increasing vertex id so i<j
  //The reversed edge (j,i) has edge_id(j,i) = edge_id(i,j)+num_edges
  //So an easy test to see if two edges are the reversals of each other is if
  // e1 % num_edges == e2 % num_edges
  int num_surf_cells = surface->GetNumberOfCells();
  int num_surf_verts = surface->GetNumberOfPoints();
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
    else x = orig_edge_lookup[ {v2, v1} ] +num_surf_edges;
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
        int wght = weight_func(foo, bar, surface);
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
  //maybe we dont have to match with integers, look into that
  //The size of the rounding heavily impacts performance
  PerfectMatching blossom_matching( 2*num_surf_edges , blwp_edge_list.size() );
  for (auto an_edge: blwp_edge_list){
    blossom_matching.AddEdge(an_edge[0], an_edge[1], an_edge[2]);
  }
  blossom_matching.Solve();
  //blossom_matching.GetMatch()

  //Build the graph to find the cut corresponding to the blossomV matching
  //Record the edges in the cut-set for easier cut computation
  vector<int> cut_set;
  auto cutted_cellgraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  cutted_cellgraph->SetNumberOfVertices(num_surf_cells);
  for(int foo=0; foo< num_surf_cells; foo++){
    vector<int> adj_tris = adjacent_triangles(foo, surface);
    for(int bar : adj_tris){
      if( foo < bar ) {
        vector<int> ptpair = cell_pair_verts_shared( foo, bar, surface);
        int ed_id = find_edge_id( ptpair[0], ptpair[1]);
        int match_id = blossom_matching.GetMatch(ed_id);
        //the graph shouldn't have any of the edges chosen by the matching
        //add in the other edges
        if( ed_id % num_surf_edges != match_id % num_surf_edges ){
          cutted_cellgraph->AddEdge(foo,bar);
        }
        else cut_set.push_back(ed_id);
      }
    }
  }
  //The connected components of the cell-graph after removing blossom's cut-set
  auto connectedComponents = vtkSmartPointer<vtkBoostConnectedComponents>::New();
  connectedComponents->SetInputData( cutted_cellgraph);
  connectedComponents->Update();
  vtkGraph* outputGraph = connectedComponents->GetOutput();
  vtkIntArray* components = vtkIntArray::SafeDownCast(outputGraph->GetVertexData()->GetArray("component"));
  surface->GetCellData()->AddArray(components);
  //Collect num_surf_cells per component. The most common decides the flip scheme.
  double comprange[2];
  components->GetRange( comprange );
  const int num_components = round(comprange[1])+1;
  //We compute the cut from the cut-set by two-coloring the graph of components.
  auto component_graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
  component_graph->SetNumberOfVertices(num_components);
  for (int an_edge : cut_set){
    array<int,2> vin_an_edge = orig_edge_list[an_edge];
    vector<int> cin_an_edge = vert_pair_cells_shared( vin_an_edge[0], vin_an_edge[1], surface);
    int comp0 = components->GetValue(cin_an_edge[0]);
    int comp1 = components->GetValue(cin_an_edge[1]);
    if( component_graph->GetEdgeId(comp0,comp1) == -1 ) component_graph->AddEdge(comp0,comp1);
  }
  //We two color by taking distance mod 2 in the graph of components from any vertex
  //Compute distance by breadth first search
  auto BFS = vtkSmartPointer<vtkBoostBreadthFirstSearch>::New();
  BFS->SetOriginVertex(0);
  BFS->SetInputData(component_graph);
  BFS->Update();
  vtkIntArray* level = vtkIntArray::SafeDownCast(
    BFS->GetOutput()->GetVertexData()->GetArray("BFS"));
  //decide the sides of the cut
  vector<bool> in_or_out;
  in_or_out.resize(num_surf_cells);
  array<int,2> side_sizes = {0,0};
  for(int foo =0; foo< num_surf_cells; foo++){
    const int comp = components->GetValue(foo);
    int bfdistfrom0 = level->GetValue(comp);
    bool is_in = ( ( bfdistfrom0 % 2 ) == 0 );
    in_or_out[foo]=is_in;
    if( is_in ) side_sizes[0] = side_sizes[0] + 1;
    else side_sizes[1] = side_sizes[1] + 1;
  }
  if( side_sizes[0] < side_sizes[1] ){
    for(int foo=0; foo<num_surf_cells; foo++) in_or_out[foo] = !in_or_out[foo];
  }
  return in_or_out;
}



vtkSmartPointer<vtkPolyData> reorient_frame (vtkSmartPointer<vtkPolyData> surface ){
  //Compute a maxcut to decide what to flip perp
  vector<bool> good_line = maxcut_surface_cells(surface, perpin_weights);
  int const num_surf_cells = surface->GetNumberOfCells();
  auto flip90 = vtkSmartPointer<vtkIntArray>::New();
  flip90->SetNumberOfTuples(num_surf_cells);
  flip90->SetName("Flip Perpindicular");
  for(int foo=0; foo<num_surf_cells; foo++){
    if( good_line[foo] ) flip90->SetValue(foo,0);
    else flip90->SetValue(foo,1);
  }
  surface->GetCellData()->AddArray(flip90);
  //Rotate the vectors that were decided to flip perpindicular
  vtkSmartPointer<vtkDataArray> euc_field = surface->GetCellData()->GetArray(euc_minor_name);
  vtkSmartPointer<vtkDataArray> normals=surface->GetCellData()->GetNormals();
  for(int foo=0; foo<num_surf_cells; foo++){
    if( !good_line[foo] ){
      double tri[3];
      euc_field->GetTuple( foo, tri);
      vtkVector3d ve(tri);
      double nortemp[3];
      normals->GetTuple( foo, nortemp);
      vtkVector3d vn(nortemp);
      ve = vn.Cross(ve);
      euc_field->SetTuple3( foo, ve[0], ve[1], ve[2]);
    }
  }
  surface->GetCellData()->AddArray( euc_field );

  //Compute a maxcut to decide what to flip parallel
  vector<bool> good_forward = maxcut_surface_cells(surface, parallel_weights);
  auto flip180 = vtkSmartPointer<vtkIntArray>::New();
  flip180->SetNumberOfTuples(num_surf_cells);
  flip180->SetName("Flip Parallel");
  for(int foo=0; foo<num_surf_cells; foo++){
    if( good_forward[foo] ) flip180->SetValue(foo,0);
    else flip180->SetValue(foo,1);
  }
  surface->GetCellData()->AddArray(flip180);

  for(int foo=0; foo<num_surf_cells; foo++){
    if( !good_forward[foo] ){
      double tri[3];
      euc_field->GetTuple( foo, tri);
      euc_field->SetTuple3( foo, -tri[0], -tri[1], -tri[2]);
    }
  }
  surface->GetCellData()->AddArray( euc_field );
  return surface;
}


//Function to align the vector field consistently
vtkSmartPointer<vtkPolyData> align_field(vtkSmartPointer<vtkPolyData> surface ){
  //Remembers the order cells were decided in the alignment
  int const num_surf_cells = surface->GetNumberOfCells();
  auto cell_order = vtkSmartPointer<vtkIntArray>::New();
  cell_order->SetNumberOfComponents(1);
  cell_order->SetNumberOfTuples( num_surf_cells );
  cell_order->SetName("Cell Order");

  //Remembers the patch cells belonged to
  auto aligned_patch = vtkSmartPointer<vtkIntArray>::New();
  aligned_patch->SetNumberOfComponents(1);
  aligned_patch->SetNumberOfTuples( num_surf_cells );
  aligned_patch->SetName("Concurrence Patch");
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
          aligned_patch->SetTuple1( current.id, patch_count);
          cell_to_patch[current.id] = patch_count;
          cell_decided[current.id] = true;
          num_cells_in_patch++;
          last_cell_seen = current.id;
          //Look for neighbors
          vector<int> ngbs = adjacent_triangles(current.id , surface);
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
        aligned_patch->SetTuple1(last_cell_seen, 0);
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
    vector<int> ngbs = adjacent_triangles(foo, surface);
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

  //Decide how to rotate patches to align moving out from the biggest patch.
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
  surface->GetCellData()->AddArray( aligned_patch );
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
  if (surface->GetCellData()->HasArray(euc_minor_name)){
    surface->GetCellData()->RemoveArray(euc_minor_name);
  }
  if (surface->GetCellData()->HasArray(curvatures_name)){
    surface->GetCellData()->RemoveArray(curvatures_name);
  }
  surface->GetCellData()->AddArray(euc_field);
  surface->GetCellData()->AddArray(curvatures);

  surface = reorient_frame(surface);
  surface = align_field(surface);

  return surface;
}



vtkSmartPointer<vtkPolyData> compute_tubular_parametrization(vtkSmartPointer<vtkPolyData> surface){
  //Contstruct a function u: {mesh points}->RR with gradient ~= minor curvature field
  int num_pts = surface->GetNumberOfPoints();
  auto cell_pts = vtkSmartPointer<vtkIdList>::New();
  int num_surf_cells = surface->GetNumberOfCells();
  //Construct the discrete gradient matrix for the surface
  //The size is slightly larger to get constraint u(0)=0
  vector< Eigen::Triplet<double> > grad_entrylist;
  grad_entrylist.reserve(4 * num_surf_cells+1);
  vector< Eigen::Triplet<double> > curv_entrylist;
  curv_entrylist.reserve(2 * num_surf_cells+1);
  double curv_entry[2];

  vtkSmartPointer<vtkDataArray> euc_field = surface->GetCellData()->GetArray(euc_minor_name);
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
  //Porting data from vtkArray to a Eigen::Matrix
  //This should probably be done when you compute it in the first place
  for(int foo=0; foo< num_surf_cells; foo++){
    surface->GetCellPoints(foo,cell_pts);
    int pt0=cell_pts->GetId(0), pt1=cell_pts->GetId(1), pt2=cell_pts->GetId(2);
    //cout<<foo<<"cell's pts"<<pt0<<pt1<<pt2<<endl;
    minor_curv_field->GetTuple(foo,curv_entry);
    double veccomp0 = curv_entry[0];
    double veccomp1 = curv_entry[1];
    double vec_norm = sqrt(veccomp0 * veccomp0 + veccomp1 * veccomp1);
    //This threshold for normalization is totally arbitrary
    //and should be informed by stats on the curvature.
    double normalizing_thresh = 1e-10;
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
  grad_entrylist.push_back( Eigen::Triplet<double> (2*num_surf_cells, 0, 1) );
  curv_entrylist.push_back( Eigen::Triplet<double> (2*num_surf_cells, 0, 0) );

  //Build Matrix
  Eigen::SparseMatrix<double> surf_gradient( 2*num_surf_cells + 1, num_pts);
  surf_gradient.setFromTriplets(grad_entrylist.begin(), grad_entrylist.end());
  //Build target vector i.e. the curvature field
  Eigen::SparseMatrix<double> curv_field( 2*num_surf_cells + 1, 1);
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


//Apply a uniform neighbor blut to a scalar field
vtkSmartPointer<vtkPolyData> blur_scalar_field (vtkSmartPointer<vtkPolyData> surface, const char* oldname , const char* newname ){
  int num_pts = surface->GetNumberOfPoints();
  auto tube_param = surface->GetPointData()->GetArray(oldname);
  auto blurred_sc = vtkSmartPointer<vtkDoubleArray>::New();
  blurred_sc->SetNumberOfComponents(1);
  blurred_sc->SetNumberOfTuples(num_pts);
  blurred_sc->SetName(newname);
  for(int foo=0; foo<num_pts; foo++){
    double blur_val = 0.0;
    double tmp[1];
    tube_param->GetTuple(foo, tmp);
    blur_val+= tmp[0];
    set<int> ngbs = adjacent_vertices(foo, surface);
    int degp1 = ngbs.size()+1;
    for (int ngb : ngbs){
      double tmp2;
      tmp2 = tube_param->GetTuple1(ngb);
      blur_val+= tmp2;
    }
    blurred_sc->SetTuple1( foo, blur_val / degp1);
  }
  surface->GetPointData()->AddArray(blurred_sc);
  return surface;
}


array<double,2> local_coordinate_vec(int at_cell,  vtkSmartPointer<vtkPolyData> surface){
  vtkSmartPointer<vtkDataArray> vec_field = surface->GetCellData()->GetArray(euc_minor_name);
  //Compute the vector field in local coordinates.
  //Compute the local coordinate representation over
  // local edge basis e0=p2-p1 and e1=p0-p2
  //Get the points as vectors
  double temp_pts[3][3];
  vtkVector3d pts[3];
  auto cell_pts = vtkSmartPointer<vtkIdList>::New();
  surface->GetCellPoints( at_cell, cell_pts);
  for(int bar=0; bar<3; bar++){
    surface->GetPoint( cell_pts->GetId(bar), temp_pts[bar] );
    pts[bar] = vtkVector3d(temp_pts[bar]);
  }
  //points are got.
  //Compute the triangle edges
  vtkVector3d edges[3];
  for(int bar=0; bar<3; bar ++){
    //Edge vectors for orientation 012
    edges[bar] = pts [ (2+bar)%3 ] - pts [ (1+bar)%3 ];
  }
  //edges are computed.
  //
  // double e00 = edges[0].Dot(edges[0]), e01 = edges[0].Dot(edges[1]), e11 = edges[1].Dot(edges[1]);
  //double edet = e00*e11-e01*e01;
  double temp[3];
  vec_field->GetTuple(at_cell,temp);
  vtkVector3d vecvect(temp);
  double e0dw = edges[0].Dot(vecvect) , e1dw = edges[1].Dot(vecvect) ;
  array<double,2> local_vec={
    e0dw,
    e1dw
  };
  return local_vec;
}



vtkSmartPointer<vtkPolyData> potential_from_vec_field(vtkSmartPointer<vtkPolyData> surface){
  //Contstruct a function u: {mesh points}->RR with gradient ~= vec_field
  int num_pts = surface->GetNumberOfPoints();
  int num_surf_cells = surface->GetNumberOfCells();
  //Construct the discrete gradient matrix for the surface
  //The size is slightly larger to get constraint u(0)=0
  vector< Eigen::Triplet<double> > grad_entrylist;
  grad_entrylist.reserve(4 * num_surf_cells+1);
  vector< Eigen::Triplet<double> > vec_field_entrylist;
  vec_field_entrylist.reserve(2 * num_surf_cells+1);

  //Porting data from vtkArray to a Eigen::Matrix
  //This should probably be done when you compute it in the first place
  for(int foo=0; foo< num_surf_cells; foo++){
    auto cell_pts = vtkSmartPointer<vtkIdList>::New();
    surface->GetCellPoints(foo,cell_pts);
    int pt0=cell_pts->GetId(0), pt1=cell_pts->GetId(1), pt2=cell_pts->GetId(2);
    array<double,2> local_vec = local_coordinate_vec(foo, surface);
    vec_field_entrylist.push_back( Eigen::Triplet<double> ( foo, 0, local_vec[0] ) );
    vec_field_entrylist.push_back( Eigen::Triplet<double> ( foo + num_surf_cells, 0, local_vec[1] ) );
    //The graditent per cell has component0 = u(p2)-u(p1)
    //and component1 = u(p0)-u(p2)
    grad_entrylist.push_back( Eigen::Triplet<double> (foo, pt2, 1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (foo, pt1, -1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (foo+ num_surf_cells, pt0, 1) );
    grad_entrylist.push_back( Eigen::Triplet<double> (foo+ num_surf_cells, pt2, -1) );
  }
  //Add constraint u(0)=0
  grad_entrylist.push_back( Eigen::Triplet<double> (2*num_surf_cells, 0, 1) );
  vec_field_entrylist.push_back( Eigen::Triplet<double> (2*num_surf_cells, 0, 0) );

  //Build Matrix
  Eigen::SparseMatrix<double> surf_gradient( 2*num_surf_cells + 1, num_pts);
  surf_gradient.setFromTriplets(grad_entrylist.begin(), grad_entrylist.end());
  //Build target vector i.e. the curvature field
  Eigen::SparseMatrix<double> vec_field_mat( 2*num_surf_cells + 1, 1);
  vec_field_mat.setFromTriplets(vec_field_entrylist.begin(), vec_field_entrylist.end());

  Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
  cout<<"Factoring surface gradient."<<endl;
  solver.compute(surf_gradient);
  Eigen::SparseMatrix<double> u_param(num_pts,1);
  cout<<"Constructing potential for vector field."<<endl;
  time_t time1,time2;
  time1 = time(0);
  u_param = solver.solve(vec_field_mat);
  time2 = time(0);
  int it_took = difftime(time2,time1);
  cout<<"Integration finished in " << it_took/60 << " min and " << it_took%60 << " sec." <<endl;
  cout<<"Required " << solver.iterations() << " iterations obtaining error " << solver.error() <<endl;
  auto potential_field = vtkSmartPointer<vtkDoubleArray>::New();
  potential_field->SetNumberOfComponents(1);
  potential_field->SetNumberOfTuples(num_pts);
  potential_field->SetName(potential_name);
  for(int foo=0; foo<num_pts; foo++){
    potential_field->SetTuple1(foo, u_param.coeffRef(foo,0));
  }
  surface->GetPointData()->AddArray(potential_field);
  return surface;
}



vtkSmartPointer<vtkPolyData> build_small_simplicial_example()
{
  string example_name ("grid"); //"simple" "simple2" "cylinder" "grid"
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
    double cyl_wobble_amp = .3, cyl_wobble_freq=2 , cylrad=1, cylheight=350;
    //Suggest parameters .3,2,1
    int xtris = 50;
    int ytris = 50;
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
      if(!is_down){
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
    // double const y_shift=sqrt(3)/2;
    //double const x_shift=1;
    for(int foo =0; foo < (xtris)*(ytris+1); foo++){
      int pt_col = foo%xtris;
      int pt_row = foo/xtris;
      param_pts[foo][0]=pt_col + 0.5*pt_row;
      param_pts[foo][1]=pt_row;
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
      pt[2] = cylheight*p2/(ytris);
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
    int gridx=100;
    int gridy=100;
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




int main(int argc, const char* argv[]){
  const char* output_name;
  auto surface = vtkSmartPointer<vtkPolyData>::New();
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

  // auto reeber = vtkSmartPointer<vtkReebGraph>::New();
  // reeber->Build(surface, potential_name);
  // reeber->Simplify(.5 , nullptr);
  //
  vtkSmartPointer<vtkMutableUndirectedGraph> g =
    vtkSmartPointer<vtkMutableUndirectedGraph>::New();

  vtkIdType v1 = g->AddVertex();
  vtkIdType v2 = g->AddVertex();
  vtkIdType v3 = g->AddVertex();

  g->AddEdge(v1, v2);
  g->AddEdge(v3, v2);
  vtkSmartPointer<vtkGraphLayoutView> graphLayoutView =
    vtkSmartPointer<vtkGraphLayoutView>::New();
  graphLayoutView->AddRepresentationFromInput(g);
  graphLayoutView->SetLayoutStrategy("Simple 2D");
  graphLayoutView->ResetCamera();
  graphLayoutView->Render();

  // vtkSimple2DLayoutStrategy::SafeDownCast(graphLayoutView->GetLayoutStrategy())->SetRandomSeed(0);
  graphLayoutView->GetInteractor()->Start();

  // vtkSmartPointer<vtkGraphWriter> writer = vtkSmartPointer<vtkGraphWriter>::New();
  // writer->SetFileName(output_name);
  // writer->SetInputData(reeber);
  // writer->Write();
}
