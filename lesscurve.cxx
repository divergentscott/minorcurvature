#include <array>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSTLReader.h>
#include "vtkSmartPointer.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include <vtkXMLPolyDataWriter.h>
#include <vtkVector.h>
#include <vtkOBJReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>

vtkVector3d operator+(vtkVector3d u, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=u[foo]+v[foo];
  }
  return w;
};
vtkVector3d operator-(vtkVector3d u, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=u[foo]-v[foo];
  }
  return w;
};
vtkVector3d operator*(double scalar, vtkVector3d v){
  vtkVector3d w;
  for(int foo=0;foo<3;foo++){
    w[foo]=scalar*v[foo];
  }
  return w;
};
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


vtkSmartPointer<vtkPolyData> build_small_simplicial_example()
{
  std::string example_name ("cylinder"); //"simple" "simple2" "cylinder"
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
      //std::cout << pt[0] << pt[1] << pt[2] << std::endl;
    }
    int cellList[4][3] = {
      {0,1,3},
      {1,3,4},
      {1,2,4},
      {3,4,5}
    };
    for(int foo = 0; foo<sizeof(cellList)/sizeof(cellList[0]); foo++){
      triangles->InsertNextCell(3);
      for(int bar = 0; bar<sizeof(cellList[0])/sizeof(cellList[0][0]); bar++){
        triangles->InsertCellPoint(cellList[foo][bar]);
      }
    }
  }

  if (example_name.compare("simple2") == 0)
  {
    double angle0 = 0.0, angle1 = 0.0 , angle2 = M_PI/4;
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
      // std::cout << points->GetNumberOfPoints() << " " << pt_pos[foo][0] << pt_pos[foo][1] << pt_pos[foo][2] << std::endl;
    }
    int cellList[4][3] = {
      {0,1,2},
      {0,3,1},
      {0,2,4},
      {1,5,2}
    };
    for(int foo = 0; foo<sizeof(cellList)/sizeof(cellList[0]); foo++){
      triangles->InsertNextCell(sizeof(cellList[0])/sizeof(cellList[0][0]));
      for(int bar = 0; bar<sizeof(cellList[0])/sizeof(cellList[0][0]); bar++){
        triangles->InsertCellPoint(cellList[foo][bar]);
      }
    }
  }
  if (example_name.compare("cylinder") == 0)
  {
    int xtris = 70;
    int ytris = 70;
    double param_pts[xtris*(ytris+1)][2];
    int const num_tris = 2 * xtris * ytris;
    int cell_list[num_tris][3];
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
    double const x_shift=1;
    for(int foo =0; foo < (xtris)*(ytris+1); foo++){
      int pt_col = foo%xtris;
      int pt_row = foo/xtris;
      param_pts[foo][0]=pt_col + 0.5*pt_row;
      param_pts[foo][1]=y_shift*pt_row;
      //std::cout<<pt_col<<pt_row<<" "<<param_pts[foo][0]<<" "<<param_pts[foo][1]<<std::endl;
    }
    //Apply parametrization to obtain points in R3
    int num_pts = sizeof(param_pts)/sizeof(param_pts[0]);
    points->SetNumberOfPoints(num_pts);
    for (int foo = 0; foo<num_pts; foo++){
      double pt[3];
      double p1 = param_pts[foo][0];
      double p2 = param_pts[foo][1];
      double cyl_wobble_amp = .1, cyl_wobble_freq=3 , cylrad=1;
      double rad_now = cyl_wobble_amp * cos(2*M_PI*cyl_wobble_freq*p2/ytris ) + cylrad;
      pt[0] = rad_now * cos(2*M_PI*p1/xtris);
      pt[1] = rad_now * sin(2*M_PI*p1/xtris);
      pt[2] = 2*M_PI*p2/ytris;
      points->SetPoint(foo, pt[0], pt[1], pt[2]);
      //std::cout << pt[0] << pt[1] << pt[2] << std::endl;
    }
    for(int foo = 0; foo<sizeof(cell_list)/sizeof(cell_list[0]); foo++){
      triangles->InsertNextCell(3);
      //std::cout<<std::endl<<"Build cell" <<foo;
      for(int bar = 0; bar<sizeof(cell_list[0])/sizeof(cell_list[0][0]); bar++){
        triangles->InsertCellPoint(cell_list[foo][bar]);
        //std::cout<<"point "<< cell_list[foo][bar];
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

vtkSmartPointer<vtkPolyData> compute_minor_curvature_field(vtkSmartPointer<vtkPolyData> surface){
  //std::cout << "Computing minor curvature field." << std::endl;
  vtkDataArray *normals=surface->GetCellData()->GetNormals();
  //Initialize memory for the minor curvautre field.
  //The field is defined per triangle of the surface.
  //The field is the direction of the shaper operator eignevalue with the lower
  //magintude curvature.
  //The field is expressed over in the barycentric coordinate basis
  // side2=p1-p0 and side1=p2-p0 for triangle /_p0p1p2
  const char* minor_name="MinorCurvatureDirectionField";
  auto minor_curv_field = vtkSmartPointer<vtkFloatArray>::New();
  minor_curv_field->SetNumberOfComponents(2);
  minor_curv_field->SetNumberOfTuples(surface->GetNumberOfCells());
  minor_curv_field->SetName(minor_name);

  //The minor curvature field here is expressed over the standar basis in RR3
  const char* euc_name="3SpaceMinorCurvatureDirectionField";
  auto euc_minor_curv_field = vtkSmartPointer<vtkFloatArray>::New();
  euc_minor_curv_field->SetNumberOfComponents(3);
  euc_minor_curv_field->SetNumberOfTuples(surface->GetNumberOfCells());
  euc_minor_curv_field->SetName(euc_name);

  //Loop through the cells to compute the shape operator locally and extract
  //the minor curvature direction
 //std::cout<<"Inside curve field, I see cells numbering: "<< surface->GetNumberOfCells() << std::endl;

  for(vtkIdType foo=0; foo< surface->GetNumberOfCells(); foo++){
    ////std::cout<<"Working Cell "<< foo << std::endl;
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
    int edge_ends[3][3] = {{2,1},{0,2},{1,0}};
    double angles[3];
    //Compute the triangle edges
    for(int bar=0; bar<3; bar ++){
      //Edge vectors for orientation 012
      edges[bar] = pts [ edge_ends[bar][0] ] - pts [ edge_ends[bar][1] ];
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
        cell_pts->GetId(edge_ends[bar][0]),
        cell_pts->GetId(edge_ends[bar][1]),
        cell_ngb);
      //If there is no neighboring cell, i.e. the edge is surface boundary,
      //then we assume the surface is simply flat in that direction.
      //Otherwise add the angle scaled projection in that direction.
      if(cell_ngb->GetNumberOfIds()>0){
        normals->GetTuple(cell_ngb->GetId(0),temp_n);
        vtkVector3d ngb_normal(temp_n);
        //The sign of the angle is taken to agree with right-hand rotation
        //about the edge vector.
        angles[bar] = asin( edges[bar].Dot( normal.Cross(ngb_normal) ) / edge_lengths[bar] ) ;
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
  double shape_small_eig;
  if(shape_trace>0){
    shape_small_eig=(shape_trace + sqrt(shape_discr))/2.0;
  }
  else{
    shape_small_eig=(shape_trace - sqrt(shape_discr))/2.0;
  }
  double ortho_curv_vect[2];
  if (fabs(shape_small_eig-shape00) > fabs(shape_small_eig-shape11)){
    ortho_curv_vect[0] = shape01;
    ortho_curv_vect[1] = shape_small_eig-shape00;
  }
  else{
    ortho_curv_vect[0] = shape_small_eig-shape11;
    ortho_curv_vect[1] = shape01;
  }
  //What's the right normalization for the field?
  if(ortho_curv_vect[0]*ortho_curv_vect[1]){
    double magni = sqrt(ortho_curv_vect[0]*ortho_curv_vect[0]+ortho_curv_vect[1]*ortho_curv_vect[1]);
    //magni *= .001*( edges[0].Cross(edges[1]) ).Norm();
    ortho_curv_vect[0]*=shape_small_eig/magni;
    ortho_curv_vect[1]*=shape_small_eig/magni;
  }
  //Store the vector in the standard R3 basis
  vtkVector3d euc_curv_vect =  ortho_curv_vect[0] * orthobase0 + ortho_curv_vect[1] * orthobase1 ;
    //Hacky attempt at local continuity
  if(euc_curv_vect[2] < 0){
    euc_curv_vect = -1 * euc_curv_vect;
  }
  euc_minor_curv_field->SetTuple3( foo , euc_curv_vect[0], euc_curv_vect[1], euc_curv_vect[2]);
  //Also store the vector in the local edge basis e0=p2-p1 and p1=p0-p2
  double e00 = edges[0].Dot(edges[0]), e01 = edges[0].Dot(edges[1]), e11 = edges[1].Dot(edges[1]);
  double edet = e00*e11-e01*e01;
  double e0dw = edges[0].Dot(euc_curv_vect)  , e1dw = edges[1].Dot(euc_curv_vect) ;
  double minor_curv_vect[2] ={
    (e11*e0dw - e01 * e1dw) /edet,
    (e00*e1dw - e01 * e0dw) /edet,
  };
  minor_curv_field->SetTuple( foo , minor_curv_vect);
  }

  //Looks noisey. Try a uniform window blur on cell neighbors.
  //There should be a fast way to do this with a convolution filter....

  const char* blur_name="BlurredField";
  auto blurred_field = vtkSmartPointer<vtkFloatArray>::New();
  blurred_field->SetNumberOfComponents(3);
  blurred_field->SetNumberOfTuples(surface->GetNumberOfCells());
  blurred_field->SetName(blur_name);
  bool apply_blur = false;
  if(apply_blur){
    for(int foo = 0; foo<surface->GetNumberOfCells(); foo++){
      double tempf[3];
      euc_minor_curv_field->GetTuple( foo , tempf);
      vtkVector3d field_here(tempf);
      int edge_ends[3][3] = {{2,1},{0,2},{1,0}};
      auto cell_pts = vtkSmartPointer<vtkIdList>::New();
      surface->GetCellPoints(foo,cell_pts);
      auto cell_ngb = vtkSmartPointer<vtkIdList>::New();
      for(int bar=0;bar<3;bar++){
        surface->GetCellEdgeNeighbors(foo,
        cell_pts->GetId(edge_ends[bar][0]),
        cell_pts->GetId(edge_ends[bar][0]),
        cell_ngb);
        if(cell_ngb->GetNumberOfIds()>0){
          euc_minor_curv_field->GetTuple(cell_ngb->GetId(0), tempf);
          vtkVector3d field_there(tempf);
          field_here=field_here+field_there;
        }
      }
      field_here.Normalize();
      double * w = &field_here[0];
      blurred_field->SetTuple(foo, w);
    }
  }

  // If an old curvature array already exists we remove it
  if (surface->GetCellData()->HasArray(minor_name))
    surface->GetCellData()->RemoveArray(minor_name);

  if (surface->GetCellData()->HasArray(euc_name))
    surface->GetCellData()->RemoveArray(euc_name);

  surface->GetCellData()->AddArray(minor_curv_field);
  surface->GetCellData()->AddArray(euc_minor_curv_field);

  if(apply_blur){
    if (surface->GetCellData()->HasArray(blur_name))
      surface->GetCellData()->RemoveArray(blur_name);
    surface->GetCellData()->AddArray(blurred_field);
  }

  return surface;
}

int main(int argc, const char* argv[])
{
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
    int g;
    surface = build_small_simplicial_example();
    // std::cout<<"After build cells numbering "<<surface->GetNumberOfCells() << std::endl;
  }

  // auto surface = vtkSmartPointer<vtkPolyData>::New();
  // vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  // reader->SetFileName(argv[0]);
  // reader->Update();
  // surface = reader->GetOutput();
  //  surface->BuildCells();
  auto normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

  //Compute normal vectors for all cells
  normalGenerator->SetInputData(surface);
  normalGenerator->ComputePointNormalsOff();
  normalGenerator->ComputeCellNormalsOn();
  normalGenerator->SetSplitting(0);
  normalGenerator->AutoOrientNormalsOn();
  normalGenerator->Update();
  surface = normalGenerator->GetOutput();
  //surface->DeleteCells();
  surface->BuildLinks();
  surface = compute_minor_curvature_field(surface);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(output_name);
  writer->SetInputData(surface);
  writer->Write();
}
