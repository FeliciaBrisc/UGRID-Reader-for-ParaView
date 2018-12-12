/*
 * UGRIDReader.cxx
 *
 *  Created on: Jun 14, 2013
 *     Authors: Felicia Brisc, Stefan Vater
 *
 */


#include "UGRIDReader.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkShortArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolygon.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkStdString.h"

#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <vtk_netcdfcpp.h>


class NcVar;
class NcFile;

using namespace std;

vtkStandardNewMacro(UGRIDReader);

#define MAX_VARS 100


// TODO: Paraview should not crash when opening file with "invalid" structure

//------------------------------------------------------------------------------
UGRIDReader::UGRIDReader()
{
  // turn on debugging messages with the following command, otherwise comment it
  //this->DebugOn();

  this->FileName        = NULL;
  this->CurrentFileName = NULL;
  this->PointsFile      = NULL;

  this->TimeSteps = NULL;
  this->NumberOfTimeSteps = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);

  InfoRequested = true;

  // setup selection callback to modify this object when array selection changes
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->NodeDataArraySelection = vtkDataArraySelection::New();
  this->SelectionObserver = vtkCallbackCommand::New();
  this->SelectionObserver->SetCallback(&UGRIDReader::SelectionCallback);
  this->SelectionObserver->SetClientData(this);
  this->CellDataArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
  this->NodeDataArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
}


//------------------------------------------------------------------------------
UGRIDReader::~UGRIDReader()
{
  this->SetFileName(NULL);
  this->SetCurrentFileName(NULL);

  if(this->PointsFile)
    {
    delete this->PointsFile;
    this->PointsFile = NULL;
    }

  if(this->TimeSteps)
    {
    delete []this->TimeSteps;
    this->TimeSteps = NULL;
    }

  if (this->CellDataArraySelection)
    {
    this->CellDataArraySelection->Delete();
    this->CellDataArraySelection = NULL;
    }

  if (this->NodeDataArraySelection)
    {
    this->NodeDataArraySelection->Delete();
    this->NodeDataArraySelection = NULL;
    }

  if (this->SelectionObserver)
    {
    this->SelectionObserver->Delete();
    this->SelectionObserver = NULL;
    }

}


// for the cases where the time dimension is named other than "time"
NcDim* UGRIDReader::IdentifyTimeDimension()
{
  NcDim* timeDimension = NULL;
  int numberOfDims = this->PointsFile->num_vars();
  for (int n=0; n<numberOfDims; n++){
	  NcVar* currDimension = this->PointsFile->get_var(n);

	  for (int a=0; a<currDimension->num_atts(); a++)		{
		  NcAtt* crtAttribute = currDimension->get_att(a);

		  if (strcmp(crtAttribute->name(), "standard_name") == 0) {
			  NcValues* crtAttValues = crtAttribute->values();

			  if (strcmp(crtAttValues->as_string(0), "time") == 0) {
				  //this is out time dimension
	     		  vtkDebugMacro(<< "The time dimension is: " << currDimension->name());
				  timeDimension = this->PointsFile->get_dim(currDimension->name());
			  }
		  }
	}
  }
  return timeDimension;
}

//------------------------------------------------------------------------------
int UGRIDReader::CanReadFile(const char* fileName)
{
  NcFile file(fileName, NcFile::ReadOnly);
  return file.is_valid();
}


//------------------------------------------------------------------------------
void UGRIDReader::SetFileName(const char* fileName)
{
  vtkDebugMacro(<<" setting FileName to " << (fileName?fileName:"(null)") );
  if ( this->FileName == NULL && fileName == NULL)
    {
    return;
    }
  if ( this->FileName && fileName && (!strcmp(this->FileName,fileName)))
    {
    return;
    }
  if(this->PointsFile)
    {
    delete this->PointsFile;
    this->PointsFile = NULL;
    }
  if (this->FileName)
    {
    delete [] this->FileName;
    this->FileName = NULL;
    }
  if (fileName)
    {
    size_t n = strlen(fileName) + 1;
    char *cp1 =  new char[n];
    const char *cp2 = (fileName);
    this->FileName = cp1;
    do
      {
      *cp1++ = *cp2++;
      }
    while ( --n );
    }
   else
    {
    this->FileName = NULL;
    }
  this->Modified();
}


//------------------------------------------------------------------------------
int UGRIDReader::RequestInformation(
  vtkInformation* vtkNotUsed(reqInfo),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  if(this->FileName == NULL) {
    vtkWarningMacro("Missing a file name.");
    return 0;
  }

  if(this->CurrentFileName != NULL && strcmp(this->CurrentFileName, this->FileName) != 0) {
    delete this->PointsFile;
    this->PointsFile = NULL;
    this->SetCurrentFileName(NULL);
  }
  if(this->PointsFile == NULL) {
    this->PointsFile = new NcFile(this->FileName, NcFile::ReadOnly);
    if(this->PointsFile->is_valid() == 0) {
      vtkErrorMacro(<< "Can't read file " << this->FileName);
      delete this->PointsFile;
      this->PointsFile = NULL;
      return 0;
    }
    this->SetCurrentFileName(this->FileName);
  }

  // GUI
  if (InfoRequested) {
    if (!BuildVarArrays()) {
      return 0;
    }
    InfoRequested = false;
  }

  // get the time dimension
  // TODO: be able to process data without time dimension

  //for now only identify the time dimension
  NcDim* timeDimension = this->IdentifyTimeDimension();


  //NcDim* timeDimension = this->PointsFile->get_dim("time");
  if(timeDimension == NULL) {
//     vtkErrorMacro("Cannot find the number of time steps (time dimension).");
//     return 0;
    vtkDebugMacro(<<"Cannot find the number of time steps (time dimension).");
    this->NumberOfTimeSteps = 0;
  } else {
    this->NumberOfTimeSteps = timeDimension->size();
  }

  // Get ParaView information pointer
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  if (this->NumberOfTimeSteps > 0) {
    if(this->TimeSteps) {
      delete []this->TimeSteps;
    }
    this->TimeSteps = new double[this->NumberOfTimeSteps];

    NcVar* timeVar = this->PointsFile->get_var(timeDimension->name());
    timeVar->get(this->TimeSteps, this->NumberOfTimeSteps);

    // tell the pipeline what steps are available
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),this->TimeSteps, this->NumberOfTimeSteps);

    // range is required to get GUI to show things  --- these values are displayed in the "Properties-Information tab"
    double tRange[2] = {this->TimeSteps[0], this->TimeSteps[this->NumberOfTimeSteps - 1]};

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), tRange, 2);
  } else {
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
  }


  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

  return 1;
}


//------------------------------------------------------------------------------
int UGRIDReader::RequestUpdateExtent(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  if(this->FileName == NULL)
    {
    vtkWarningMacro("Missing a file name.");
    return 0;
    }
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  return 1;
}


//------------------------------------------------------------------------------
int UGRIDReader::RequestData(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  if(this->FileName == NULL )
    {
    vtkWarningMacro("Missing a file name.");
    return 0;
    }

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *output =
    vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro(<<"Reading NetCDF UGRID file.");
  this->SetProgress(0);

  // set the NetCDF error handler to not kill the application.
  // upon exiting this method the error handler will be restored
  // to its previous state.
  NcError ncError(NcError::verbose_nonfatal);

  // parse file for variables and their associated mesh topologies
  NcValues* vals;
  std::vector<NcVar*> meshes;
  std::vector<NcVar*> nodevars, facevars; //Felicia 2018

  // get all face variables
  for(int i=0; i < this->NumberOfCellVars; i++) {
    NcVar* variable = this->PointsFile->get_var(this->GetCellArrayName(i));

    NcAtt* apm = variable->get_att("mesh");
    if (apm == 0)
      {
      vtkErrorMacro("mesh attribute not found!");
      return 0;
      }

    vals = apm->values();
    // check if mesh is already loaded
    int j = 0;
    while (j < meshes.size() && strcmp(vals->as_string(0), meshes[j]->name()) != 0)
      j++;
    if (j == meshes.size())
      {
      // mesh is not yet loaded
      // load mesh variable and check if it is valid
      NcVar* meshvar = this->PointsFile->get_var(vals->as_string(0));
      // TODO: check if mesh variable exists
      NcAtt* ap = meshvar->get_att("cf_role");
      if (ap != 0)
        {
        vals = ap->values();
        if (strcmp(vals->as_string(0), "mesh_topology") == 0)
          meshes.push_back(meshvar);
        else
          {
          vtkErrorMacro("This is not a valid mesh variable!");
          return 0;
          }
        }
      else
        {
        vtkErrorMacro("This is not a valid mesh variable!");
        return 0;
        }
      }

    // TODO: we probably have to associate each variable with the corresponding mesh?!
    // verify that the data is specified per cell or node
    NcAtt* apv = variable->get_att("location");
    NcValues* vals = apv->values();

    if (strcmp(vals->as_string(0), "face") == 0)
      facevars.push_back(variable);
    else
      {
      vtkErrorMacro("data is not specified per face");
      return 0;
      }
  }

  // get all node variables
  for(int i=0; i < this->NumberOfNodeVars; i++) {
    NcVar* variable = this->PointsFile->get_var(this->GetPointArrayName(i));

    NcAtt* apm = variable->get_att("mesh");
    if (apm == 0)
      {
      vtkErrorMacro("mesh attribute not found!");
      return 0;
      }

    vals = apm->values();
    // check if mesh is already loaded
    int j = 0;
    while (j < meshes.size() && strcmp(vals->as_string(0), meshes[j]->name()) != 0)
      j++;
    if (j == meshes.size())
      {
      // mesh is not yet loaded
      // load mesh variable and check if it is valid
      NcVar* meshvar = this->PointsFile->get_var(vals->as_string(0));
      // TODO: check if mesh variable exists
      // check if mesh is valid, i.e. has attribute "cf_role" with value "mesh_topology"
      NcAtt* ap = meshvar->get_att("cf_role");
      if (ap != 0)
        {
        vals = ap->values();
        if (strcmp(vals->as_string(0), "mesh_topology") == 0)
          meshes.push_back(meshvar);
        else
          {
          vtkErrorMacro("This is not a valid mesh variable!");
          return 0;
          }
        }
      else
        {
        vtkErrorMacro("This is not a valid mesh variable!");
        return 0;
        }
      }

    // TODO: we probably have to associate each variable with the corresponding mesh?!
    // verify that the data is specified per cell or node
    NcAtt* apv = variable->get_att("location");
    NcValues* vals = apv->values();

    if (strcmp(vals->as_string(0), "node") == 0)
      nodevars.push_back(variable);
    else
      {
      vtkErrorMacro("data is not specified per node");
      return 0;
      }
  }


  // some debug output
  string dbgstrm;
  for(size_t i=0; i!=meshes.size(); ++i)
    dbgstrm = dbgstrm+" "+meshes[i]->name();
  vtkDebugMacro(<< "meshes found: " << dbgstrm);

  string dbgstrfv;
  for(size_t i=0; i!=facevars.size(); ++i)
    dbgstrfv = dbgstrfv+" "+facevars[i]->name();
  vtkDebugMacro(<< "face variables loaded: " << dbgstrfv);

  string dbgstrnv;
  for(size_t i=0; i!=nodevars.size(); ++i)
    dbgstrnv = dbgstrnv+" "+nodevars[i]->name();
  vtkDebugMacro(<< "node variables loaded: " << dbgstrnv);



  // get mesh variable, TODO: we have to iterate over meshes!
  // maybe we have to use vtkMultiBlockDataSet,
  // see VTK/IO/LSDyna/vtkLSDynaReader.cxx for an example
  NcVar* meshvar = meshes[0];
  NcAtt* ap;

  // get information about coordinate variable names
  ap   = meshvar->get_att("node_coordinates");
  vals = ap->values();
  stringstream ss(vals->as_string(0));
  string item;
  std::vector<string> coords;
  while (getline(ss, item, ' '))
    coords.push_back(item);

  string dbgstrc;
  for(size_t i=0; i!=coords.size(); ++i)
    dbgstrc = dbgstrc+" "+coords[i];

  vtkDebugMacro(<< "Coordinate variables: " << dbgstrc);

  NcVar* node_x = this->PointsFile->get_var((char*)coords[0].c_str());
  NcVar* node_y = this->PointsFile->get_var((char*)coords[1].c_str());

  if (node_x == NULL || node_y == NULL)
    {
    vtkErrorMacro("Cannot get coordinate variables.");
    return 0;
    }

  // get information about connectivity
  ap   = meshvar->get_att("face_node_connectivity");
  vals = ap->values();
  NcVar* connectivity = this->PointsFile->get_var(vals->as_string(0));

  if(connectivity == NULL) {
    vtkErrorMacro("Cannot find mesh (face-node) connectivity.");
    return 0;
  }

  // get dimension for nodes, face-node-connectivity in NetCDF structure
  NcDim* node_dim            = node_x->get_dim(0);
  NcDim* face_node_conn_dim0 = connectivity->get_dim(0);
  NcDim* face_node_conn_dim1 = connectivity->get_dim(1);

  int nodes_num = node_dim->size();
  int faces_num, max_nodes_per_face, stride_face, stride_npfc;

  // get information about number of faces and max nodes per face
  ap = meshvar->get_att("face_dimension");
  if (ap == 0) {
    vtkWarningMacro("Face dimension not defined by mesh variable. Assuming face dimension is the first in connectivity array.");
    faces_num          = face_node_conn_dim0->size();
    max_nodes_per_face = face_node_conn_dim1->size();
    stride_face        = max_nodes_per_face;
    stride_npfc        = 1;
  } else {
    vals = ap->values();
    if (strcmp(vals->as_string(0), face_node_conn_dim0->name()) == 0) {
      faces_num          = face_node_conn_dim0->size();
      max_nodes_per_face = face_node_conn_dim1->size();
      stride_face        = max_nodes_per_face;
      stride_npfc        = 1;
    } else {
      faces_num          = face_node_conn_dim1->size();
      max_nodes_per_face = face_node_conn_dim0->size();
      stride_face        = 1;
      stride_npfc        = faces_num;
    }
  }

  // if max_nodes_per_face will be >3, the mesh contains quadrilaterals
  NcAtt* connFillValue;
  if (max_nodes_per_face > 3) {
    // we have a mixed mesh

    // we need to get the _FillValue of the connectivity array:
    connFillValue = connectivity->get_att("_FillValue");
    if (!connFillValue) {
      vtkErrorMacro("_FillValue attribute missing - The connectivity variable has to specify a _FillValue attribute");
      return 0;
    }
    vtkDebugMacro(<< "_FillValue is: " << connFillValue->as_int(0));
  }

  // get node coordinates
  vtkSmartPointer<vtkPoints> node_points = vtkSmartPointer<vtkPoints>::New();

  if(node_x->type() == ncDouble)
    {
    node_points->SetDataTypeToDouble();
    node_points->SetNumberOfPoints(nodes_num);
    std::vector<double> array_x(nodes_num), array_y(nodes_num);

    if(!node_x->get(&array_x[0], nodes_num))
      return 0;

    if(!node_y->get(&array_y[0], nodes_num))
      return 0;

    for(long i=0;i<nodes_num;i++)
      node_points->SetPoint(i, array_x[i], array_y[i], 1);

    }
  else
    {
    node_points->SetDataTypeToFloat();
    node_points->SetNumberOfPoints(nodes_num);

    std::vector<float> array_x(nodes_num), array_y(nodes_num);

    if(!node_x->get(&array_x[0], nodes_num))
      return 0;

    if(!node_y->get(&array_y[0], nodes_num))
      return 0;

    for(long i=0;i<nodes_num ;i++)
      node_points->SetPoint(i, array_x[i], array_y[i], 1);

    }

  // educated guess for progress
  this->SetProgress(.25);

  // get start index
  ap   = connectivity->get_att("start_index");
  vals = ap->values();
  int start_index = vals->as_int(0);
  // TODO: check that start_index is either 0 or 1

  vtkDebugMacro(<< "start index of face-node connectivity: " << start_index);

  // get connectivity array
  std::vector<int> faceConnectivity(max_nodes_per_face*faces_num);
  bool b = connectivity->get(&(faceConnectivity[0]), face_node_conn_dim0->size(),
                             face_node_conn_dim1->size());



  node_points->Modified();
  output->SetPoints(node_points); //this might be problematic with the _FillValue from the connectivity array in the case of mixed meshes

  // educated guess for progress
  this->SetProgress(.5);

  // collect the time step requested
  vtkInformationDoubleKey* timeKey =  static_cast<vtkInformationDoubleKey*>
                                      (vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  double dTime = 0.0;
  if (outInfo->Has(timeKey))
    {
    dTime = outInfo->Get(timeKey);
    }

  // actual time for the time step
  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), dTime);

  // index of the time step to request
  int timeStep = 0;
  while (timeStep < this->NumberOfTimeSteps && this->TimeSteps[timeStep] < dTime)
    {
    timeStep++;
    }

  // get all variables associated with faces (cells)
  for(size_t i=0; i!=facevars.size(); ++i)
    {
    NcVar* variable = facevars[i];

    //------GUI--------
    if ( !CellDataArraySelection->ArrayIsEnabled(variable->name()) )
      continue;
    //------------------

    if(variable->type() == ncDouble)
      {
      vtkDoubleArray* doubleArray = vtkDoubleArray::New();
      doubleArray->SetName(variable->name());
      doubleArray->SetNumberOfComponents(1);

      doubleArray->SetNumberOfTuples(faces_num);
      double* dataBlock = doubleArray->GetPointer(0);
      variable->get(dataBlock, 1, faces_num, 1);
      output->GetCellData()->AddArray(doubleArray);
      doubleArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(doubleArray->GetPointer(0), 1, faces_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncFloat)
      {
      vtkFloatArray* floatArray = vtkFloatArray::New();
      floatArray->SetName(variable->name());
      floatArray->SetNumberOfComponents(1);

      floatArray->SetNumberOfTuples(faces_num);
      float* dataBlock = floatArray->GetPointer(0);
      variable->get(dataBlock, 1, faces_num, 1);
      output->GetCellData()->AddArray(floatArray);
      floatArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(floatArray->GetPointer(0), 1, faces_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncInt)
      {
      vtkIntArray* intArray = vtkIntArray::New();
      intArray->SetName(variable->name());
      intArray->SetNumberOfComponents(1);

      intArray->SetNumberOfTuples(faces_num);
      int* dataBlock = intArray->GetPointer(0);
      variable->get(dataBlock, 1, faces_num, 1);
      output->GetCellData()->AddArray(intArray);
      intArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(intArray->GetPointer(0), 1, faces_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncShort)
      {
      vtkShortArray* shortArray = vtkShortArray::New();
      shortArray->SetName(variable->name());
      shortArray->SetNumberOfComponents(1);

      shortArray->SetNumberOfTuples(faces_num);
      short* dataBlock = shortArray->GetPointer(0);
      variable->get(dataBlock, 1, faces_num, 1);
      output->GetCellData()->AddArray(shortArray);
      shortArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(shortArray->GetPointer(0), 1, faces_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }

    }

  // get all variables associated with nodes
  for(size_t i=0; i!=nodevars.size(); ++i)
    {
    NcVar* variable = nodevars[i];

    //------GUI--------
    if ( !NodeDataArraySelection->ArrayIsEnabled(variable->name()) )
      continue;
    //------------------

    if(variable->type() == ncDouble)
      {
      vtkDoubleArray* doubleArray = vtkDoubleArray::New();
      doubleArray->SetName(variable->name());
      doubleArray->SetNumberOfComponents(1);

      doubleArray->SetNumberOfTuples(nodes_num);
      double* dataBlock = doubleArray->GetPointer(0);
      variable->get(dataBlock, 1, nodes_num, 1);
      output->GetPointData()->AddArray(doubleArray);
      doubleArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(doubleArray->GetPointer(0), 1, nodes_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncFloat)
      {
      vtkFloatArray* floatArray = vtkFloatArray::New();
      floatArray->SetName(variable->name());
      floatArray->SetNumberOfComponents(1);

      floatArray->SetNumberOfTuples(nodes_num);
      float* dataBlock = floatArray->GetPointer(0);
      variable->get(dataBlock, 1, nodes_num, 1);
      output->GetPointData()->AddArray(floatArray);
      floatArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(floatArray->GetPointer(0), 1, nodes_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncInt)
      {
      vtkIntArray* intArray = vtkIntArray::New();
      intArray->SetName(variable->name());
      intArray->SetNumberOfComponents(1);

      intArray->SetNumberOfTuples(nodes_num);
      int* dataBlock = intArray->GetPointer(0);
      variable->get(dataBlock, 1, nodes_num, 1);
      output->GetPointData()->AddArray(intArray);
      intArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(intArray->GetPointer(0), 1, nodes_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }
    else if (variable->type() == ncShort)
      {
      vtkShortArray* shortArray = vtkShortArray::New();
      shortArray->SetName(variable->name());
      shortArray->SetNumberOfComponents(1);

      shortArray->SetNumberOfTuples(nodes_num);
      short* dataBlock = shortArray->GetPointer(0);
      variable->get(dataBlock, 1, nodes_num, 1);
      output->GetPointData()->AddArray(shortArray);
      shortArray->Delete();

      variable->set_cur(timeStep, 0);

      if(!variable->get(shortArray->GetPointer(0), 1, nodes_num))
        {
        vtkErrorMacro("Problem getting NetCDF variable " << variable->name());
        return 0;
        }
      }

    }

  this->SetProgress(.75);  // guess for progress

  // create the actual cells --- this is the allocate nedeed before the InsertNextCell operation
  output->Allocate(faces_num);

  bool isQuad     = false;
  bool isTriangle = false;

  for(long face_i=0; face_i<faces_num; face_i++) {
    vtkSmartPointer<vtkIdTypeArray> pointIds = vtkSmartPointer<vtkIdTypeArray>::New();

    for(int j=0; j<max_nodes_per_face; j++) {
      vtkIdType nextId = faceConnectivity[j*stride_npfc + face_i*stride_face];
      if (max_nodes_per_face>3 && nextId == connFillValue->as_int(0)) {
        // we found an "empty" node with a _FillValue in this mixed mesh,
        // therefore this cell is a triangle
        isQuad     = false;
        isTriangle = true;
        // skip this empty node:
        continue;
      } else if (max_nodes_per_face > 3) {
        pointIds->InsertNextValue(nextId-start_index);
        isQuad = true;
      } else { //this is a simple triangle mesh
        pointIds->InsertNextValue(nextId-start_index);
        isTriangle=true;
      }

    }

    int ids_size = pointIds->GetNumberOfTuples();
    vtkIdType *pointIdsArray = new vtkIdType[ids_size];
    for (int c=0; c<ids_size; c++)
      pointIdsArray[c] = pointIds->GetValue(c);

    if (isQuad && !isTriangle) {
      output->InsertNextCell(VTK_QUAD, 4, pointIdsArray );
    } else {
      output->InsertNextCell(VTK_TRIANGLE, 3, pointIdsArray);
    }
    isQuad     = false;
    isTriangle = false;
    pointIds = NULL;
    delete [] pointIdsArray;
  }

  vtkDebugMacro(<< "Read " << output->GetNumberOfPoints() << " points,"
                << output->GetNumberOfCells() << " cells.\n");

  return 1;
}


//------------------------------------------------------------------------------
//                    GUI  functions
//------------------------------------------------------------------------------

// build the selection arrays for points and cells in the GUI
int UGRIDReader::BuildVarArrays() {
  // set the NetCDF error handler to not kill the application.
  // upon exiting this method the error handler will be restored
  // to its previous state.
  NcError ncError(NcError::verbose_nonfatal);

  // parse file for variables and their associated mesh topologies
  for(int i=0; i < this->PointsFile->num_vars(); i++) {
    NcVar* variable = this->PointsFile->get_var(i);

    NcAtt* apm = variable->get_att("mesh");
    if (apm != 0) {
      // TODO: we probably have to associate each variable with the corresponding mesh?!
      // verify that the data is specified per cell or node
      NcAtt* apv = variable->get_att("location");
      NcValues* vals = apv->values();

      if (strcmp(vals->as_string(0), "face") == 0) {
        this->CellDataArraySelection->AddArray(variable->name());
      } else if ( strcmp(vals->as_string(0), "node") == 0) {
        this->NodeDataArraySelection->AddArray(variable->name());
      } else {
        vtkErrorMacro("Cannot identify if data is specified per face or per node");
        return 0;
      }
    } else {
      NcAtt* aunits = variable->get_att("units");
      if (aunits != 0) {
        NcValues* vunits = aunits->values();
        vtkStdString units = vunits->as_string(0);
        vtkDebugMacro(<< "Variable " << variable->name() << " has unit " << units);
      }
    }
  }

  this->NumberOfCellVars = this->GetNumberOfCellArrays();
  this->NumberOfNodeVars = this->GetNumberOfPointArrays();

  string dbgstrfv;
  for(size_t i=0; i!=this->NumberOfCellVars; ++i)
    dbgstrfv = dbgstrfv+" "+this->GetCellArrayName(i);
  vtkDebugMacro(<< this->NumberOfCellVars << " face variables found: " << dbgstrfv);

  string dbgstrnv;
  for(size_t i=0; i!=this->NumberOfNodeVars; ++i)
    dbgstrnv = dbgstrnv+" "+this->GetPointArrayName(i);
  vtkDebugMacro(<< this->NumberOfNodeVars << " node variables found: " << dbgstrnv);


  vtkDebugMacro( << "NumberOfCellVars: " << this->NumberOfCellVars << ", "
                 << "NumberOfNodeVars: " << this->NumberOfNodeVars << endl);

  // figure out what variables to visualize -
  for (int var = 0; var < this->NumberOfNodeVars; var++) {
    this->NodeDataArraySelection->EnableArray(this->NodeDataArraySelection->GetArrayName(var));
  }

  for (int var = 0; var < this->NumberOfCellVars; var++) {
    this->CellDataArraySelection->EnableArray(this->CellDataArraySelection->GetArrayName(var));
  }


  return(1);
}


//------------------------------------------------------------------------------
//callback if the user selects a variable
void UGRIDReader::SelectionCallback(vtkObject*, unsigned long vtkNotUsed(eventid),
                                    void* clientdata, void* vtkNotUsed(calldata))
{
  static_cast<UGRIDReader*>(clientdata)->Modified();
}


//------------------------------------------------------------------------------
// return the output
vtkUnstructuredGrid* UGRIDReader::GetOutput()
{
  return this->GetOutput(0);
}


//------------------------------------------------------------------------------
//returns the output given an id
vtkUnstructuredGrid* UGRIDReader::GetOutput(int idx)
{
  if (idx)
    {
    return NULL;
    }
  else
    {
    return vtkUnstructuredGrid::SafeDownCast( this->GetOutputDataObject(idx) );
    }
}


//------------------------------------------------------------------------------
void UGRIDReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}


//------------------------------------------------------------------------------
void UGRIDReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}


//------------------------------------------------------------------------------
int UGRIDReader::GetNumberOfPointArrays()
{
  return this->NodeDataArraySelection->GetNumberOfArrays();
}


//------------------------------------------------------------------------------
int UGRIDReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}

//------------------------------------------------------------------------------
void UGRIDReader::EnableAllPointArrays()
{
  this->NodeDataArraySelection->EnableAllArrays();
}


//------------------------------------------------------------------------------
void UGRIDReader::DisableAllPointArrays()
{
  this->NodeDataArraySelection->DisableAllArrays();
}


//------------------------------------------------------------------------------
void UGRIDReader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}


//------------------------------------------------------------------------------
void UGRIDReader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}


//------------------------------------------------------------------------------
const char* UGRIDReader::GetPointArrayName(int index)
{
  return this->NodeDataArraySelection->GetArrayName(index);
}


//------------------------------------------------------------------------------
int UGRIDReader::GetPointArrayStatus(const char* name)
{
  return this->NodeDataArraySelection->ArrayIsEnabled(name);
}


//------------------------------------------------------------------------------
void UGRIDReader::SetPointArrayStatus(const char* name, int status)
{
  if (status)
    {
    this->NodeDataArraySelection->EnableArray(name);
    }
  else
    {
    this->NodeDataArraySelection->DisableArray(name);
    }
}


//------------------------------------------------------------------------------
const char* UGRIDReader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}


//------------------------------------------------------------------------------
int UGRIDReader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}


//------------------------------------------------------------------------------
void UGRIDReader::SetCellArrayStatus(const char* name, int status)
{
  if (status)
    {
    this->CellDataArraySelection->EnableArray(name);
    }
  else
    {
    this->CellDataArraySelection->DisableArray(name);
    }
}


//------------------------------------------------------------------------------
/*
void UGRIDReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(NULL)") << endl;

  if(this->PointsFile)
    {
    os << indent << "PointsFile: " << this->PointsFile << endl;
    }
  else
    {
    os << indent << "PointsFile: (NULL)" << endl;
    }

}
*/
