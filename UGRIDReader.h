/*
 * UGRIDReader.h
 *
 *  Created on: Jun 14, 2013
 *     Authors: Felicia Brisc, Stefan Vater
 *
 */


#ifndef UGRIDREADER_H
#define UGRIDREADER_H

//#include "vtkIONetCDFModule.h" // For export macro
#include "vtkUnstructuredGridAlgorithm.h"

#include <vtk_netcdfcpp.h>
#include <vector>

class NcFile;
class NcVar;

class vtkDataArraySelection;
class vtkCallbackCommand;
class vtkUGRIDReaderInternal;

//class VTKIONETCDF_EXPORT UGRIDReader : public vtkUnstructuredGridAlgorithm
class UGRIDReader : public vtkUnstructuredGridAlgorithm
{
public:
  static UGRIDReader *New();
  vtkTypeMacro(UGRIDReader,vtkUnstructuredGridAlgorithm);
 // void PrintSelf(ostream& os, vtkIndent indent);


  // returns 1 if this file can be read and 0 if the file cannot be read.
  static int CanReadFile(const char* fileName);

  void SetFileName(const char* fileName);
  vtkGetStringMacro(FileName);

  // Description (GUI functions):
  // The following methods allow selective reading of solutions fields.
  // By default, ALL data fields on the nodes are read, but this can
  // be modified.
  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);
  void DisableAllPointArrays();
  void EnableAllPointArrays();

  int GetNumberOfCellArrays();
  const char* GetCellArrayName(int index);
  int GetCellArrayStatus(const char* name);
  void SetCellArrayStatus(const char* name, int status);
  void DisableAllCellArrays();
  void EnableAllCellArrays();

  void Enable(const char* name);
  void Disable(const char* name);

  // Get the reader's output
  vtkUnstructuredGrid *GetOutput();
  vtkUnstructuredGrid *GetOutput(int index);


protected:
  UGRIDReader();
  ~UGRIDReader();

  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*);

  virtual int RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *);

  virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                                  vtkInformationVector *);

  static void SelectionCallback(vtkObject* caller, unsigned long eid, void* clientdata, void* calldata);

  vtkDataArraySelection* NodeDataArraySelection;
  vtkDataArraySelection* CellDataArraySelection;

  vtkCallbackCommand* SelectionObserver;

//   int GetNcVars ();
  int BuildVarArrays();
  NcDim* IdentifyTimeDimension();


private:
  UGRIDReader(const UGRIDReader&);  // not implemented.
  void operator=(const UGRIDReader&);  // not implemented.

  // The file name of the file that contains all of the data
  // (mesh and variable fields).
  char* FileName;
  char* CurrentFileName;
  vtkSetStringMacro(CurrentFileName);

  double* TimeSteps;
  long NumberOfTimeSteps;

  int NumberOfCellVars;
  int NumberOfNodeVars;

  // The NetCDF file descriptors.  NULL indicates they haven't been opened.
  NcFile* PointsFile;

  bool InfoRequested;

};


#endif /* UGRIDREADER_H */
