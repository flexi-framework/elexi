/*
!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
*/

#include <visuReader.h>
#include <../../plugin_visu.h>

#include <hdf5.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkLine.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkVertex.h>
#if USE_MPI
#include <vtkPDistributedDataFilter.h>
#endif /* USE_MPI */

#include <libgen.h>
#include <unistd.h>
#include <algorithm>
#include <sstream>

// MPI
#include <vtkMultiProcessController.h>
vtkStandardNewMacro(visuReader);
vtkCxxSetObjectMacro(visuReader, Controller, vtkMultiProcessController);

#define SWRITE(x) {if (ProcessId == 0) std::cout << "@@@ " << this << " " << x << "\n";};

/*
 * Construtor of State Reader
 */
visuReader::visuReader()
{
   SWRITE("visuReader");
   this->FileName = NULL;
   this->NVisu = 0;
   this->NCalc = 0;
   this->NodeTypeVisu = NULL;
   this->Avg2d = 0;
   this->DGonly = 0;
   this->ParameterFileOverwrite = NULL;
   this->MeshFileOverwrite = NULL;
   this->SetNumberOfInputPorts(0);
#if USE_PARTICLES
   this->SetNumberOfOutputPorts(4);
#else
	 this->SetNumberOfOutputPorts(2);
#endif

   // Setup the selection callback to modify this object when an array
   // selection is changed.
   // Used to tell the visuReader, that we (un)selected a state,primite or derived quantity
   // and that the 'Apply' button becomes clickable to reload the data (load the selected quantities)

   this->SelectionObserver     = vtkCallbackCommand::New();
   this->SelectionObserver->SetCallback(&visuReader::SelectionModifiedCallback);
   this->SelectionObserver->SetClientData(this);
   // create array for state,primitive and derived quantities
   this->VarDataArraySelection = vtkDataArraySelection::New();
   this->BCDataArraySelection  = vtkDataArraySelection::New();
#if USE_PARTICLES
   this->VarParticleDataArraySelection = vtkDataArraySelection::New();
#endif
   // add an observer
   this->VarDataArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
   this->BCDataArraySelection ->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
#if USE_PARTICLES
   this->VarParticleDataArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
#endif
}

/*
 * This function is called when a file is inserted into the Pipeline-browser.
 * Here we load the variable names (state,primitive, derived).
 */
int visuReader::RequestInformation(vtkInformation *,
      vtkInformationVector **,
      vtkInformationVector *outputVector)
{
   // We take the first state file and use it to read the varnames
   SWRITE("RequestInformation");

   // Set up MPI communicator
   this->Controller = NULL;
   this->SetController(vtkMultiProcessController::GetGlobalController());
   if (this->Controller == NULL) {
      NumProcesses = 1;
      ProcessId    = 0;
   } else {
      NumProcesses = this->Controller->GetNumberOfProcesses();
      ProcessId    = this->Controller->GetLocalProcessId();
   }

   vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
   mpiComm = MPI_COMM_NULL;
   if (communicator) {
      mpiComm = *(communicator->GetMPIComm()->GetHandle());
   }
   // get the info object of the output ports
   vtkSmartPointer<vtkInformation> outInfoVolume  = outputVector->GetInformationObject(0);
   vtkSmartPointer<vtkInformation> outInfoSurface = outputVector->GetInformationObject(1);
#if USE_PARTICLES
   vtkSmartPointer<vtkInformation> outInfoPart    = outputVector->GetInformationObject(2);
   vtkSmartPointer<vtkInformation> outInfoImpact = outputVector->GetInformationObject(3);
#endif /*USE_PARTICLES*/

   // sets the number of pieces to the number of processsors
   outInfoVolume ->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
   outInfoSurface->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
#if USE_PARTICLES
   outInfoPart   ->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
   outInfoImpact ->Set(CAN_HANDLE_PIECE_REQUEST(), 1);
#endif /*USE_PARTICLES*/

   // RequestInformation may be called before AddFileName, thus the arrays with timesteps and
   // file names may be empty.  In this case, simply leave the function. It will be called again
   // and with the files loaded. This does not make any sense...
   if (Timesteps.empty()) {
      std::cout << "No Filenames given, skipping...\n";
      return 1;
   }

   // if we have more then one file loaded at once (timeseries)
   // we have to set the number and range of the timesteps
   double timeRange[2] = {Timesteps[0], Timesteps[Timesteps.size()-1]};
   outInfoVolume ->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &Timesteps[0], Timesteps.size());
   outInfoSurface->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &Timesteps[0], Timesteps.size());
#if USE_PARTICLES
   outInfoPart   ->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &Timesteps[0], Timesteps.size());
   outInfoImpact ->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &Timesteps[0], Timesteps.size());
#endif /*USE_PARTICLES*/
   outInfoVolume ->Set (vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
   outInfoSurface->Set (vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

   // Change to directory of state file (path of mesh file is stored relative to path of state file)
   char* dir = strdup(FileNames[0].c_str());
   dir = dirname(dir);
   int ret = chdir(dir);
   if (ret != 0) {
      SWRITE("Directory of state file not found: " << dir);
      return 0;
   }

   // convert the MPI communicator to a Fortran communicator
   int fcomm;
   fcomm = MPI_Comm_c2f(mpiComm);
   MPI_Barrier(mpiComm);

   struct CharARRAY varnames;
   struct CharARRAY bcnames;
#if USE_PARTICLES
   struct CharARRAY partnames;
#endif /*USE_PARTICLES*/
   // Call Posti-function requestInformation:
   // This function returns the varnames of state, primitive and derived quantities
   int strlen_state = strlen(FileNames[0].c_str());
   int strlen_mesh;
   if (MeshFileOverwrite == NULL ) {
      strlen_mesh = 0;
   } else {
      strlen_mesh = strlen(MeshFileOverwrite);
   }
   __mod_visu_cwrapper_MOD_visu_requestinformation(&fcomm, &strlen_state, FileNames[0].c_str(), &strlen_mesh, MeshFileOverwrite, &varnames, &bcnames
#if USE_PARTICLES
                                                  ,&partnames
#endif /*USE_PARTICLES*/
                                                  );

   MPI_Barrier(mpiComm);

   // We copy the varnames to the corresponding DataArraySelection objects.
   // These objects are used to build the gui.
   // (see the functions below:
   //   DisableAllVarArrays, EnableAllVarArrays, GetNumberOfVarArrays, GetVarArrayName, GetVarArrayStatus, SetVarArrayStatus,
   // )
   for (int iVar=0; iVar<varnames.len/255; iVar++) {
      char tmps[255];
      strncpy(tmps, varnames.data+iVar*255, 255);
      std::string varname(tmps);
      varname = varname.substr(0,varname.find(" "));

      if (!this->VarDataArraySelection->ArrayExists(varname.c_str())) {
         // function DisableArray deselects the varname in the gui and if not existend inserts the varname
         this->VarDataArraySelection->DisableArray(varname.c_str());


         // Select Density, FV_Elems by default
         if (varname.compare("Density") == 0) this->VarDataArraySelection->EnableArray(varname.c_str());
#if FV_ENABLED == 1
         if (varname.compare("ElemData:FV_Elems") == 0) this->VarDataArraySelection->EnableArray(varname.c_str());
#elif FV_ENABLED == 2
         if (varname.compare("ElemData:FV_alpha") == 0) this->VarDataArraySelection->EnableArray(varname.c_str());
#endif
      }
   }
   for (int iVar=0; iVar<bcnames.len/255; iVar++) {
      char tmps[255];
      strncpy(tmps, bcnames.data+iVar*255, 255);
      std::string bcname(tmps);
      bcname = bcname.substr(0,bcname.find(" "));

      if (!this->BCDataArraySelection->ArrayExists(bcname.c_str())) {
         // function DisableArray deselects the bcname in the gui and if not existend inserts the bcname
         this->BCDataArraySelection->DisableArray(bcname.c_str());
      }
   }
#if USE_PARTICLES
	 for (int iVar=0; iVar<partnames.len/255; iVar++) {
      char tmps[255];
      strncpy(tmps, partnames.data+iVar*255, 255);
      std::string partname(tmps);
      partname = partname.substr(0,partname.find(" "));

      if (!this->VarParticleDataArraySelection->ArrayExists(partname.c_str())) {
         // function DisableArray deselects the partname in the gui and if not existend inserts the partname
         this->VarParticleDataArraySelection->DisableArray(partname.c_str());
      }
   }
#endif
   return 1;
}

/*
 * This function is called whenever a filename is loaded.
 * If multiple files are selected in the file-dialog, this function is called multiple times.
 * Attention: For multiple files, we assume a timeseries.
 */
void visuReader::AddFileName(const char* filename_in) {
   SWRITE("AddFileName "<<filename_in);
   // append the filename to the list of filenames
   this->FileNames.push_back(filename_in);
   this->Modified();

   // open the file with HDF5 and read the attribute 'time' to build a timeseries
   hid_t state = H5Fopen(filename_in, H5F_ACC_RDONLY, H5P_DEFAULT);
   double time;

   // check if attribute time exits
   htri_t exists = H5Aexists(state, "Time");

   // only access the attribute if it exists
   if (exists > 0){
      hid_t attr = H5Aopen(state, "Time", H5P_DEFAULT);
      SWRITE("attribute Time "<<attr);
      // only write back valid values
      if (attr > -1){
         hid_t attr_type = H5Aget_type( attr );
         H5Aread(attr, attr_type, &time);
         Timesteps.push_back(time);
      }
      H5Aclose(attr);
   } else {
      Timesteps.push_back(0.);
   }
   H5Fclose(state);
}

void visuReader::RemoveAllFileNames() {
   this->FileNames.clear();
   this->Timesteps.clear();
}

// Get number of state file in the filenames list closest to the requested time
int visuReader::FindClosestTimeStep(double requestedTimeValue)
{
   int ts = 0;
   double mindist = fabs(Timesteps[0] - requestedTimeValue);

   for (unsigned int i = 0; i < Timesteps.size(); i++) {
      double tsv = Timesteps[i];
      double dist = fabs(tsv - requestedTimeValue);
      if (dist < mindist) {
         mindist = dist;
         ts = i;
      }
   }
   return ts;
}

vtkStringArray* visuReader::GetNodeTypeVisuList() {
   vtkStringArray* arr = vtkStringArray::New();
   arr->InsertNextValue("VISU");
   arr->InsertNextValue("GAUSS");
   arr->InsertNextValue("GAUSS-LOBATTO");
   arr->InsertNextValue("VISU_INNER");
   return arr;
}

/*
 * This function is called, when the user presses the 'Apply' button.
 * Here we call the Posti and load all the data.
 */

int visuReader::RequestData(
      vtkInformation *vtkNotUsed(request),
      vtkInformationVector **vtkNotUsed(inputVector),
      vtkInformationVector *outputVector)
{
   RequestInformation(NULL,NULL,outputVector);
   SWRITE("RequestData");

   // get the index of the timestep to load
   int timestepToLoad = 0;
   std::string FileToLoad;
   vtkSmartPointer<vtkInformation> outInfoVolume  = outputVector->GetInformationObject(0);
   vtkSmartPointer<vtkInformation> outInfoSurface = outputVector->GetInformationObject(1);
#if USE_PARTICLES
   vtkSmartPointer<vtkInformation> outInfoPart    = outputVector->GetInformationObject(2);
   vtkSmartPointer<vtkInformation> outInfoImpact  = outputVector->GetInformationObject(3);
#endif
   if (outInfoVolume->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
      // get the requested time
      double requestedTimeValue = outInfoVolume->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      timestepToLoad = FindClosestTimeStep(requestedTimeValue);
   }
   if (timestepToLoad==0 && outInfoSurface->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
      // get the requested time
      double requestedTimeValue = outInfoSurface->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      timestepToLoad = FindClosestTimeStep(requestedTimeValue);
   }
#if USE_PARTICLES
   if (timestepToLoad==0 && outInfoPart->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
      // get the requested time
      double requestedTimeValue = outInfoPart->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      timestepToLoad = FindClosestTimeStep(requestedTimeValue);
   }
   if (timestepToLoad==0 && outInfoImpact->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
      // get the requested time
      double requestedTimeValue = outInfoImpact->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      timestepToLoad = FindClosestTimeStep(requestedTimeValue);
   }
#endif
   FileToLoad = FileNames[timestepToLoad];
   SWRITE("File to load "<<FileToLoad);


   // convert the MPI communicator to a fortran communicator
   int fcomm = MPI_Comm_c2f(mpiComm);
   MPI_Barrier(mpiComm); // all processes should call the Fortran code at the same time

   // get all variables selected for visualization
   // and check if they are the same as before (bug workaround, explanation see below)
   int nVars = VarDataArraySelection->GetNumberOfArrays();
   VarNames_selected.resize(nVars);
   for (int i = 0; i< nVars; ++i)
   {
      const char* name = VarDataArraySelection->GetArrayName(i);
      if ((name != NULL) && (name[0] != '\0')) {
        VarNames_selected[i] = VarDataArraySelection->ArrayIsEnabled(name);
      } else {
        VarNames_selected[i] = 0;
      }
   }

   int nBCs = BCDataArraySelection->GetNumberOfArrays();
   BCNames_selected.resize(nBCs);
   for (int i = 0; i< nBCs; ++i)
   {
      const char* name = BCDataArraySelection->GetArrayName(i);
      if ((name != NULL) && (name[0] != '\0')) {
        BCNames_selected[i] = BCDataArraySelection->ArrayIsEnabled(name);
      } else {
        BCNames_selected[i] = 0;
      }
   }

#if USE_PARTICLES
	 int nParts = VarParticleDataArraySelection->GetNumberOfArrays();
   PartNames_selected.resize(nParts);
   for (int i = 0; i< nParts; ++i)
   {
      const char* name = VarParticleDataArraySelection->GetArrayName(i);
      if ((name != NULL) && (name[0] != '\0')) {
        PartNames_selected[i] = VarParticleDataArraySelection->ArrayIsEnabled(name);
      } else {
        PartNames_selected[i] = 0;
      }
   }
#endif

   // Change to directory of state file (path of mesh file is stored relative to path of state file)
   char* dir = strdup(FileToLoad.c_str());
   dir = dirname(dir);
   int ret = chdir(dir);
   if (ret != 0) {
      SWRITE("Directory of state file not found: " << dir);
      return 0;
   }

   // Write temporary parameter file for Posti tool
   // TODO: only root-proc
   //if (ProcessId == 0) {

   // get temporary file for the posti parameter files
   setlocale(LC_ALL, "C");
   char posti_filename[] = "/tmp/f2p_posti_XXXXXX.ini";
   int posti_unit = mkstemps(posti_filename,4);

   // write settings to Posti parameter file
   dprintf(posti_unit, "NVisu = %d\n", NVisu); // insert NVisu
   dprintf(posti_unit, "NCalc = %d\n", NCalc); // insert NCalc
   dprintf(posti_unit, "NodeTypeVisu = %s\n", NodeTypeVisu); // insert NodeType
   dprintf(posti_unit, "Avg2D = %s\n", (this->Avg2d ? "T" : "F"));
   dprintf(posti_unit, "DGonly = %s\n", (this->DGonly ? "T" : "F"));
   if (MeshFileOverwrite != NULL ) {
     if (strlen(MeshFileOverwrite) > 0) {
        dprintf(posti_unit, "MeshFile = %s\n", MeshFileOverwrite);
     }
   }

   // write selected state varnames to the parameter file
   // if no varnames are selected, write flag to suppress standard vars to be visualized
   int noVisuVars = true;
   for (int i = 0; i< nVars; ++i)
   {
      if (VarNames_selected[i]) {
         noVisuVars = false;
         const char* name = VarDataArraySelection->GetArrayName(i);
         dprintf(posti_unit, "VarName = %s\n", name) ;
      }
   }
   if (noVisuVars) {
      dprintf(posti_unit, "noVisuVars = T\n") ;
   }
   for (int i = 0; i< nBCs; ++i)
   {
      if (BCNames_selected[i]) {
         const char* name = BCDataArraySelection->GetArrayName(i);
         dprintf(posti_unit, "BoundaryName = %s\n", name) ;
      }
   }
#if USE_PARTICLES
	 for (int i = 0; i< nParts; ++i)
   {
      if (PartNames_selected[i]) {
         const char* name = VarParticleDataArraySelection->GetArrayName(i);
         dprintf(posti_unit, "VarName = %s\n", name) ;
      }
   }
#endif
   close(posti_unit);

   MPI_Barrier(mpiComm); // all processes should call the Fortran code at the same time

   // call Posti tool (Fortran code)
   // the arrays coords_*, values_* and nodeids_* are allocated in the Posti tool
   // and contain the vtk data
   int strlen_prm;
   if (ParameterFileOverwrite == NULL) {
     strlen_prm = 0;
   } else {
     strlen_prm = strlen(ParameterFileOverwrite);
   }
   int strlen_posti = strlen(posti_filename);
   int strlen_state = strlen(FileToLoad.c_str());
   __mod_visu_cwrapper_MOD_visu_cwrapper(&fcomm,
         &this->HighOrder,&this->UseCurveds,
         &strlen_prm,   ParameterFileOverwrite,
         &strlen_posti, posti_filename,
         &strlen_state, FileToLoad.c_str(),
         &coords_DG    ,&values_DG    ,&nodeids_DG,
         &coords_FV    ,&values_FV    ,&nodeids_FV,
         &varnames,
         &coordsSurf_DG,&valuesSurf_DG,&nodeidsSurf_DG,
         &coordsSurf_FV,&valuesSurf_FV,&nodeidsSurf_FV,
         &varnamesSurf
#if USE_PARTICLES
        ,&coords_Part,&values_Part,&nodeids_Part,&varnames_Part,&components_Part,
         &coords_Impact,&values_Impact,&nodeids_Impact,&varnames_Impact,&components_Impact
#endif /*USE_PARTICLES*/
   );

   MPI_Barrier(mpiComm); // wait until all processors returned from the Fortran Posti code

   // get the MultiBlockDataset
   vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(outInfoVolume->Get(vtkDataObject::DATA_OBJECT()));
   if (!mb) {
      std::cout << "DownCast to MultiBlockDataset Failed!" << std::endl;
      return 0;
   }

   // adjust the number of Blocks in the MultiBlockDataset to 4 (DG and FV)
   SWRITE("Number of Blocks in MultiBlockDataset : " << mb->GetNumberOfBlocks());
   if (mb->GetNumberOfBlocks() < 2) {
      SWRITE("Create new DG and FV output blocks");
      mb->SetBlock(0, vtkUnstructuredGrid::New());
      mb->SetBlock(1, vtkUnstructuredGrid::New());
      //mb->SetBlock(2, vtkUnstructuredGrid::New());
      //mb->SetBlock(3, vtkUnstructuredGrid::New());
   }

   // Insert DG data into output
   InsertData(mb, 0, &coords_DG, &values_DG, &nodeids_DG, &varnames);

   // Insert FV data into output
   InsertData(mb, 1, &coords_FV, &values_FV, &nodeids_FV, &varnames);

   // get the MultiBlockDataset
   mb = vtkMultiBlockDataSet::SafeDownCast(outInfoSurface->Get(vtkDataObject::DATA_OBJECT()));
   if (!mb) {
      std::cout << "DownCast to MultiBlockDataset Failed!" << std::endl;
      return 0;
   }

   // adjust the number of Blocks in the MultiBlockDataset to 4 (DG and FV)
   SWRITE("Number of blocks in MultiBlockDataset : " << mb->GetNumberOfBlocks());
   if (mb->GetNumberOfBlocks() < 2) {
      SWRITE("Create new DG and FV output Blocks");
      mb->SetBlock(0, vtkUnstructuredGrid::New());
      mb->SetBlock(1, vtkUnstructuredGrid::New());
      //mb->SetBlock(2, vtkUnstructuredGrid::New());
      //mb->SetBlock(3, vtkUnstructuredGrid::New());
   }

   // Insert Surface DG data into output
   InsertData(mb, 0, &coordsSurf_DG, &valuesSurf_DG, &nodeidsSurf_DG, &varnamesSurf);

   // Insert Surface FV data into output
   InsertData(mb, 1, &coordsSurf_FV, &valuesSurf_FV, &nodeidsSurf_FV, &varnamesSurf);

#if USE_PARTICLES
	 // write PartData
   /* vtkPolyData* mb_part = vtkPolyData::SafeDownCast(outInfoPart->Get(vtkDataObject::DATA_OBJECT())); */
   vtkMultiBlockDataSet* mb_part = vtkMultiBlockDataSet::SafeDownCast(outInfoPart->Get(vtkDataObject::DATA_OBJECT()));
   if (!mb_part) {
      std::cout << "DownCast to MultiBlockDataset Failed!" << std::endl;
      return 0;
   }

   SWRITE("Number of Blocks in MultiBlockDataset : " << mb_part->GetNumberOfBlocks())
   if (mb_part->GetNumberOfBlocks() < 2) {
     SWRITE("Create new part output Block");
     /* mb_part->SetBlock(0, vtkUnstructuredGrid::New()); */
     mb_part->SetBlock(0, vtkPolyData::New());
   }

   // Insert particle data into output
   InsertPartData(mb_part,0, &coords_Part, &values_Part, &nodeids_Part, &varnames_Part, &components_Part);

   /* vtkMultiBlockDataSet* mb_impact = vtkMultiBlockDataSet::SafeDownCast(outInfoImpact->Get(vtkDataObject::DATA_OBJECT())); */
   vtkMultiBlockDataSet* mb_impact = vtkMultiBlockDataSet::SafeDownCast(outInfoImpact->Get(vtkDataObject::DATA_OBJECT()));
   if (!mb_impact) {
      std::cout << "DownCast to MultiBlockDataset Failed!" << std::endl;
      return 0;
   }

   SWRITE("Number of Blocks in MultiBlockDataset : " << mb_impact->GetNumberOfBlocks())
   if (mb_impact->GetNumberOfBlocks() < 2) {
     SWRITE("Create new impact output Block");
     /* mb_impact->SetBlock(0, vtkUnstructuredGrid::New()); */
     mb_impact->SetBlock(0, vtkPolyData::New());
   }

   // Insert impact data into output
   InsertPartData(mb_impact,0, &coords_Impact, &values_Impact, &nodeids_Impact, &varnames_Impact, &components_Impact);
#endif /*USE_PARTICLES*/

   __mod_visu_cwrapper_MOD_visu_dealloc_nodeids();

   // tell paraview to render data
   this -> Modified();

   MPI_Barrier(mpiComm); // synchronize again (needed?)
   SWRITE("RequestData finished");
   return 1;
}

/*
 * This function inserts the data, loaded by the Posti tool, into a ouput
 */
void visuReader::InsertData(vtkMultiBlockDataSet* mb    ,int blockno
                           ,struct DoubleARRAY* coords  ,struct DoubleARRAY* values
                           ,struct IntARRAY*    nodeids ,struct CharARRAY* varnames) {
    vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(mb->GetBlock(blockno));

    // create a 3D double array (must be 3D even we use a 2D Posti tool, since paraview holds the data in 3D)
    vtkSmartPointer <vtkDoubleArray> pdata = vtkSmartPointer<vtkDoubleArray>::New();
    pdata->SetNumberOfComponents(3); // 3D
    pdata->SetNumberOfTuples(coords->len/3);
    // copy coordinates
    double* ptr = pdata->GetPointer(0);
    for (long i = 0; i < coords->len; ++i)
    {
      *ptr++ = coords->data[i];
    }

    // create points array to be used for the output
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetData(pdata);
    output->SetPoints(points);

    // create cellarray
    vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();

    int CellLength;
    int CellType;
    if (this->HighOrder && blockno == 0){ // blockno: 0 = DG, 1 = FV
      if (coords->dim == 1) {
        // Use the nodeids to build lines
        CellLength = NVisu+1;
        CellType   = VTK_LAGRANGE_CURVE;
      } else if (coords->dim == 2) {
        // Use the nodeids to build quads
        CellLength = (NVisu+1)*(NVisu+1);
        CellType   = VTK_LAGRANGE_QUADRILATERAL;
      } else if (coords->dim == 3) {
        // Use the nodeids to build hexas
        CellLength = (NVisu+1)*(NVisu+1)*(NVisu+1);
        CellType   = VTK_LAGRANGE_HEXAHEDRON;
      } else {
        exit(1);
      }
    } else{
      if (coords->dim == 1) {
        // Use the nodeids to build lines
        CellLength = 2;
        CellType   = VTK_LINE;
      } else if (coords->dim == 2) {
        // Use the nodeids to build quads
        CellLength = 4;
        CellType   = VTK_QUAD;
      } else if (coords->dim == 3) {
        // Use the nodeids to build hexas
        CellLength = 8;
        CellType   = VTK_HEXAHEDRON;
      } else {
        exit(1);
      }
    }

    // (here we must copy the nodeids, we can not just assign the array of nodeids to some vtk-structure)
    int gi = 0;
    // loop over all cells
    for (int iCell=0; iCell<nodeids->len/CellLength; iCell++) {
      cellarray->InsertNextCell(CellLength);
      for (int i=0; i<CellLength; i++) {
        cellarray->InsertCellPoint(nodeids->data[gi]);
        gi++;
      }
    }

    // Use the nodeids to build all cells at once
    output->SetCells(CellType,cellarray);

    // assign the actual data, loaded by the Posti tool, to the output
    unsigned int nVar = varnames->len/255;

    if (nVar > 0) {
      unsigned int sizePerVar = values->len/nVar;
      int          dataPos    = 0;
      // loop over all loaded variables
      for (unsigned int iVar = 0; iVar < nVar; iVar++) {
        // For each variable, create a new array and set the number of components to 1
        // Each variable is loaded separately.
        // One might implement vector variables (velocity), but then must set number of componenets to 2/3
        vtkSmartPointer <vtkDoubleArray> vdata = vtkSmartPointer<vtkDoubleArray>::New();
        vdata->SetNumberOfComponents(1);
        vdata->SetNumberOfTuples(sizePerVar);
        // copy values
        double* ptr = vdata->GetPointer(0);
        for (long i = 0; i < sizePerVar; ++i)
        {
          *ptr++ = values->data[dataPos+i];
        }
        dataPos += sizePerVar;
        // set name of variable
        char tmps[255];
        strncpy(tmps, varnames->data+iVar*255, 255);
        std::string varname(tmps);
        varname = varname.substr(0,varname.find(" "));
        vdata->SetName(varname.c_str());
        // insert array of variable into the output
        output->GetPointData()->AddArray(vdata);
      }
   }
}

#if USE_PARTICLES
/*
 * This function inserts the data, loaded by the Posti tool, into a ouput
 */
/* void visuReader::InsertPartData(vtkPolyData* mb_part, int blockno, struct DoubleARRAY* coords, */
/*     struct DoubleARRAY* values, struct IntARRAY* nodeids, struct CharARRAY* varnames, struct IntARRAY* components) { */
void visuReader::InsertPartData(vtkMultiBlockDataSet* mb_part,int blockno
                               ,struct DoubleARRAY* coords   ,struct DoubleARRAY* values
                               ,struct IntARRAY* nodeids     ,struct CharARRAY* varnames
                               ,struct IntARRAY* components) {
   SWRITE("Insert particle data");
	 /* vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(mb_part->GetBlock(blockno)); */
   vtkSmartPointer<vtkPolyData> output = vtkPolyData::SafeDownCast(mb_part->GetBlock(blockno));

   // create points(array)
   vtkSmartPointer <vtkDoubleArray> pdata = vtkSmartPointer<vtkDoubleArray>::New();
	 pdata->SetNumberOfComponents(3); // 3D
   pdata->SetNumberOfTuples(coords->len/3);
   // copy coordinates
   double* ptr = pdata->GetPointer(0);
   for (long i = 0; i < coords->len; ++i)
   {
      *ptr++ = coords->data[i];
   }

	 // create points array to be used for the output
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
   points->SetData(pdata);
   output->SetPoints(points);

   vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
   vtkSmartPointer<vtkVertex>    vertex    = vtkSmartPointer<vtkVertex>::New();
   // Use the nodeids to build vertices (a cell that represents a 3D point)
   // (here we must copy the nodeids, we can not just assign the array of nodeids to some vtk-structure)
   // loop over all nodeids
   for (int i=0; i<nodeids->len; i++) {
      vertex   ->GetPointIds()->SetId(0,nodeids->data[i]);
      cellarray->InsertNextCell(vertex);
   }
   output->SetVerts(cellarray);

   // assign the actual data, loaded by the Posti tool, to the output
   unsigned int nVarCombine = varnames->len/255;
   unsigned int nVar = 0;
   for (unsigned int iVar=0;iVar<nVarCombine;iVar++) {
     nVar += components->data[iVar];
   }

	 if (nVar > 0) {
	   unsigned int sizePerVar = values->len/nVar;
	   int          dataPos    = 0;
	   // loop over all loaded variables
	   for (unsigned int iVar = 0; iVar < nVarCombine; iVar++) {
        vtkSmartPointer <vtkDoubleArray> vdata = vtkSmartPointer<vtkDoubleArray>::New();
        vdata->SetNumberOfComponents(components->data[iVar]);
        vdata->SetNumberOfTuples(sizePerVar);
        // copy coordinates
        double* ptr = vdata->GetPointer(0);
				for (long j = 0; j < sizePerVar; ++j)
				{
          for (long i = 0; i < components->data[iVar]; ++i)
          {
            *ptr++ = values->data[dataPos+i+j*nVar];
         // std::cout << "Data " << values->data[dataPos+i+j*nVar] << "\n";
          }
        }
        dataPos += components->data[iVar];
        // set name of variable
        char tmps[255];
        strncpy(tmps, varnames->data+iVar*255, 255);
        std::string varname(tmps);
        varname = varname.substr(0,varname.find(" "));
        // std::cout << "Varname " << varname << "\n";
        vdata->SetName(varname.c_str());
        // insert array of variable into the output
        output->GetPointData()->AddArray(vdata);
    }
  }
}
#endif /*USE_PARTICLES*/

#if USE_MPI
void visuReader::DistributeData(vtkMultiBlockDataSet* mb, int blockno) {
    // Check if there are Cells in the DataSet
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::SafeDownCast(mb->GetBlock(blockno));
    if (grid->GetNumberOfCells() == 0) {
      return;
    }

    int nGhosts = 1;

    // Apply D3 filter to create ghost cells
    SWRITE("Distributing data with minimum level of ghost cells : " << nGhosts);

    // Apply D3 filter to create ghost cells
    vtkPDistributedDataFilter * d3 = vtkPDistributedDataFilter::New();
    d3->SetInputData(mb->GetBlock(blockno));
    d3->SetMinimumGhostLevel(nGhosts);
    d3->SetUseMinimalMemory(false);
    d3->SetIncludeAllIntersectingCells(true);
    d3->SetRetainKdtree(true);
    d3->SetClipCells(false);

    // Information must be duplicated on all procs, ASSIGN_TO_ALL_INTERSECTING_REGIONS = 1
    d3->SetBoundaryMode(1);
    d3->Update();

    // Copy information pointer back
    mb->GetBlock(blockno)->ShallowCopy(d3->GetOutput());
    d3->Delete();
}
#endif /* USE_MPI */

visuReader::~visuReader(){
   SWRITE("~visuReader");
   delete [] FileName;
   this->VarDataArraySelection->Delete();
#if USE_PARTICLES
   this->VarParticleDataArraySelection->Delete();
#endif
}

/*
 * The following functions create the interaction between this Reader and
 * the gui, which is defined in the visu2DReader.xml
 * They return the number of available variables (state,primitive, derived)
 * and return the names of the variables, ....
 */

void visuReader::DisableAllVarArrays() {
   this->VarDataArraySelection->DisableAllArrays();
}

void visuReader::EnableAllVarArrays() {
   this->VarDataArraySelection->EnableAllArrays();
}

int visuReader::GetNumberOfVarArrays() {
   return this->VarDataArraySelection->GetNumberOfArrays();
}

const char* visuReader::GetVarArrayName(int index) {
   if (index >= ( int ) this->GetNumberOfVarArrays() || index < 0) {
      return NULL;
   }
   else {
      return this->VarDataArraySelection->GetArrayName(index);
   }
}

int visuReader::GetVarArrayStatus(const char* name) {
   return this->VarDataArraySelection->ArrayIsEnabled(name);
}

void visuReader::SetVarArrayStatus(const char* name, int status) {
   if (status) {
      this->VarDataArraySelection->EnableArray(name);
   }
   else {
      this->VarDataArraySelection->DisableArray(name);
   }
}

void visuReader::DisableAllBCArrays() {
   this->BCDataArraySelection->DisableAllArrays();
}

void visuReader::EnableAllBCArrays() {
   this->BCDataArraySelection->EnableAllArrays();
}

int visuReader::GetNumberOfBCArrays() {
   return this->BCDataArraySelection->GetNumberOfArrays();
}

const char* visuReader::GetBCArrayName(int index) {
   if (index >= ( int ) this->GetNumberOfBCArrays() || index < 0) {
      return NULL;
   }
   else {
      return this->BCDataArraySelection->GetArrayName(index);
   }
}

int visuReader::GetBCArrayStatus(const char* name) {
   return this->BCDataArraySelection->ArrayIsEnabled(name);
}

void visuReader::SetBCArrayStatus(const char* name, int status) {
   if (status) {
      this->BCDataArraySelection->EnableArray(name);
   }
   else {
      this->BCDataArraySelection->DisableArray(name);
   }
}

void visuReader::SelectionModifiedCallback(vtkObject*, unsigned long,
      void* clientdata, void*) {
   static_cast<visuReader*>(clientdata)->Modified();
}

int visuReader::FillOutputPortInformation(
      int vtkNotUsed(port), vtkInformation* info) {
   info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
   return 1;
}

#if USE_PARTICLES
void visuReader::DisableAllVarParticleArrays() {
   this->VarParticleDataArraySelection->DisableAllArrays();
}

void visuReader::EnableAllVarParticleArrays() {
   this->VarParticleDataArraySelection->EnableAllArrays();
}

int visuReader::GetNumberOfVarParticleArrays() {
   return this->VarParticleDataArraySelection->GetNumberOfArrays();
}

const char* visuReader::GetVarParticleArrayName(int index) {
   if (index >= ( int ) this->GetNumberOfVarParticleArrays() || index < 0) {
      return NULL;
   }
   else {
      return this->VarParticleDataArraySelection->GetArrayName(index);
   }
}

int visuReader::GetVarParticleArrayStatus(const char* name) {
   return this->VarParticleDataArraySelection->ArrayIsEnabled(name);
}

void visuReader::SetVarParticleArrayStatus(const char* name, int status) {
   if (status) {
      this->VarParticleDataArraySelection->EnableArray(name);
   }
   else {
      this->VarParticleDataArraySelection->DisableArray(name);
   }
}
#endif /*USE_PARTICLES*/
