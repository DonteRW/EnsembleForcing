//**********************************************************************************************
//
//
//********************************************************************************************** 

#include "ens_land.h"
#include <iterator>

//2.25.18 writes vector to a 2D nc at time step    //--without-mpi
//__device__ __host__
//int Write1DVector_to2DNC(const char* FileName, const char* varName, int tIndx, int Nz_dim, float* var_inp) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, v_varid = 0;
//
//	const int  NDIMS = 2;
//	// The start and count arrays will tell the netCDF library where to write our data.    
//	size_t start[NDIMS], count[NDIMS];
//	// Open netcdf file.  
//	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))             //| NC_MPIIO, inpComm, inpInfo
//		ERR(retncval);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, varName, &v_varid)))
//		ERR(retncval);
//
//	start[0] = tIndx;	//0;
//	start[1] = 0;   // z_dim
//	count[0] = 1;   // Nt_dim;
//	count[1] = Nz_dim;
//
//	//put variable values
//	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0]))
//		ERR(retncval);
//
//	//close files
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//
//	//std::cout << "Sucess storing dimension vars in: " << FileName << std::endl;
//	//fflush(stdout); 
//	return 0;
//}
////--without-mpi
////1.25.18 writes a 2D array at time step
////__device__ __host__
//int Write2DSlub_to3DNC(const char* FileName, const char* varName, int tIndx, int Np_dim, int Nz_dim, float** var_inp) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, v_varid = 0;
//
//	const int  NDIMS = 3;
//	// The start and count arrays will tell the netCDF library where to write our data.    
//	size_t start[NDIMS], count[NDIMS];
//	// Open netcdf file.  
//	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))             //| NC_MPIIO, inpComm, inpInfo
//		ERR(retncval);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, varName, &v_varid)))
//		ERR(retncval);
//
//	start[0] = tIndx;	//0;
//	start[1] = 0;   // point
//	start[2] = 0;   // z_dim;
//	count[0] = 1;   // Nt_dim;
//	count[1] = Np_dim;
//	count[2] = Nz_dim;
//
//	//put variable values
//	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0][0]))
//		ERR(retncval);
//
//	//close files
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//
//	//std::cout << "Sucess storing dimension vars in: " << FileName << std::endl;
//	//fflush(stdout); 
//	return 0;
//}
////no-mpi: 2.25.18 for ens and da outputs 
////creates 2D netcdf and stores dimension variables; called once for a given output netcdf 
////__device__ __host__
//int Create2DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
//	const char* zName, int tDim, int zDim, float* t_inp, float* fillVal)           // MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	// IDs for the netCDF file, dimensions, and variables. 
//	int ncid = 0, t_dimid = 0, z_dimid = 0;
//	int v_varid = 0, t_varid = 0;
//	const char* vunits = "UNITS";
//	const char* fillValue = "_FillValue";
//	//float missVal = -9999;
//	const int  NDIMS = 2;
//	int dimids[NDIMS];
//	int oldFill = 0;
//	// The start and count arrays will tell the netCDF library where to write our data.    
//	size_t start[NDIMS], count[NDIMS];
//	// Error handling.  
//	int retncval = 0;
//	// Create netcdf file.  
//	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
//		ERR(retncval);
//	//set fill on
//	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
//		ERR(retncval);
//	/* Define the dimensions. record dim can be unlimited*/
//	if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
//		ERR(retncval);
//	if ((retncval = nc_def_dim(ncid, zName, zDim, &z_dimid)))
//		ERR(retncval);
//	/* The dimids array is used to pass the dimids of the dimensions of
//	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
//	dimids[0] = t_dimid;
//	dimids[1] = z_dimid;
//	// Define the netCDF variables 
//	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))
//		ERR(retncval);
//	//assign missing
//	if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
//		ERR(retncval);
//	// Assign units attributes to the netCDF variables.  
//	if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
//		ERR(retncval);
//	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
//		ERR(retncval);
//	// Assign units attributes to the netCDF variables.  
//	if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
//		ERR(retncval);
//	if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
//		ERR(retncval);
//	if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
//		ERR(retncval);
//
//	//put values to dim variables
//	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0])))
//		ERR(retncval);
//
//	//close file
//	if ((retncval = nc_close(ncid)))
//		ERR(retncval);
//	//delte 3D array	
//	std::cout << "Sucess creating and storing dimension vars in: " << FileName << std::endl;
//	//fflush(stdout); 
//	return 0;
//}
////no-mpi: 1.18.18 for ens and da outputs 
////creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
////__device__ __host__
//int Create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
//	const char* yxName, const char* zName, int tDim, int yxDim, int zDim, float* t_inp, float* fillVal)           // MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	// IDs for the netCDF file, dimensions, and variables. 
//	int ncid = 0, t_dimid = 0, yx_dimid = 0, z_dimid = 0;
//	int v_varid = 0, t_varid = 0;
//	const char* vunits = "UNITS";
//	const char* fillValue = "_FillValue";
//	//float missVal = -9999;
//	const int  NDIMS = 3;
//	int dimids[NDIMS];
//	int oldFill = 0;
//	// The start and count arrays will tell the netCDF library where to write our data.    
//	size_t start[NDIMS], count[NDIMS];
//	// Error handling.  
//	int retncval = 0;
//	// Create netcdf file.  
//	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
//		ERR(retncval);
//	//set fill on
//	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
//		ERR(retncval);
//	/* Define the dimensions. record dim can be unlimited*/
//	if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
//		ERR(retncval);
//	if ((retncval = nc_def_dim(ncid, yxName, yxDim, &yx_dimid)))
//		ERR(retncval);
//	if ((retncval = nc_def_dim(ncid, zName, zDim, &z_dimid)))
//		ERR(retncval);
//	/* The dimids array is used to pass the dimids of the dimensions of
//	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
//	dimids[0] = t_dimid;
//	dimids[1] = yx_dimid;
//	dimids[2] = z_dimid;
//	// Define the netCDF variables 
//	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))
//		ERR(retncval);
//	//assign missing
//	if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
//		ERR(retncval);
//	// Assign units attributes to the netCDF variables.  
//	if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
//		ERR(retncval);
//	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
//		ERR(retncval);
//	// Assign units attributes to the netCDF variables.  
//	if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
//		ERR(retncval);
//	if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
//		ERR(retncval);
//	if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
//		ERR(retncval);
//
//	//put values to dim variables
//	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0])))
//		ERR(retncval);
//
//	//close file
//	if ((retncval = nc_close(ncid)))
//		ERR(retncval);
//	//delte 3D array	
//	std::cout << "Sucess creating and storing dimension vars in: " << FileName << std::endl;
//	//fflush(stdout); 
//	return 0;
//}
//
////1.18.18 for forcing at obs points
////  read wole matrix at once
////__device__ __host__
//int read2DNC_Contigious(const char* FILE_NAME, const char* VAR_NAME, double ** &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
//	//Open the netcdf file.  
//	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
//		ERR(retncval);  // , FILE_NAME);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
//		ERR(retncval);
//
//	//read variable (input data)	
//	if (retncval = nc_get_var_double(ncid, pvarid, &pvar_in[0][0]))
//		ERR(retncval);
//
//	//close netcdf file			
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//
//	return 0;
//}
////10.6.17 read multiple arrays (vectors of vals at points) for given time range
////__device__ __host__
//int readNC_vector(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, int tEnd, float** &pvar_in, int &nrecords)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//float* pvarin_temp = NULL;
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
//											//variable data type
//	nc_type varType;
//	size_t pdim_sizes;
//	//array of dimensions
//	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
//					//dimension names 
//	char pdim_Names[80];
//	size_t start[2], count[2];
//	//Open the netcdf file.  
//	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
//		ERR(retncval);   // , FILE_NAME);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
//		ERR(retncval);
//	//var information, checking the dimension array
//	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
//		ERR(retncval);
//
//	//check dimension info and set start and count arrays; 
//	for (int i = 0; i < 2; i++) {
//		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
//			ERR(retncval);
//		if (strcmp(pdim_Names, tcor_NAME) == 0) {
//			start[i] = tStart;
//			if (tEnd < pdim_sizes) {
//				count[i] = tEnd - tStart;
//				tStart += count[i];                //new start point
//			}
//			else {
//				count[i] = pdim_sizes - tStart;        //take the lower of the tEnd/count[i] to guarantee against going out of bounds
//													   //numNc++;                               //next time go to the next netcdf file
//				tStart = 0;
//			}
//			nrecords = count[i];
//		}
//		else {
//			start[i] = 0;
//			count[i] = pdim_sizes;
//			//yxDim *= count[i];
//		}
//	}
//	/*if (pvar_in != NULL)
//	delete3DArrayblock_Contiguous(pvar_in);
//	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
//	//read var data
//	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0][0]))
//		ERR(retncval);
//
//	//close netcdf file			
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//	return 0;
//}
////2.19.18 read array at an index (time step)
////__device__ __host__
//int readNC_Array_atIndex(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, float* pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//float* pvarin_temp = NULL;
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
//											//variable data type
//	nc_type varType;
//	size_t pdim_sizes;
//	//array of dimensions
//	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
//					//dimension names 
//	char pdim_Names[80];
//	size_t start[2], count[2];
//	//Open the netcdf file.  
//	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
//		ERR(retncval);    // , FILE_NAME);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
//		ERR(retncval);
//	//var information, checking the dimension array
//	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
//		ERR(retncval);
//
//	//check dimension info and set start and count arrays; 
//	for (int i = 0; i < 2; i++) {
//		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
//			ERR(retncval);
//		if (strcmp(pdim_Names, tcor_NAME) == 0) {
//			start[i] = tIndex;
//			count[i] = 1;
//		}
//		else {
//			start[i] = 0;
//			count[i] = pdim_sizes;
//			//yxDim *= count[i];
//		}
//	}
//	/*if (pvar_in != NULL)
//	delete3DArrayblock_Contiguous(pvar_in);
//	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
//	//read var data
//	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0]))
//		ERR(retncval);
//
//	//close netcdf file			
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//	return 0;
//}
////for single point //2.19.18 read array at an index (time step)
////__device__ __host__
//int readNC_Array_atIndices(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, int pIndex, float* &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//{
//	//float* pvarin_temp = NULL;
//	//ids for variable, axes,...
//	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
//											//variable data type
//	nc_type varType;
//	size_t pdim_sizes;
//	//array of dimensions
//	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
//					//dimension names 
//	char pdim_Names[80];
//	size_t start[2], count[2];
//	//Open the netcdf file.  
//	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
//		ERR(retncval);  // , FILE_NAME);
//	// get variable id
//	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
//		ERR(retncval);
//	//var information, checking the dimension array
//	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
//		ERR(retncval);
//
//	//check dimension info and set start and count arrays; 
//	for (int i = 0; i < 2; i++) {
//		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
//			ERR(retncval);
//		if (strcmp(pdim_Names, tcor_NAME) == 0) {
//			start[i] = tIndex;
//			count[i] = 1;
//		}
//		else {
//			start[i] = pIndex;
//			count[i] = 1;    // pdim_sizes;
//			//yxDim *= count[i];
//		}
//	}
//	/*if (pvar_in != NULL)
//	delete3DArrayblock_Contiguous(pvar_in);
//	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
//	//read var data
//	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0]))
//		ERR(retncval);
//
//	//close netcdf file			
//	if (retncval = nc_close(ncid))
//		ERR(retncval);
//	return 0;
//}

int ReadVectorLength(const char* static_file_name, const char* VAR_NAME, int& vector_length)
{
	
	int retncval = 0, ncid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
												//variable data type
	//nc_type varType;
	//size_t pdim_sizes;
	//array of dimensions
	int pdimid; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
					//dimension names 
	//char pdim_Names[80];
	size_t i_len;  // , start[2], count[2];
	//Open the netcdf file.  
	if ((retncval = nc_open(static_file_name, NC_NOWRITE, &ncid)))
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_dimid(ncid, VAR_NAME, &pdimid)))
		ERR(retncval);
	if ((retncval = nc_inq_dimlen(ncid, pdimid, &i_len)))
		ERR(retncval);
	vector_length = int(i_len);

	//close netcdf file			
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	return 0;
}
int read_latlon_vec(const char* FILE_NAME, const char* lat_NAME, const char* lon_NAME,
	int v_len, std::vector<float> &RLA_land, std::vector<float> &RLO_land)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);  // , FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, lat_NAME, &pvarid)))
		ERR(retncval);
	float* pvar_in = new float[v_len];
		//read variable (input data)	
	if ((retncval = nc_get_var_float(ncid, pvarid, &pvar_in[0])))
		ERR(retncval);
	RLA_land.assign(pvar_in, pvar_in + v_len);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, lon_NAME, &pvarid)))
		ERR(retncval);

	//read variable (input data)	
	if ((retncval = nc_get_var_float(ncid, pvarid, &pvar_in[0])))
		ERR(retncval);
	RLO_land.assign(pvar_in, pvar_in + v_len);

	//close netcdf file			
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	delete[] pvar_in;
	pvar_in = NULL;

	return 0;
}

int read1DNC(const char* FILE_NAME, const char* VAR_NAME, double* &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);  // , FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	//read variable (input data)	
	if ((retncval = nc_get_var_double(ncid, pvarid, &pvar_in[0])))
		ERR(retncval);

	//close netcdf file			
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	return 0;
}

int read2DNC(const char* FILE_NAME, const char* VAR_NAME, double ** &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);  // , FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);

	//read variable (input data)	
	if ((retncval = nc_get_var_double(ncid, pvarid, &pvar_in[0][0])))
		ERR(retncval);

	//close netcdf file			
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	return 0;
}

int write_random_samples(int grid_size, int ens_size, int n_vars,
	std::vector<Eigen::Matrix<float, Dynamic, Dynamic>>  //, RowMajor
	multivarNormalDistSamplesForc,
	std::string FileName, std::string var_names[8], float* fillVal)
			 // MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, yx_dimid = 0, z_dimid = 0;   // , t_dimid = 0, 
	int v_varid = 0;     // , t_varid = 0;
	//const char* vunits = "UNITS";
	const char* fillValue = "_FillValue";

	int retncval = 0;

	float** ensForcingMultiplier = create2DArray_Contiguous(grid_size, ens_size); // new float*[grid_size];
	/*for (int id = 0; id < grid_size; id++)
	{
		ensForcingMultiplier[id] = new float[ens_size];
	}
	float* ensForcingMultiplier = new float[grid_size];*/
	

	const int  NDIMS = 2;
	int dimids[NDIMS];
	int oldFill = 0;
	
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	start[0] = 0;
	start[1] = 0;
	count[0] = 1;  // grid_size;
	count[1] = ens_size;

	// Create netcdf file.  
	if ((retncval = nc_create((const char*)FileName.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid)))
	{
		std::cout << " create file " << FileName;
		ERR(retncval);
	}
	//set fill on
	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
	{
		std::cout << " set fill " << FileName;
		ERR(retncval);
	}
	/* Define the dimensions. record dim can be unlimited*/
	/*if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
	ERR(retncval);*/
	if ((retncval = nc_def_dim(ncid, "location", grid_size, &yx_dimid)))
	{
		std::cout << " create dim location " << FileName;
		ERR(retncval);
	}
	if ((retncval = nc_def_dim(ncid, "ensemble_member", ens_size, &z_dimid)))
	{
		std::cout << " create dim ens " << FileName;
		ERR(retncval);
	}
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	dimids[0] = yx_dimid;
	dimids[1] = z_dimid;
	//std::cout << " temp anom " << FileName;
	//std::cout << multivarNormalDistSamplesForc[0];
	for (int is = 0; is < n_vars; is++) {
		// Define the netCDF variables 
		if ((retncval = nc_def_var(ncid, (const char*)var_names[is].c_str(), 
			                       NC_FLOAT, NDIMS, dimids, &v_varid)))
		{
			std::cout << " define var " << var_names[is];
			ERR(retncval);
		}
		//assign missing
		if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
		{
			std::cout << " put att fill " << var_names[is];
			ERR(retncval);
		}
		// Assign units attributes to the netCDF variables.  
		/*if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
		ERR(retncval);
		if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
		ERR(retncval);*/
		// Assign units attributes to the netCDF variables.  
		/*if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
		ERR(retncval);
		if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
		ERR(retncval);
		if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
		ERR(retncval);*/
		for (int id = 0; id < grid_size; id++)
		{
			for (int ie = 0; ie < ens_size; ie++) {
				ensForcingMultiplier[id][ie] = multivarNormalDistSamplesForc[is](id, ie);
			}
			//start[0] = id;
			//start[1] = 0;
			//count[0] = 1;  // grid_size;
			//count[1] = ens_size;
		}
		//put values tovariables
		if ((retncval = nc_put_var_float(ncid, v_varid, &ensForcingMultiplier[0][0])))
		{
			std::cout << " put variable " << var_names[is];
			ERR(retncval);
		}
		
	}

	//close file
	if ((retncval = nc_close(ncid)))
	{
		std::cout << " close file " << FileName;
		ERR(retncval);
	}
	//delte 3D array	
	std::cout << "Sucess creating and writing " << FileName << std::endl;
	//fflush(stdout); 

	delete2DArray_Contiguous(ensForcingMultiplier);
	/*for (int id = 0; id < grid_size; id++)
	{
		delete[] ensForcingMultiplier[id];
		ensForcingMultiplier[id] = NULL;

	}
	delete[] ensForcingMultiplier;
	ensForcingMultiplier = NULL;*/

	return 0;
}

int Read_random_samples(int grid_size, int ens_size, int n_vars,
	std::vector<Eigen::Matrix<float, Dynamic, Dynamic>> //, RowMajor
	&multivarNormalDistSamplesForc_in,
	std::string FileName, std::string var_names[8])   //, const char* yxName, float* fillVal)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid;   // = 0;   // , yx_dimid = 0, z_dimid = 0;   // , t_dimid = 0, 
	int v_varid;   // = 0;     // , t_varid = 0;
	// Error handling.  
	int retncval = 0;
	const int NDIMS = 2;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	start[0] = 0;
	start[1] = 0;
	count[0] = grid_size;
	count[1] = ens_size;

	float** ensForcingMultiplier = create2DArray_Contiguous(grid_size, ens_size); //new float*[grid_size];
	/*for (int id = 0; id < grid_size; id++)
	{
		ensForcingMultiplier[id] = new float[ens_size];
	}*/

	//Open the netcdf file.  
	if ((retncval = nc_open((const char*)FileName.c_str(), NC_NOWRITE, &ncid)))
	{
		std::cout << " open file " << FileName;
		ERR(retncval);
	}
	
	for (int is = 0; is < n_vars; is++) {
		// get variable id
		if ((retncval = nc_inq_varid(ncid, (const char*)var_names[is].c_str(), &v_varid)))
		{
			std::cout << " inq variable " << var_names[is];
			ERR(retncval);
		}

		//read variable (input data)	
		if ((retncval = nc_get_var_float(ncid, v_varid, &ensForcingMultiplier[0][0])))
		{
			std::cout << " get variable " << var_names[is];
			ERR(retncval);
		}

		for (int id = 0; id < grid_size; id++) 
		{
			for (int ie = 0; ie < ens_size; ie++) {
				multivarNormalDistSamplesForc_in[is](id, ie) = ensForcingMultiplier[id][ie];
			}
		}
	}
	//close file
	if ((retncval = nc_close(ncid)))
	{
		std::cout << " close file " << FileName;
		ERR(retncval);
	}

	std::cout << "Sucess reading " << FileName << std::endl;
	//fflush(stdout); 
	delete2DArray_Contiguous(ensForcingMultiplier);
	/*for (int id = 0; id < grid_size; id++)
	{
		delete[] ensForcingMultiplier[id];
		ensForcingMultiplier[id] = NULL;
	}
	delete[] ensForcingMultiplier;
	ensForcingMultiplier = NULL;*/

	return 0;
}

int write_random_samples(int grid_size, int ens_size,
	Eigen::Matrix<float, Dynamic, Dynamic> multivarNormalDistSamplesForc, //, RowMajor
	const char* FileName, const char* VarName, const char* yxName, float* fillVal)
	//const char *varUnits, 
	//const char* tName, int tDim, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
			 // MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, yx_dimid = 0, z_dimid = 0;   // , t_dimid = 0, 
	int v_varid = 0;     // , t_varid = 0;
	//const char* vunits = "UNITS";
	const char* fillValue = "_FillValue";

	float** ensForcingMultiplier = new float*[ens_size];
	for (int ie = 0; ie < ens_size; ie++)
	{
		ensForcingMultiplier[ie] = new float[grid_size];
	}
	//float missVal = -9999;
	const int  NDIMS = 2;
	int dimids[NDIMS];
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	//size_t start[NDIMS], count[NDIMS];

	// Error handling.  
	int retncval = 0;
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
		ERR(retncval);
	//set fill on
	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
		ERR(retncval);
	/* Define the dimensions. record dim can be unlimited*/
	/*if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
	ERR(retncval);*/
	if ((retncval = nc_def_dim(ncid, yxName, grid_size, &yx_dimid)))
		ERR(retncval);
	if ((retncval = nc_def_dim(ncid, "ensemble_member", ens_size, &z_dimid)))
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */

	dimids[0] = z_dimid;
	dimids[1] = yx_dimid;
	// Define the netCDF variables 
	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))
		ERR(retncval);
	//assign missing
	if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	/*if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
	ERR(retncval);
	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
	ERR(retncval);*/
	// Assign units attributes to the netCDF variables.  
	/*if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
	ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
	ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
	ERR(retncval);*/
	for (int ie = 0; ie < ens_size; ie++)
	{
		for (int id = 0; id < grid_size; id++)
		{
			ensForcingMultiplier[ie][id] = multivarNormalDistSamplesForc(id, ie);
		}
	}
	//put values tovariables
	if ((retncval = nc_put_var_float(ncid, v_varid, &ensForcingMultiplier[0][0])))
		ERR(retncval);

	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);
	//delte 3D array	
	std::cout << "Sucess creating and writing " << FileName << std::endl;
	//fflush(stdout); 

	for (int ie = 0; ie < ens_size; ie++)
	{
		delete[] ensForcingMultiplier[ie];
		ensForcingMultiplier[ie] = NULL;

	}
	delete[] ensForcingMultiplier;
	ensForcingMultiplier = NULL;

	return 0;
}

int Read_random_samples(int grid_size, int ens_size,
	Eigen::Matrix<float, Dynamic, Dynamic> &multivarNormalDistSamplesForc, //,RowMajor
	const char* FileName, const char* VarName)   //, const char* yxName, float* fillVal)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid;   // = 0;   // , yx_dimid = 0, z_dimid = 0;   // , t_dimid = 0, 
	int v_varid;   // = 0;     // , t_varid = 0;

	float** ensForcingMultiplier = new float*[ens_size];
	for (int ie = 0; ie < ens_size; ie++)
	{
		ensForcingMultiplier[ie] = new float[grid_size];
	}

	// Error handling.  
	int retncval = 0;


	//Open the netcdf file.  
	if ((retncval = nc_open(FileName, NC_NOWRITE, &ncid)))
		ERR(retncval);  // , FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VarName, &v_varid)))
		ERR(retncval);

	//read variable (input data)	
	if ((retncval = nc_get_var_float(ncid, v_varid, &ensForcingMultiplier[0][0])))
		ERR(retncval);

	for (int ie = 0; ie < ens_size; ie++)
	{
		for (int id = 0; id < grid_size; id++)
		{
			multivarNormalDistSamplesForc(id, ie) = ensForcingMultiplier[ie][id];
		}
	}

	//close file
	if ((retncval = nc_close(ncid)))
	{
		std::cout << " close file " << FileName;
		ERR(retncval);
	}

	std::cout << "Sucess reading " << FileName << std::endl;
	//fflush(stdout); 
	for (int ie = 0; ie < ens_size; ie++)
	{
		delete[] ensForcingMultiplier[ie];
		ensForcingMultiplier[ie] = NULL;

	}
	delete[] ensForcingMultiplier;
	ensForcingMultiplier = NULL;

	return 0;
}

int ReadTileInfo(int myrank, const char* filename, int vector_length, int &beg_indx, int &end_indx,
	//std::vector<int> tile_xy, std::vector<int> Idim_xy, std::vector<int> Jdim_xy,
	//std::vector<double> OROG, std::vector<double> VETFCS,
	std::vector<float> &RLA, std::vector<float> &RLO) //! vegetation_category)
{

	int ncid, dimid, varid;
	// Error handling.  
	int retncval = 0;
	size_t vector_length_in;

	if ((retncval = nc_open(filename, NC_NOWRITE, &ncid)))
		ERR(retncval);

	if ((retncval = nc_inq_dimid(ncid, "location", &dimid)))
		ERR(retncval);
	if ((retncval = nc_inq_dimlen(ncid, dimid, &vector_length_in)))
		ERR(retncval);

	if(vector_length != int(vector_length_in))
	{
		std::cout << "wrong vector size, stop" << std::endl;
		exit(1);
	}
	int *tile_xy_in = new int[vector_length_in];

	if ((retncval = nc_inq_varid(ncid, "cube_tile", &varid)))
	// (status /= nf90_noerr) call handle_err(status)
	   ERR(retncval);     //(status, 'reading cube_tile var id')
	if ((retncval = nc_get_var(ncid, varid, tile_xy_in)))
		ERR(retncval);     //(status, 'reading cube tile value')

	/*auto itr = std::find(tile_xy_in, tile_xy_in + vector_length_in, indx);
	if (itr != std::end(tile_xy_in)
		indx = std::distance(tile_xy_in, itr);*/
	std::vector<int> vec(tile_xy_in, tile_xy_in + vector_length_in);

	std::vector<int>::iterator itr = std::find(vec.begin(), vec.end(), myrank+1);
	if (itr != vec.cend())
		beg_indx = std::distance(vec.begin(), itr);
	else
	{
		std::cout << "error proc "<< myrank << " vector index not found; exit" << std::endl;
		MPI::Finalize();
		return 1;
	}
	itr = std::find(vec.begin(), vec.end(), myrank + 2);
	if (itr != vec.cend())
		end_indx = std::distance(vec.begin(), itr) - 1;
	else
		end_indx = vec.size() - 1;

	float *latlon_Array = new float [vector_length_in];

	if ((retncval = nc_inq_varid(ncid, "latitude", &varid)))
		ERR(retncval);
	if ((retncval = nc_get_var(ncid, varid, latlon_Array)))
		ERR(retncval);

	RLA.assign(latlon_Array+beg_indx, latlon_Array+ end_indx+1);

	if ((retncval = nc_inq_varid(ncid, "longitude", &varid)))
		ERR(retncval);
	if ((retncval = nc_get_var(ncid, varid, latlon_Array)))
		ERR(retncval);    
	RLO.assign(latlon_Array + beg_indx, latlon_Array + end_indx + 1);

	/*if ((retncval = nc_inq_varid(ncid, "elevation", &varid)
	if ((retncval = nc_get_var(ncid, &varid, OROG, &
		start = (/ 1 / ), count = (/ vector_length / ))
		ERR(retncval);    

	if ((retncval = nc_inq_varid(ncid, "cube_i", &varid)
		ERR(retncval);     
	if ((retncval = nc_get_var(ncid, &varid, Idim_xy, &
		start = (/ 1 / ), count = (/ vector_length / ))
		ERR(retncval);   

	if ((retncval = nc_inq_varid(ncid, "cube_j", &varid)
		ERR(retncval);     
	if ((retncval = nc_get_var(ncid, &varid, Jdim_xy, &
		start = (/ 1 / ), count = (/ vector_length / ))
		ERR(retncval);    

	if ((retncval = nc_inq_varid(ncid, "vegetation_category", &varid)
		ERR(retncval);     
	if ((retncval = nc_get_var(ncid, &varid, VETFCS, &
		start = (/ 1 / ), count = (/ vector_length / ))
		ERR(retncval);     

	!if ((retncval = nc_inq_varid(ncid, "land_mask", &varid)
	!if ((retncval = nc_get_var(ncid, &varid, this%land_mask, &
		!start = (/ namelist % begsub / ), count = (/ namelist % lensub / ))
		ERR(retncval);*/

	if ((retncval = nc_close(ncid)))
		ERR(retncval);     //(status, 'closing file')

	delete[] tile_xy_in;
	tile_xy_in = NULL;
	delete[] latlon_Array;
	latlon_Array = NULL;

	return 0;
}

float*** create3DArrayblock_Contiguous(int nt, int nr, int nc)      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
{
	float*** myMatrix = new float**[nt];
	float* ptrMemory = new float[nt*nr*nc];
	for (int t = 0; t < nt; t++)
	{
		//this looks suspicious 
		//#??_TBC 6.18.13
		myMatrix[t] = new float*[nr];
		for (int r = 0; r < nr; r++)
		{
			myMatrix[t][r] = ptrMemory; //new float*[nr];  
			ptrMemory += nc;
		}
	}
	return myMatrix;

}//float*** create3DArrayblock

//delets a 3D array (frees memory) allocated contiguously
//__host__ __device__ 
void delete3DArrayblock_Contiguous(float*** myMatrix) // int nr, int nc)// input: 3D array
{
	/*for (int t = 0; t< nt; t++)
	{
	/*for (int r=0; r<nr; r++)
	delete [] myMatrix[t][r];

	}*/
	delete[] myMatrix[0];
	delete[] myMatrix;
	return;

}//double** CreateMatrix
//******Warning the following code seems to have memo leak------ 5.7.15 
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
//__host__ __device__ 
float** create2DArray_Contiguous(int nr, int nc)      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
{
	float** myMatrix = new float*[nr];
	float* ptrMemory = new float[nr*nc];
	for (int r = 0; r < nr; r++)
	{
		myMatrix[r] = ptrMemory; //new float*[nr];  
		ptrMemory += nc;
	}
	return myMatrix;
}//float*** create3DArrayblock

//delets a 2D array (frees memory allocated) contiguously
//__host__ __device__ 
void delete2DArray_Contiguous(float** myMatrix) // int nr, int nc)// input: 2D array
{
	delete[] myMatrix[0];
	delete[] myMatrix;

	return;
}//double** 

//Creates a matrix. The inputs are matrix dimensions. It allocates a memory block of size nrows*ncols* (size of float)
//and returns an array of pointers to the allocated memory block
//__host__ __device__ 
float*** Create3DArray(int nt, int nr, int nc)       //inputs: no. of rows and no. of colos (dimensions) of matrix
{
	float*** myMatrix = new float**[nt];
	for (int i = 0; i < nt; i++)
	{
		myMatrix[i] = new float*[nr];
		for (int j = 0; j < nr; j++)
			myMatrix[i][j] = new float[nc];
	}
	return myMatrix;

}//float** CreateMatrix

//The following program deletes a matrix passed to it (it frees up the memory block allocated to the matrix)
//__host__ __device__ 
void Delete3DArray(float ***A, int nt, int nr, int nc)            //input: A matrix
{
	for (int i = 0; i < nt; i++)
	{
		for (int j = 0; j < nr; j++)
			delete[] A[i][j];
		delete[] A[i];
	}
	delete[] A;

	return;
}//void DeleteMatrix

