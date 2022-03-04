#pragma once
#ifndef ENSLAND_H
#define ENSLAND_H

#include <ctime>
#include "Eigen/Dense"
#include<iomanip>
#include<iostream>
#include <fstream>
#include <cstring>
#include<sstream>
//#include <stdio.h>
#include "mpi.h"
#include <netcdf.h>
#include <netcdf_par.h>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <random>

using namespace Eigen;

#include "eigenmvn.h"
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"
//#pragma warning(disable : 4996)
//#include <algorithm>
//using namespace std;
//error handling for netcdf

#define  ERR(e) {std::cout<<"Error: "<< nc_strerror(e)<<std::endl; return 2; }

//int ERR(int e) { 
//	std::cout << "Error: " << nc_strerror(e) << std::endl; 
//	return 2; 
//}
//error handling for cuda
//define  cuda_checkERR(err) { if (err != cudaSuccess) std::cout << "Error: "<<cudaGetErrorString(err)<< std::endl; exit(EXIT_FAILURE); }
////__device__ __host__
//void cuda_checkERR(cudaError_t err);

struct ncdaOutput {
	char outfName[256];
	char symbol[256];
	char units[256];
	int outvarIndex;
};

class modelGridCells {

public:
	//default contr.
	modelGridCells();
	modelGridCells(int modGridCells, int numObsPoints, int ensSize, 
		int num_tim_cyc, int ini_tim, float modDT, float pert_facts[7]);

	~modelGridCells()
	{
		//clean up
	};

	//da setttings
	double daTime;

	int ns_statSize;       // for now 1 state assimilated

	int mod_gridSize;     //model domain number of grids
	int num_obs_points;
	int tot_gridSize;	  // total number of grid cells includeing obs points
	
	int ens_size;
	int particle_size;         
	// particle size = enssize * es_pfSize_Fact

	float spat_corLength;     //correlation length for forcing
	float temp_decorrLength;     //temporal decorrelation length (hrs, dt units)

	float obsErrStdev;		  //observed (state or equivalent) var Standard deviation
	float forc_stddev;

	float modelDT;
	int num_time_cycle;  //number of time steps in (daily) cycle
	int ini_time;  //very initial time (no stored random pattern)

	float corrFact1;   // = 1.0 - modelDT / decorLength;
	float corrFact2;	// = sqrtf(1.0 - corrFact1 * corrFact1);    //1.0 - corrFact1;  // 

	std::vector<float> ensForcingMultiplier;

	//for ens run with historical forcing
	int ModelStartDateEns[4];
	int ModelEndDateEns[4];

	/*//__host__ __device__*/
	void initDAMatrices(int myrank, std::vector<float> lat_arry, std::vector<float> lon_arry,
		    double &interm_time, bool use_seed, int seed_in = 0);

	void Generate_Samples(int myrank, double &interm_time, int istep = 0);

	//void Generate_Samples_and_Cov(int myrank, int istep, double &interm_time);
	//	// std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>> 
	//	//&multivarNormalDistSamplesForc);
	//void Create_Ensemble_Members(int myrank, int beg_indx,
	//	std::string FileDir, std::string forc_inp_file, double &interm_time);

	void Generate_Samples_and_Cov(int myrank, int istep,
		int num_samples, std::string FileName, double &interm_time);

	void Create_Ensemble_Samples(int myrank, int num_samples,
		std::string FileDir, std::string forc_inp_file, double &interm_time);

	

	//private:
	Eigen::EigenMultivariateNormal<float> std_norm_dist_Forc_Default;
	//Eigen::EigenMultivariateNormal<float> norm_dist_0Mean_Tempr;     //for temperature
	/*Eigen::EigenMultivariateNormal<float> norm_dist_0Mean_Default;
	Eigen::EigenMultivariateNormal<float> norm_dist_1Mean_Default;*/

	//Get choleski L to apply forcing correlation 
	Eigen::Matrix<float, Dynamic, Dynamic> corrFs_L;  //, RowMajor

	std::vector<Eigen::Matrix<float, Dynamic, Dynamic>>  //, RowMajor
		_multivarNormalDistSamplesForc;

	//Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> multivarNormalDistSamplesForc;
	
	/*Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> NormalDistSamplesForc;
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> NormalDistSamplesForc_mean1;*/
	//Eigen::VectorXf Z_obs;
	//Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> R_obsErrCov;
	//Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> P_stateCov;
	//std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > P_stateCovBackground;

	void ncERR(int e, std::string err_message) {
		std::cout << "Error: " << nc_strerror(e) <<" "<< err_message <<  std::endl;
		MPI::COMM_WORLD.Abort(2);
	}

	//3.20.18 ens random multipliers copied as class members
	//__host__ __device__
	//void runEnsembles(int nEns) {};
	///*//__host__ __device__*/
 //   //set background states for next time step to filter updated
	////__host__ __device__
	//void updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr) {};
	////4.11.18 update by weight
	////__host__ __device__
	//void updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr, std::vector<int> weightIndices) {};

	////__host__ __device__
	//void setNextStepStates(int nEns) {};

	////********UPDATEtime ()  Update time for each time step
	////__host__ __device__
	//void  UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT) {};

private:
	//int inputforcIndex;    // the forcing to perturb for ensembles--using the string lists below
	std::string var_names[7] = {"precipitation", "solar_radiation", "longwave_radiation",
		"temperature", "wind_speed", "specific_humidity", "surface_pressure" }; //"precipitation_bilinear","precipitation_conserve", not perturbed atm

	const char* var_units[7] = { "mm/s", "W/m2", "W/m2", "K", "m/s", "kg/kg", "Pa" };
	// ????? only do ens for precip, temp, SWRad, LWRad ?
	float std_dev[7] = { 0.1, 0.05, 20.0, 2.0, 0.05, 0.05, 0.01 };
	float mean_var[7] = { 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
	const int n_vars = 7;
	//bool multi_bool[8] = { "True", "True", "False", "False", "True", "True", "True", "True" };
	int multi_bool[7] = { 1, 1, 0, 0, 1, 1, 1 };

};
int read_latlon_vec(const char* FILE_NAME, const char* lat_NAME, const char* lon_NAME,
	int v_len, std::vector<float> &RLA_land, std::vector<float> &RLO_land);
int read1DNC(const char* FILE_NAME, const char* VAR_NAME, double* &pvar_in);
int ReadVectorLength(const char* static_file_name, const char* VAR_NAME, int& vector_length);
int read2DNC(const char* FILE_NAME, const char* VAR_NAME, double ** &pvar_in);

int write_random_samples(int grid_size, int ens_size, 
	Eigen::Matrix<float, Dynamic, Dynamic> multivarNormalDistSamplesForc, //,RowMajor
	const char* FileName, const char* VarName, const char* yxName, float* fillVal);

int Read_random_samples(int grid_size, int ens_size, 
	Eigen::Matrix<float, Dynamic, Dynamic> &multivarNormalDistSamplesForc, //, RowMajor
	const char* FileName, const char* VarName);

int ReadTileInfo(int myrank, const char* filename, int vector_length, int &beg_indx, int &end_indx,
	//std::vector<int> tile_xy, std::vector<int> Idim_xy, std::vector<int> Jdim_xy,
	//std::vector<double> OROG, std::vector<double> VETFCS,
	std::vector<float> &RLA, std::vector<float> &RLO); //! vegetation_category)

int Read_random_samples(int grid_size, int ens_size, int n_vars,
	std::vector<Eigen::Matrix<float, Dynamic, Dynamic>> //, RowMajor
	&multivarNormalDistSamplesForc_in,
	std::string FileName, std::string var_names[8]);

int write_random_samples(int grid_size, int ens_size, int n_vars,
	std::vector<Eigen::Matrix<float, Dynamic, Dynamic>>  //, RowMajor
	multivarNormalDistSamplesForc,
	std::string FileName, std::string var_names[8], float* fillVal);

//int read2DNC_Contigious(const char* FILE_NAME, const char* VAR_NAME, double** &pvar_in);  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
////__host__ __device__
//int Write1DVector_to2DNC(const char* FileName, const char* varName, int tIndx, int Nz_dim, float* var_inp); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
////--without-mpi //1.25.18 writes a 2D array at time step
////__host__ __device__
//int Write2DSlub_to3DNC(const char* FileName, const char* varName, int tIndx, int Np_dim, int Nz_dim, float** var_inp); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
////no-mpi: 2.25.18 for ens and da outputs  //creates 2D netcdf and stores dimension variables; called once for a given output netcdf 
////__host__ __device__
//int Create2DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
//	const char* zName, int tDim, int zDim, float* t_inp, float* fillVal);           // MPI::Intracomm inpComm, MPI::Info inpInfo)
////no-mpi: 1.18.18 for ens and da outputs //creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
////__host__ __device__
//int Create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
//	const char* yxName, const char* zName, int tDim, int yxDim, int zDim, float* t_inp, float* fillVal);

//create 3D array and allocate contiguous memory block this enbales a full block read of netcdf
//__host__ __device__ 
float*** create3DArrayblock_Contiguous(int nt, int nr, int nc);

//delets a 3D array (frees memory) allocated contiguously
//__host__ __device__ 
void delete3DArrayblock_Contiguous(float*** myMatrix);
//******Warning the following code seems to have memo leak------ 5.7.15 
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
//__host__ __device__ 
float** create2DArray_Contiguous(int nr, int nc);

//delets a 2D array (frees memory allocated) contiguously
//__host__ __device__ 
void delete2DArray_Contiguous(float** myMatrix);

//Creates a matrix. The inputs are matrix dimensions. It allocates a memory block of size nrows*ncols* (size of float)
//and returns an array of pointers to the allocated memory block
//__host__ __device__ 
float*** Create3DArray(int nt, int nr, int nc);

//The following program deletes a matrix passed to it (it frees up the memory block allocated to the matrix)
//__host__ __device__ 
void Delete3DArray(float ***A, int nt, int nr, int nc);

#endif

