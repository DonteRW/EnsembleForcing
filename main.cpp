#include "ens_land.h"

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem/operations.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

//using namespace std;

int main(int argc, char* argv[])
{

	char conFile[256];  // , static_file_name[256];
	//clock_t beginTime, endTime;
	int size, rank;
	//int irank = 0, jrank = 0;
	double startTime, interm_time, elapsed_time, end_time;
	double input_time, obj_time, dist_time, corr_time, cov_time, rand_setup_time,
		sample_time, copy_time, ens_gen_time;
	//beginTime = clock();
	MPI::Init(argc, argv);
	//how many processes
	size = MPI::COMM_WORLD.Get_size(); //	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//which rank is yours? 
	rank = MPI::COMM_WORLD.Get_rank(); //_Comm_rank(MPI_COMM_WORLD,&rank);
	//std::cout << "\n rank "<< rank << " of "<< size << " processes has started\n" << std::endl;	
	MPI::Intracomm worldComm = MPI::COMM_WORLD;
	MPI::Info worldInfo = MPI::INFO_NULL;
	if (rank == 0)
	{
		//microsecond wall time: to time block of work
		interm_time = MPI::Wtime();
		startTime = MPI::Wtime();		
	}
	
	//  Input Arguments		
	if (argc > 1)
	{
		//conFile = new char[sizeof(argv[0])];
		strcpy(conFile, argv[1]);
	}
	else
	{
		if (rank == 0)
			std::cout << "file not found exiting" << std::endl;
		MPI::Finalize();
		return 1;
		//cin >> conFile;
	}
	std::cout << std::endl << "program started on rank " << rank << std::endl;
	/*if (rank == 0) {
		//void copy(const path&, const path&);
		//bool copy_file(const path&, const path&, copy_options = copy_options::none);
		//bool create_directory(const path&);

		//auto old_path = std::experimental::filesystem::v1::current_path(); 
		boost::filesystem::path exec_path(argv[1]);
		auto old_path = exec_path.parent_path();
		std::cout << "current path of exec " << old_path << std::endl;

		old_path = boost::filesystem::current_path();
		std::cout << "current path of exec " << old_path << std::endl;
		//std::experimental::filesystem::v1::path old_path
		//	 = std::experimental::filesystem::v1::current_path(); //get outer 
		// set rank dir
		std::string ens_subdir = "ens" + std::to_string(rank);
		//std::experimental::filesystem::v1::path new_path = old_path / ens_subdir;
		boost::filesystem::path new_path = old_path / ens_subdir;
		//copy_file(old_path / "ens_err_all.nc", new_path / "ens_err_all.nc");
		boost::filesystem::copy_file(old_path / "ens_err_all.nc", new_path / "ens_err_all.nc");
		//, boost::filesystem::copy_option::overwrite_if_exists);
		std::cout << "files copied to " << new_path << std::endl;
		//std::experimental::filesystem::v1::current_path(new_path);
		boost::filesystem::current_path(new_path);

		old_path = boost::filesystem::current_path();
		std::cout << "current path of exec " << old_path << std::endl;

		std::cout << "copy files time " << MPI::Wtime() - interm_time << std::endl;
		interm_time = MPI::Wtime();
	}*/
	
	char static_file_name[256], forc_dir_name[256], forc_file_name[256];
	char lat_name[256], lon_name[256];
	int num_ens, time_step, num_tim_cyc, ini_time;
	int vector_length;
	float modDT;
	FILE* pconFile = fopen(conFile, "rt");
	//fgets(static_file_name, 256, pconFile); fgets(forc_file_name, 256, pconFile);
	/*datetime_str = new_date.strftime('%Y-%m-%d');
	forcingFile = proc_dir + forc_prefix + datetime_str + ".nc";*/
	fscanf(pconFile, "%s\n %s\n %s\n", static_file_name, forc_dir_name, forc_file_name);
	fscanf(pconFile, "%s\n %s\n ", lat_name, lon_name);
	fscanf(pconFile, "%d %d %d %d %d %f\n", &vector_length, &num_ens, 
		&time_step, &num_tim_cyc, &ini_time, &modDT);

	float pert_fact[7]; 
	for (int i = 0; i < 7; i++)
	{
		fscanf(pconFile, "%f ", &pert_fact[i]);
	}
	//close control file
	fclose(pconFile);
	std::string FileDir(forc_dir_name);
	std::string forc_inp_file(forc_file_name);
	std::string lat_name_s(lat_name), lon_name_s(lon_name);
	
	if (rank == 0)
	{
		std::cout << "static file " << static_file_name << std::endl;
		std::cout << "forcing dir " << FileDir << std::endl;
		std::cout << "forcing file " << forc_inp_file << std::endl;
		std::cout << "lat_name " << lat_name_s << std::endl;
		std::cout << "lon_name " << lon_name_s << std::endl;

		std::cout << "vector size " << vector_length << std::endl;
		std::cout << "ensemble size " << num_ens << std::endl;
		std::cout << "cur time step " << time_step << std::endl;
		std::cout << "number of daily cyc " << num_tim_cyc << std::endl;
		std::cout << "ini_time " << ini_time << std::endl;
		std::cout << "modDT " << modDT << std::endl;
		std::cout << "pert factors " << std::endl;
		for (int i = 0; i < 7; i++)
		{
			std::cout << pert_fact[i]<< " ";
		}
		std::cout << std::endl;
	}
	//ReadVectorLength(static_file_name, "location", vector_length);	
	/*if (rank == 0)
	{
		std::cout << "copy time " << MPI::Wtime() - interm_time << std::endl;
		interm_time = MPI::Wtime();
	}*/

	std::vector<float> RLA_land, RLO_land; // , OROG_land;
	int beg_indx = 0, end_indx = vector_length - 1;
	read_latlon_vec(static_file_name, lat_name, lon_name, vector_length, RLA_land, RLO_land);
	/*ReadTileInfo(rank, static_file_name, vector_length, beg_indx, end_indx,
		         RLA_land, RLO_land); */
	std::cout << "proc "<<rank<< " vector  begin "<< beg_indx << " end "<< end_indx << std::endl;
	std::cout << "proc " << rank << " lat vec len " << RLA_land.size() << 
		" lon vec len " << RLO_land.size() << std::endl;

	//Where(RLO_Tile > 180) RLO_Tile = RLO_Tile - 360
	for (int i = 0; i < RLO_land.size(); i++) {
		if (RLO_land[i] > 180)
			RLO_land[i] = RLO_land[i] - 360;
	}
	std::cout << "proc " << rank << " read vector " << std::endl;
	/*if (rank == 0)
	{
		std::cout << "read vector time " << MPI::Wtime() - interm_time << std::endl;
		interm_time = MPI::Wtime();
	}*/
	/*if (rank == 2) {
		std::cout << "lat array" << std::endl;
		for (int i = 0; i < RLA_land.size(); i++)
			printf(" %.6f ", RLA_land[i]);
		std::cout << std::endl;

		std::cout << "lon array" << std::endl;
		for (int i = 0; i < RLO_land.size(); i++)
			printf(" %.6f ", RLO_land[i]);
		std::cout << std::endl;
	}*/
	modelGridCells modelgObj(vector_length, 0, num_ens, num_tim_cyc, ini_time, modDT, pert_fact);  //modelgObj(vector_length, 0, num_ens);
	std::cout << "proc " << rank << " set obj " << std::endl;
	/*if (rank == 0) {
		std::cout << "obj time " << MPI::Wtime() - interm_time << std::endl;
		interm_time = MPI::Wtime();
	}*/
	modelgObj.initDAMatrices(rank, RLA_land, RLO_land, interm_time, false); //seed_in = 0
	std::cout << "proc " << rank << " set matrices " << std::endl;
	//if (rank == 0) {
	//	//second wall time: to time block of work
	//	std::cout << "init matrices time " << MPI::Wtime() - interm_time << std::endl;
	//}
    /*for (time_step = 0; time_step < end_time; time_step++)
		modelgObj.Generate_Samples(interm_time, time_step);*/
	//modelgObj.Create_Ensemble_Members(rank, beg_indx, FileDir, forc_inp_file, interm_time);
	modelgObj.Create_Ensemble_Samples(rank, num_ens, FileDir, forc_inp_file, interm_time);
	std::cout << "proc " << rank << " generate sample " << std::endl;
	MPI::COMM_WORLD.Barrier();

	if (rank == 0){
		//second wall time: to time block of work
		std::cout << "total time elapsed endTime - startTime " << MPI::Wtime() - startTime << std::endl;
	}

	MPI::Finalize();

	return 0;
}