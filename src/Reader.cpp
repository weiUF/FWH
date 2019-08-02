// Reading all the input files

#include "Reader.h"
#include "param_reader.h"


using namespace std;

Reader::Reader(char *input_fnameC)
{
	input_fnameS.assign(input_fnameC);
}

void Reader::read()
{
	//read input parameters
	read_input();
	//decide file type
	if(input.fwh_surf_fname.compare(input.fwh_surf_fname.size()-2,2,"h5") == 0)
		Reader::readh5();
	else if(input.fwh_surf_fname.compare(input.fwh_surf_fname.size()-4,4,"cgns") == 0)
		Reader::readCGNS();
	else
		Fatal_Error("input file name is not supported. h5 and cgns are currently supported.");
}

void Reader::readh5()
{
	hid_t fid, attr_id, dataspace_id;
	hid_t str_ftype,str_type;
	char **temp_field;
	hsize_t n_fields;
	ndarray<string> fields;

	//read input parameters
	//read_input();

	//open hdf5 data file
	fid = H5Fopen(input.fwh_surf_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (fid < 0)
		Fatal_Error("Can't open HDF5 data file");

	//read geometry
	read_face_ctr(fid);
	read_face_normal(fid);
	read_face_area(fid);
	faces.n_eles=faces.A.get_len();

	//read fields
	attr_id = H5Aopen(fid, "fields", H5P_DEFAULT);
	if (attr_id < 0)
		Fatal_Error("Can't open fields");
	dataspace_id = H5Aget_space(attr_id);
	H5Sget_simple_extent_dims(dataspace_id, &n_fields, NULL);
	fields.setup(n_fields); 
	temp_field = new char *[n_fields];
    str_ftype = H5Aget_type(attr_id);
    str_type = H5Tget_native_type(str_ftype, H5T_DIR_ASCEND);
	H5Aread(attr_id, str_type, temp_field);
	for (size_t i = 0; i < (size_t)n_fields; i++)
	{
		fields(i).assign(temp_field[i]);
		delete[] temp_field[i];
	}
	delete[] temp_field;
	H5Tclose(str_ftype);
    H5Tclose(str_type);
	H5Sclose(dataspace_id);
	H5Aclose(attr_id);
	//test if fields are 'rho u v w pressure'
	if (fields(0) != "rho" || fields(1) != "u" || fields(2) != "v" || fields(3) != "w" || fields(4) != "pressure")
		Fatal_Error("Fields in the fwh surface file are not in the order of \"rho u v w pressure\"");

	//read dt and set up source time
	attr_id = H5Aopen(fid, "dt", H5P_DEFAULT);
	if (attr_id < 0)
		Fatal_Error("Can't open dt");
	H5Aread(attr_id, H5T_NATIVE_DOUBLE, &faces.dt);
	H5Aclose(attr_id);

	//read field variables
	read_data(fid);

	//close handles
	H5Fclose(fid);
}


void Reader::readCGNS()
{
	hid_t fid, attr_id, dataspace_id;
	hid_t str_ftype,str_type;
	char **temp_field;
	hsize_t n_fields;
	ndarray<string> fields;

	int fn,B,nuserdata,i,j,k,ifwh,nsteps,narrays;
	char* Name;
	cgsize_t isize[3][1];
	char zonename[33],bitername[33],zitername[33],arrayname[33];
	int id1;
	cgsize_t idims[2],id2,irmin[3],irmax[3];
	CGNS_ENUMT(DataType_t) idatatype;
	double* temp_data;
	int  normalindex[3];
	int normallist;
	char bcname[33];
	CGNS_ENUMT(BCType_t) bctype;
	CGNS_ENUMT(PointSetType_t) iptset;
	CGNS_ENUMT(DataType_t) normaldatatype;
	cgsize_t npts,normallistflag;
	cgsize_t* ipnts;

	FWH_surf allpts;	//all faces in the CGNS file

	//read input parameters
	//read_input();

	//open cgns file for read-only
	if (cg_open("grid_c.cgns",CG_MODE_READ,&fn)) cg_error_exit();	
	B = 1;
	//read geometry
	if (cg_zone_read(fn,B,1,CG_Null,isize[0])) cg_error_exit;
	if (cg_nuser_data(&nuserdata)) cg_error_exit();	
	if (cg_goto(fn,B,"Zone_t", 1,"UserDefinedData_t",1,"end")) cg_error_exit();	
		
	
	// Allocating the memory for face center, area and normal
	allpts.center.setup({(size_t)3, (size_t)isize[1][0]});
	allpts.normal.setup({(size_t)3, (size_t)isize[1][0]});
	allpts.A.setup((size_t)isize[1][0]);

	// Reading the face Centers file and storing the data in global variables.
	if (cg_array_read_as(2,CGNS_ENUMV(RealDouble),allpts.center.get_ptr())) cg_error_exit();	
	// Reading the face Aera file and storing the data in global variables.
	if (cg_array_read_as(3,CGNS_ENUMV(RealDouble),allpts.A.get_ptr())) cg_error_exit();	
	// Reading the face normal file and storing the data in global variables.
	if (cg_array_read_as(4,CGNS_ENUMV(RealDouble),allpts.normal.get_ptr())) cg_error_exit();	

	allpts.n_eles=allpts.A.get_len();

	ifwh = 1; // TEMPORARY: use the 1st fwh surface for now
	// read biter data  - source time
	cg_biter_read(fn,B,bitername,&nsteps);	
	cg_goto(fn,B,"BaseIterativeData_t",1,"end");
	cg_narrays(&narrays);
	if(narrays != 1){
		Fatal_Error("Biter narrays != 1");
	}
 	cg_array_info(1,arrayname,&idatatype,&id1,&id2);
	if(nsteps != id2){
		Fatal_Error("nsteps != id2");
	}

	faces.tau.setup((size_t)id2);
	cg_array_read_as(1,CGNS_ENUMV(RealDouble),faces.tau.get_ptr());

	// read zone iterative data  - solution names
	cg_goto(fn,B,"Zone_t",1,"ZoneIterativeData_t",1,"end");
	cg_narrays(&narrays);
	if(narrays != 1){
		Fatal_Error("Ziter narrays != 1");
	}
	cg_array_info(1,arrayname,&idatatype,&id1,idims);	
	if (id1 != 2 || idims[0] != 32 || idims[1] != nsteps){
		Fatal_Error("id1 != 2 || idims[0] != 32");
	}

	new char *solname[ idims[1] ][ idims[0] ];
	cg_array_read_as(1,CGNS_ENUMV(Character),solname); 

	// read flow solutions
	n_fields = 5; //TEMPORARY: number of fields can change
	allpts.data.setup({(size_t)nsteps, (size_t)size[1][0], (size_t)n_fields});//time*face*field
	temp_data = new double[nsteps];
	irmin[0]=1; // lower range
	irmin[1]=1;
	irmin[2]=1;
	irmax[0]=isize[1][0]; // upper range - using cell index
	irmax[1]=isize[1][0]; 
	irmax[2]=isize[1][0]; 

	for (i=0; i<nsteps; i++){
		field_read(fn,B,1,i,"Density", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
		for (j=0; j<(int)size[1][0]; j++){allpts.data[i,j,0] = temp_data[j];}
		field_read(fn,B,1,i,"VelocityX", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
		for (j=0; j<(int)size[1][0]; j++){allpts.data[i,j,1] = temp_data[j];}
		field_read(fn,B,1,i,"VelocityY", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
		for (j=0; j<(int)size[1][0]; j++){allpts.data[i,j,2] = temp_data[j];}
		field_read(fn,B,1,i,"VelocityZ", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
		for (j=0; j<(int)size[1][0]; j++){allpts.data[i,j,3] = temp_data[j];}
		field_read(fn,B,1,i,"Pressure", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
		for (j=0; j<(int)size[1][0]; j++){allpts.data[i,j,4] = temp_data[j];}
		//field_read(fn,B,1,i,"Temperature", CGNS_ENUMV(RealDouble),irmin,irmax,temp_data[0]);
	}
	
	delete[] temp_data;
	
	// find the data for the FWH surface
	cg_goto(fn,B,"Zone_t",1,"ZoneBC_t",1,"BC_t",ifwh,"end");
	cg_boco_info(fn,B,1,ifwh,bcname,&bctype,&iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
	if (iptset != CGNS_ENUMV(PointList)){Fatal_Error("Error.  For this program, BCs must be set up as PointList type.");}
	// Allocating the memory for face center, area and normal
	faces.center.setup({(size_t)3, (size_t)npts});
	faces.normal.setup({(size_t)3, (size_t)npts});
	faces.A.setup((size_t)npts);
	allpts.n_eles=allpts.A.get_len();
	//read point list from BC
	ipnts = new int[npts];	
	cg_boco_read(fn,B,1,ifwh,ipnts,&normallist);
	//store data on fwh surface to faces
	faces.data.setup({(size_t)nsteps, (size_t)npts, (size_t)n_fields});//time*face*field
	for(k=0;k<n_fields;k++){
		for(j=0;i<npts;j++){
			for(i=0;i<nsteps;i++){
				face.data[i,j,k] = allpts.data[i,ipnts[j],k];
			}
		}	
	}
	for(j=0;i<npts;j++){
		face.A[j] = allpts.A[ipnts[j]];
		for(k=0;k<;k++){
			face.center[k,ipnts[j]];
			face.normal[k,ipnts[j]];
		}
	}

	//close cgns file
	if (cg_close(fn)) cg_error_exit();	
}

void Reader::read_face_ctr(hid_t fid)
{
	hid_t ctr_id, dataspace_id;
	hsize_t dim[2];
	cout << " -> Reading the face Centers data from the file... " << flush;

	//open dataset
	ctr_id = H5Dopen2(fid, "coord", H5P_DEFAULT);
	if (ctr_id < 0)
		Fatal_Error("Can't open dataset \"coord\"");
	dataspace_id = H5Dget_space(ctr_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	// Allocating the memory for face center coordinates
	faces.center.setup({(size_t)dim[1], (size_t)dim[0]});

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(ctr_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.center.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(ctr_id);

	cout << "Done." << endl;
}

void Reader::read_face_normal(hid_t fid)
{
	hid_t normal_id, dataspace_id;
	hsize_t dim[2];
	cout << " -> Reading the face Normals data from the file... " << flush;

	//open dataset
	normal_id = H5Dopen2(fid, "normal", H5P_DEFAULT);
	if (normal_id < 0)
		Fatal_Error("Can't open dataset \"normal\"");
	dataspace_id = H5Dget_space(normal_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	// Allocating the memory for face normals
	faces.normal.setup({(size_t)dim[1], (size_t)dim[0]});

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(normal_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.normal.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(normal_id);

	cout << "Done." << endl;
}

void Reader::read_face_area(hid_t fid)
{
	hid_t area_id, dataspace_id;
	hsize_t dim;
	cout << " -> Reading the face Areas data from the file... " << flush;

	//open dataset
	area_id = H5Dopen2(fid, "area", H5P_DEFAULT);
	if (area_id < 0)
		Fatal_Error("Can't open dataset \"area\"");
	dataspace_id = H5Dget_space(area_id);
	H5Sget_simple_extent_dims(dataspace_id, &dim, NULL);

	// Allocating the memory for face normals
	faces.A.setup((size_t)dim);

	// Reading the face Centers file and storing the data in global variables.
	H5Dread(area_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, faces.A.get_ptr());

	H5Sclose(dataspace_id);
	H5Dclose(area_id);

	cout << "Done." << endl;
}

void Reader::read_data(hid_t fid)
{
	hid_t dataset_id, dataspace_id, memspace_id;
	hsize_t dim[3], offset[3], count[3];

	dataset_id = H5Dopen2(fid, "data", H5P_DEFAULT);
	if (dataset_id < 0)
		Fatal_Error("Can't open dataset");
	dataspace_id = H5Dget_space(dataset_id);
	H5Sget_simple_extent_dims(dataspace_id, dim, NULL);

	//setup source time
	faces.tau.setup((size_t)dim[2]);
	for (size_t i = 0; i < faces.tau.get_len(); i++)
		faces.tau(i) = i * faces.dt;

	//initialize data array
	faces.data.setup({(size_t)dim[1], (size_t)dim[0],(size_t)dim[2]});//face,field,time
	count[0]=dim[0];offset[0]=0;
	count[1]=dim[1];offset[1]=0;
	count[2]=1;//time=1
	memspace_id = H5Screate_simple(3, count, NULL);
	//read data
	for (size_t i = 0; i < (size_t)dim[2]; i++)//for each time step
	{
		cout << " -> Reading field variables on the surface file... time step:" << i + 1 << "of" << dim[2] << "\r" << flush;
		offset[2] = i;
		if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0)
			Fatal_Error("Failed to get hyperslab");
		H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, faces.data.get_ptr((size_t)dim[1] * (size_t)dim[0] * i));
	}
	cout << " -> Reading field variables on the surface file... done.        " << endl;
	//close handles
	H5Sclose(dataspace_id);
	H5Sclose(memspace_id);
	H5Dclose(dataset_id);
}

void Reader::read_input()
{
	cout << " -> Reading input parameters... " << flush;

	param_reader pr(input_fnameS);
	pr.openFile();

	//read ambient parameters
	pr.getScalarValue("T_static", input.T_static, 300.);
	pr.getScalarValue("gamma", input.gamma, 1.4);
	pr.getScalarValue("R_gas", input.R_gas, 286.9);
	pr.getScalarValue("p_static",input.p_static,101325.);
	pr.getScalarValue("fwh_surf_fname", input.fwh_surf_fname);
	pr.getScalarValue("output_fname", input.output_fname);
	//read observer positions
	pr.getVectorValue("observer_x", microphone.x);
	pr.getVectorValue("observer_y", microphone.y);
	pr.getVectorValue("observer_z", microphone.z);
	microphone.n_oberver = microphone.x.get_len();

	if (microphone.y.get_len() != microphone.n_oberver || microphone.z.get_len() != microphone.n_oberver)
		Fatal_Error("Number of observers not agree in each dimension");

	//endcap averaging
	pr.getScalarValue("endcap_avg", input.endcap_avg);
	if (input.endcap_avg)
	{
		pr.getVectorValue("endcap_x", input.endcap_x);
		input.n_endcaps = input.endcap_x.get_len();
	}
	pr.closeFile();
	cout << "Done.\n";
	if (input.endcap_avg)
		cout << "Number of endcaps: " << input.n_endcaps << endl;
}

Reader::~Reader()
{
	
}
