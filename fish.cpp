#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sift3d/sift.h>
#include <sift3d/reg.h>
#include <sift3d/imtypes.h>
#include <sift3d/imutil.h>

#include "highfive/H5File.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5Easy.hpp"

using namespace HighFive;
using namespace H5Easy;

using namespace std;

class Model
{
public:
	vector<vector<vector<vector<uint8_t > > > >  data;
	//const unsigned long int treshold = 5;

	Model()
	{

	}

	void drawQuad(int offsetx, int offsety, int offsetz, int sizex, int sizey, int sizez)
	{
		for(size_t a = offsetz; a < offsetz+sizez; a++) //z dimension
		{
			for(size_t b = offsety; b < offsety+sizey; b++) //y dimension
			{
				for(size_t c = offsetx; c < offsetx+sizex; c++) //x dimension
				{
					data[0][a][b][c] = 255;
				}
			}
		}
	}

	void makeAllVoxelsBlack()
	{
		for(size_t a = 0; a < data[0].size(); a++) //z dimension
		{
			for(size_t b = 0; b < data[0][0].size(); b++) //y dimension
			{
				for(size_t c = 0; c < data[0][0][0].size(); c++) //x dimension
				{
					data[0][a][b][c] = 0;
				}
			}
		}
	}

	void createTestCaseBase()
	{
		makeAllVoxelsBlack();
		drawQuad(10,10,10,50,40,5);
		drawQuad(100,100,10,80,20,8);

	}

	void createTestCaseMoved()
	{
		makeAllVoxelsBlack();
		drawQuad(40,10,10,50,40,5);
		drawQuad(130,100,10,80,20,8);
	}


	void initImage(Image& _target) const
	{

		size_t z = data[0].size();
		size_t y = data[0][0].size();
		size_t x = data[0][0][0].size();

		/*
		_target.data = (float*) malloc(sizeof(float) * z * y * x);
		_target.size = x*y*z;


		_target.xs = 1;
		_target.ys = 1;
		_target.zs = 1;

		_target.nx = x;
		_target.ny = y;
		_target.nz = z;

		_target.nc = 1;
		*/



		init_im_with_dims(&_target,x,y,z,1);

	}

	void toImage(Image& _target) const
	{
		size_t width_y = data[0][0].size();
		size_t width_x = data[0][0][0].size();

		for(int z = 0; z < _target.nz; z++)
		{
			for(int y = 0; y < _target.ny; y++)
			{
				for(int x = 0; x < _target.nx; x++)
				{
					_target.data[z* width_y*width_x  + y * width_x +  x] = this->data[0][z][y][x];
					//sum += this->data[0][z][y][x];
				}
			}
		}
	}

	std::vector<std::vector<double> > MatchToStartingModel(const Model& _starting_model)
	{
		Image source;
		this->initImage(source);
		this->toImage(source);

		std::cout << "Image source done" << endl;

		Image target;
		_starting_model.initImage(target);
		_starting_model.toImage(target);
		std::cout << "Image target done" << endl;

		Affine affine;
		Reg_SIFT3D reg;

		std::vector<std::vector<double> > resultMatrix;

		if (init_tform(&affine, AFFINE))
		{
			exit(1);
		}

		if (init_Reg_SIFT3D(&reg))
		{
			cleanup_tform(&affine);
			exit(1);
		}

		// Set the images
		if (set_src_Reg_SIFT3D(&reg, &source) || set_ref_Reg_SIFT3D(&reg, &target))
		{
			exit(1);
		}

		// Match features and solve for an affine transformation
		if (register_SIFT3D(&reg, &affine))
		{
			exit(1);
		}

		/*
		std::cout << "Matrix Type double:" << (affine.A.type == SIFT3D_DOUBLE) << std::endl
				  << "Matrix Type float:" << (affine.A.type == SIFT3D_FLOAT) << std::endl
				  << "Matrix Type int:" << (affine.A.type == SIFT3D_INT) << std::endl;
		*/

		for(int a = 0; a < affine.A.num_cols; a++)
		{
			resultMatrix.push_back(std::vector<double>());

			for(int b = 0; b < affine.A.num_rows; b++)
			{
				resultMatrix[a].push_back(affine.A.u.data_double[a*affine.A.num_rows+b]);
			}
		}

		return resultMatrix;


	}

	Model(const std::vector<long unsigned int>& dimensions)
	{
		data.resize(1); //t axis always 1, because our model is always the voxel data from one timestep
		data[0].resize(dimensions[1]); //z axis 21

		for (size_t i = 0; i < data[0].size(); i++)
		{
			data[0][i].resize(dimensions[2]); //y axis 512
			for (size_t x = 0; x < data[0][i].size(); x++)
			{
				data[0][i][x].resize(dimensions[3]); //x axis 512
			}
		}
	}

	void exportToXYZ(const std::string& filename)
	{
		std::ofstream outputFile(filename);


		//data[0] is just the beginning of Z. We need it since it contains almost one timestep and it matches the hdf5 file format
		/*for(size_t a = 0; a < data[0].size(); a++) //going through z
		{
			for(size_t b = 0; b < data[0][a].size(); b++) //going through y
			{
				for(size_t c = 0; c < data[0][a][b].size(); c++) //going through x
				{
					//if(data[0][a][b][c] > 30)
					if(isSurface(a,b,c))
					{
						outputFile << c << " " << b << " " << a * 10 << std::endl;
					}
				}
			}
		}*/

		for(size_t a = 0; a < data[0].size(); a++) //going through z
		{
			for(size_t b = 0; b < data[0][a].size(); b++) //going through y
			{
				for(size_t c = 0; c < data[0][a][b].size(); c++) //going through x
				{
					outputFile << c << " " << b << " " << a * 10 << " " << (int) data[0][a][b][c] << " " << (int) data[0][a][b][c] << " " << (int) data[0][a][b][c] << std::endl;
				}
			}
		}

		outputFile.close();
	}

	bool isSurface(const int z, const int y, const int x)
	{
		//check if voxel is empty
		if(data[0][z][y][x] < 5)
		{
			return false;
		}

		//check for edge of the image voxel

		if(
				x == 0 || x == data[0][0][0].size()-1 ||
				y == 0 || y == data[0][0].size()-1 ||
				z == 0 || z == data[0].size()-1
				)
		{
			return false;
		}


		const int start_x = x - 1 >= 0 ? x-1 : 0;
		const int start_y = y - 1 >= 0 ? y-1 : 0;
		const int start_z = z - 1 >= 0 ? z-1 : 0;

		const int end_x = x + 1 < data[0][0][0].size() ? x+1 : x;
		const int end_y = y + 1 < data[0][0].size() ? y+1 : y;
		const int end_z = z + 1 < data[0].size() ? z+1 : z;

		//std::vector<unsigned long int> voxelvalues;

		//check for empty space next to the voxel
		for(int a = start_z; a <= end_z; a++)
		{
			for(int b = start_y; b <= end_y; b++)
			{
				for(int c = start_x; c <= end_x; c++)
				{
					if(!(a == z && b == y && c == x))
					{
						if(z == a && y == b && x == c)
						{
							continue;
						}

						if(data[0][a][b][c] < 15)
						{
							return true;
						}
						//voxelvalues.push_back(data[0][a][b][c]);
					}
				}
			}
		}
		/*
		unsigned long int max = std::max(voxelvalues.begin(),voxelvalues.end());

		if(max == 0) //remove points which are surrounded by emptyness
		{
			return false;
		}
		*/

		return false;

	}

	void setModelToTimeStep(Selection selection)
	{
		selection.read(data);
	}
};

void printMatrix(std::vector<std::vector<double> > _matrix)
{
	for(int a = 0; a < _matrix.size(); a++)
	{
		for(int b = 0; b < _matrix[0].size(); b++)
		{
			if( fabs(_matrix[a][b]) > 0.000001)
			{
				cout << _matrix[a][b] << " ";
			}
			else
			{
				cout << 0 << " ";
			}


		}
		cout << endl;
	}
}

int main()
{
	File file("fish17_post_aligned.hdf5", File::ReadOnly);
	Model model;
	Model referenceModel;

	DataSet dataset =file.getDataSet("TZYX");

	cout << file.getName() << " Dimensions " << endl;

	const std::vector<long unsigned int> dimensions = dataset.getDimensions();

	for(size_t i = 0; i < dimensions.size(); i++)
	{
		cout << dimensions[i] << endl;
	}
	cout << endl;

	const long unsigned int timesteps = dimensions[0];
	const long unsigned int size_z = dimensions[1];
	const long unsigned int size_y = dimensions[2];
	const long unsigned int size_x = dimensions[3];

	vector<long unsigned int> sample_dimensions;
	sample_dimensions.push_back(1);
	sample_dimensions.push_back(21);
	sample_dimensions.push_back(256);
	sample_dimensions.push_back(256);

	if(dimensions.size() == 4)
	{
		//model = Model(dimensions);
		//referenceModel = Model(dimensions);

		model = Model(sample_dimensions);
		referenceModel = Model(sample_dimensions);


		//referenceModel.setModelToTimeStep(dataset.select({0,0,0,0}, {1,size_z,size_y,size_x}));
		referenceModel.createTestCaseBase();
		referenceModel.exportToXYZ("base.xyz");


		for(long unsigned int timestep = 1; timestep < timesteps-timesteps+2; timestep++)
		{
			//model.setModelToTimeStep(dataset.select({timestep,0,0,0}, {1,size_z,size_y,size_x}));
			model.createTestCaseMoved();
			model.exportToXYZ("moved.xyz");
			//model.exportToXYZ("fish.xyz");
			cout << "Timestep: " << timestep << endl;
		}

		std::vector<std::vector<double> > matrix = model.MatchToStartingModel(referenceModel);
		printMatrix(matrix);


		//dataset.select()
	}



	return 0;
}
