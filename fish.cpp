#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sift3d/sift.h>
#include <sift3d/reg.h>

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

	void toImage(Image& _target)
	{

	}

	std::vector<std::vector<double> > MatchToStartingModel(const Model& _starting_model)
	{
		Image source;
		this->toImage(source);

		Image target;
		_starting_model.toImage(target);

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
		if (set_src_Reg_SIFT3D(&reg, &src) || set_ref_Reg_SIFT3D(&reg, &ref))
		{
			exit(1);
		}

		// Match features and solve for an affine transformation
		if (register_SIFT3D(&reg, &affine))
		{
			exit(1);
		}

		std::cout << "Matrix Type double:" << affine.A.type == SIFT3D_DOUBLE << std::endl
					 << "Matrix Type float:" << affine.A.type == SIFT3D_FLOAT << std::endl
						<< "Matrix Type int:" << affine.A.type == SIFT3D_INT << std::endl;

		for(int a = 0; a < affine.A.num_cols; a++)
		{
			resultMatrix.push_back(std::vector<double>);

			for(int b = 0; b < affine.A.num_rows; b++)
			{
				resultMatrix[a].push_back(affine.A.u[a*affine.A.num_rows+b]);
			}
		}

		return resultMatrix;


	}

	Model(const std::vector<long unsigned int>& dimensions)
	{
		data.resize(1); //t axis always 1, because our model is always the voxel data from one timestep
		data[0].resize(dimensions[1]); //z axis 21

		for (size_t i = 0; i < data.size(); i++)
		{
			data[0][i].resize(dimensions[2]); //y axis 512
			for (size_t x = 0; x < data[i].size(); x++)
			{
				data[0][i][x].resize(dimensions[3]); //x axis 512
			}
		}
	}

	void exportToXYZ(const std::string& filename)
	{
		std::ofstream outputFile(filename);


		//data[0] is just the beginning of Z. We need it since it contains almost one timestep and it matches the hdf5 file format
		for(size_t a = 0; a < data[0].size(); a++) //going through z
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

int main()
{
	File file("fish17_post_aligned.hdf5", File::ReadOnly);
	Model model;

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


	if(dimensions.size() == 4)
	{
		model = Model(dimensions);


		for(long unsigned int timestep = 0; timestep < timesteps-timesteps+1; timestep++)
		{
			model.setModelToTimeStep(dataset.select({timestep,0,0,0}, {1,size_z,size_y,size_x}));
			model.exportToXYZ("fish.xyz");
			cout << "Timestep: " << timestep << endl;
		}
		//dataset.select()
	}



	return 0;
}
