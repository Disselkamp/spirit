#include "Manifoldmath.h"
#include "Vectormath.h"


namespace Utility
{
	namespace Manifoldmath
	{
		void Tangent(std::vector<std::vector<double>> & configurations, std::vector<double> energies, int idx_img, std::vector<double> & tangent)
		{
			int noi = configurations.size();
			int nos = configurations[0].size() / 3;

			// First Image
			if (idx_img == 0)
			{
				//tangents = IMAGES_LAST(idx_img + 1, :, : ) - IMAGES_LAST(idx_img, :, : );
				for (int dim = 0; dim < 3; ++dim)
				{
					int di = dim*nos;
					for (int i = 0; i < nos; ++i)
					{
						tangent[i + di] = configurations[idx_img + 1][i + di] - configurations[idx_img][i + di];
					}
				}
			}
			// Last Image
			else if (idx_img == noi - 1)
			{
				//tangents = IMAGES_LAST(idx_img, :, : ) - IMAGES_LAST(idx_img - 1, :, : );
				for (int dim = 0; dim < 3; ++dim)
				{
					int di = dim*nos;
					for (int i = 0; i < nos; ++i)
					{
						tangent[i + di] = configurations[idx_img][i + di] - configurations[idx_img - 1][i + di];
					}
				}
			}
			// Images Inbetween
			else
			{
				// Energies
				double E_mid = 0, E_plus = 0, E_minus = 0;
				E_mid = energies[idx_img];
				E_plus = energies[idx_img + 1];
				E_minus = energies[idx_img - 1];

				// Vectors to neighbouring images
				std::vector<double> t_plus(3 * nos), t_minus(3 * nos);
				for (int dim = 0; dim < 3; ++dim)
				{
					int di = dim*nos;
					for (int i = 0; i < nos; ++i)
					{
						t_plus[i + di] = configurations[idx_img + 1][i + di] - configurations[idx_img][i + di];
						t_minus[i + di] = configurations[idx_img][i + di] - configurations[idx_img - 1][i + di];
					}
				}

				// Near maximum or minimum
				if ((E_plus < E_mid && E_mid > E_minus) || (E_plus > E_mid && E_mid < E_minus))
				{
					// Get a smooth transition between forward and backward tangent
					double E_max = std::fmax(std::abs(E_plus - E_mid), std::abs(E_minus - E_mid));
					double E_min = std::fmin(std::abs(E_plus - E_mid), std::abs(E_minus - E_mid));

					if (E_plus > E_minus)
					{
						//tangents = t_plus*E_max + t_minus*E_min;
						for (int dim = 0; dim < 3; ++dim)
						{
							int di = dim*nos;
							for (int i = 0; i < nos; ++i)
							{
								tangent[i + di] = t_plus[i + di] * E_max + t_minus[i + di] * E_min;
							}
						}
					}
					else
					{
						//tangents = t_plus*E_min + t_minus*E_max;
						for (int dim = 0; dim < 3; ++dim)
						{
							int di = dim*nos;
							for (int i = 0; i < nos; ++i)
							{
								tangent[i + di] = t_plus[i + di] * E_min + t_minus[i + di] * E_max;
							}
						}
					}
				}
				// Rising slope
				else if (E_plus > E_mid && E_mid > E_minus)
				{
					//tangents = t_plus;
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							tangent[i + di] = t_plus[i + di];
						}
					}
				}
				// Falling slope
				else if (E_plus < E_mid && E_mid < E_minus)
				{
					//tangents = t_minus;
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							tangent[i + di] = t_minus[i + di];
						}
					}
				}
				// No slope(constant energy)
				else
				{
					//tangents = t_plus + t_minus;
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							tangent[i + di] = t_plus[i + di] + t_minus[i + di];
						}
					}
				}

			}

			// Project tangents onto normal planes of spin vectors to make them actual tangents
			//Project_Orthogonal(tangent, configurations[idx_img]);
			double v1v2 = 0.0;
			int dim;
			for (int i = 0; i < nos; ++i)
			{
				// Get the scalar product of the vectors
				v1v2 = 0;
				for (dim = 0; dim < 3; ++dim)
				{
					v1v2 += tangent[i + dim*nos] * configurations[idx_img][i + dim*nos];
				}
				// Take out component in direction of v2
				for (int dim = 0; dim < 3; ++dim)
				{
					/*for (int i = 0; i < nos; ++i)
					{*/
					tangent[i + dim*nos] = tangent[i + dim*nos] - v1v2 * configurations[idx_img][i + dim*nos];
					//}
				}
			}

			// Normalise in 3N - dimensional space
			Normalise(tangent);

		}

		/*
			Calculates the 'tangent' vectors, i.e.in crudest approximation the difference between an image and the neighbouring
		*/
		void Tangents(std::vector<std::vector<double>> & configurations, std::vector<double> energies, std::vector<std::vector<double>> & tangents)
		{
			int noi = configurations.size();
			int nos = configurations[0].size()/3;

			for (int idx_img = 0; idx_img < noi; ++idx_img)
			{
				// First Image
				if (idx_img == 0)
				{
					//tangents = IMAGES_LAST(idx_img + 1, :, : ) - IMAGES_LAST(idx_img, :, : );
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							tangents[idx_img][i + di] = configurations[idx_img+1][i + di] - configurations[idx_img][i + di];
						}
					}
				}
				// Last Image
				else if (idx_img == noi - 1)
				{
					//tangents = IMAGES_LAST(idx_img, :, : ) - IMAGES_LAST(idx_img - 1, :, : );
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							tangents[idx_img][i + di] = configurations[idx_img][i + di] - configurations[idx_img-1][i + di];
						}
					}
				}
				// Images Inbetween
				else
				{
					// Energies
					double E_mid = 0, E_plus = 0, E_minus = 0;
					E_mid = energies[idx_img];
					E_plus = energies[idx_img+1];
					E_minus = energies[idx_img-1];

					// Vectors to neighbouring images
					std::vector<double> t_plus(3 * nos), t_minus(3 * nos);
					for (int dim = 0; dim < 3; ++dim)
					{
						int di = dim*nos;
						for (int i = 0; i < nos; ++i)
						{
							t_plus[i + di] = configurations[idx_img+1][i + di] - configurations[idx_img][i + di];
							t_minus[i + di] = configurations[idx_img][i + di] - configurations[idx_img-1][i + di];
						}
					}

					// Near maximum or minimum
					if ((E_plus < E_mid && E_mid > E_minus) || (E_plus > E_mid && E_mid < E_minus))
					{
						// Get a smooth transition between forward and backward tangent
						double E_max = std::fmax(std::abs(E_plus - E_mid), std::abs(E_minus - E_mid));
						double E_min = std::fmin(std::abs(E_plus - E_mid), std::abs(E_minus - E_mid));

						if (E_plus > E_minus)
						{
							//tangents = t_plus*E_max + t_minus*E_min;
							for (int dim = 0; dim < 3; ++dim)
							{
								int di = dim*nos;
								for (int i = 0; i < nos; ++i)
								{
									tangents[idx_img][i + di] = t_plus[i + di] * E_max + t_minus[i + di] * E_min;
								}
							}
						}
						else
						{
							//tangents = t_plus*E_min + t_minus*E_max;
							for (int dim = 0; dim < 3; ++dim)
							{
								int di = dim*nos;
								for (int i = 0; i < nos; ++i)
								{
									tangents[idx_img][i + di] = t_plus[i + di] * E_min + t_minus[i + di] * E_max;
								}
							}
						}
					}
					// Rising slope
					else if (E_plus > E_mid && E_mid > E_minus)
					{
						//tangents = t_plus;
						for (int dim = 0; dim < 3; ++dim)
						{
							int di = dim*nos;
							for (int i = 0; i < nos; ++i)
							{
								tangents[idx_img][i + di] = t_plus[i + di];
							}
						}
					}
					// Falling slope
					else if (E_plus < E_mid && E_mid < E_minus)
					{
						//tangents = t_minus;
						for (int dim = 0; dim < 3; ++dim)
						{
							int di = dim*nos;
							for (int i = 0; i < nos; ++i)
							{
								tangents[idx_img][i + di] = t_minus[i + di];
							}
						}
					}
					// No slope(constant energy)
					else
					{
						//tangents = t_plus + t_minus;
						for (int dim = 0; dim < 3; ++dim)
						{
							int di = dim*nos;
							for (int i = 0; i < nos; ++i)
							{
								tangents[idx_img][i + di] = t_plus[i + di] + t_minus[i + di];
							}
						}
					}

				}

				// Project tangents onto normal planes of spin vectors to make them actual tangents
				//Project_Orthogonal(tangents[idx_img], configurations[idx_img]);
				double v1v2 = 0.0;
				int dim;
				for (int i = 0; i < nos; ++i)
				{
					// Get the scalar product of the vectors
					v1v2 = 0;
					for (dim = 0; dim < 3; ++dim)
					{
						v1v2 += tangents[idx_img][i+dim*nos] * configurations[idx_img][i+dim*nos];
					}
					// Take out component in direction of v2
					for (int dim = 0; dim < 3; ++dim)
					{
						/*for (int i = 0; i < nos; ++i)
						{*/
							tangents[idx_img][i + dim*nos] = tangents[idx_img][i + dim*nos] - v1v2 * configurations[idx_img][i + dim*nos];
						//}
					}
				}

				// Normalise in 3N - dimensional space
				Normalise(tangents[idx_img]);

			}// end for idx_img
		}// end Tangents


		void Normalise(std::vector<double> & v)
		{
			double norm = 0.0;
			for (unsigned int i = 0; i < v.size(); ++i)
			{
				norm += v[i] * v[i];
			}
			norm = sqrt(norm);
			for (unsigned int i = 0; i < v.size(); ++i)
			{
				v[i] = v[i] / norm;
			}
		}

		
		void Project_Orthogonal(std::vector<double> & v1, std::vector<double> & v2)
		{
			// Get the scalar product of the vectors
			double v1v2 = 0.0;
			for (unsigned int i = 0; i < v1.size(); ++i)
			{
				v1v2 += v1[i] * v2[i];
			}

			// Take out component in direction of v2
			int nos = v1.size() / 3;
			for (int dim = 0; dim < 3; ++dim)
			{
				int di = dim*nos;
				for (int i = 0; i < nos; ++i)
				{
					v1[i + di] = v1[i + di] - v1v2 * v2[i + di];
				}
			}
		}
		
		
		void Project_Reverse(std::vector<double> & v1, std::vector<double> & v2)
		{
			// Get the scalar product of the vectors
			double v1v2 = 0.0;
			for (unsigned int i = 0; i < v1.size(); ++i)
			{
				v1v2 += v1[i] * v2[i];
			}

			// Take out component in direction of v2
			int nos = v1.size() / 3;
			for (int dim = 0; dim < 3; ++dim)
			{
				int di = dim*nos;
				for (int i = 0; i < nos; ++i)
				{
					v1[i + di] = v1[i + di] - 2.0 * v1v2 * v2[i + di];
				}
			}
		}


		double Dist_Geodesic(std::vector<double> s1, std::vector<double> s2)
		{
			double dist = 0;
			for (unsigned int i = 0; i < s1.size()/3; ++i)
			{
				dist = dist + pow(Dist_Greatcircle(s1, s2, i), 2);
			}
			return sqrt(dist);
		}


		double Dist_Greatcircle(std::vector<double> s1, std::vector<double> s2, int idx_spin)
		{
			int nos = s1.size() / 3;
			double r = 0, norm1 = 0, norm2 = 0;

			// Norm of vectors
			for (int i = 0; i < 3; ++i)
			{
				norm1 += pow(s1[idx_spin + i*nos], 2);
				norm2 += pow(s2[idx_spin + i*nos], 2);
			}
			norm1 = sqrt(norm1);
			norm2 = sqrt(norm2);

			// Scalar product
			for (int i = 0; i < 3; ++i)
			{
				r += s1[idx_spin + i*nos] / norm1 * s2[idx_spin + i*nos] / norm2;
			}
			// Prevent NaNs from occurring
			r = std::fmax(-1.0, std::fmin(1.0, r));
			
			// Greatcircle distance
			return std::acos(r);
		}

		void Rotate_Spin(std::vector<double> v, std::vector<double> axis, double angle, std::vector<double> & v_out)
		{
			//v_out = v*cos(angle) + cross_product(axis, v)*sin(angle);
			Vectormath::Cross_Product(axis, v, v_out);
			for (int dim = 0; dim < 3; ++dim)
			{
				v_out[dim] = v[dim] * std::cos(angle) + v_out[dim] * std::sin(angle);
			}
		}

	}// end namespace Manifoldmath
}// end namespace Utility