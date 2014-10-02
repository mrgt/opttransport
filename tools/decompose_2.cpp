#include <OT/Decomposition_2.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "misc.hpp"
#include <OT/Decomposition_2.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

  namespace po = boost::program_options;

template <class Decomposition>
void
decompose_2 (const std::vector<Point> &positions,
	     const std::vector<double> &masses, 
	     Decomposition &decomposition,
	     size_t L, size_t k,
	     double threshold = 1e-4)
{
  OT::Discrete_measure_2<K> meas (masses, positions);
  decomposition = Decomposition (meas);
  
  for (size_t l = 0; l < L; ++l)
    {
      std::cerr << "\tgenerating level " << (l+1);
      decomposition.precompute_next_level(k, threshold);
      std::cerr << " (N = " << decomposition.level(l).size() << ")\n";
    }
  
  decomposition.precompute_levels(L, k, threshold);
}


template <class Density, class Decomposition>
void
decompose_density_2 (const Density &density, 
		     Decomposition &decomposition,
		     size_t N, size_t L, size_t k,
		     double threshold = 1e-4)
{     
  std::vector<Point> positions;
  std::vector<double> masses;
  std::cerr << "\tgenerating level 0 (N = " << N << ")\n";
  OT::k_means<K>(N, density,
		 std::back_inserter(positions), 
		 std::back_inserter(masses),
		 OT::null_output_iterator(),
		 0.001, 30);
  OT_DEBUG_SHOW(positions.size());
  OT_DEBUG_SHOW(masses.size());
#if 1
  std::ofstream os("/tmp/test.diracs");
  for (size_t i = 0; i < positions.size(); ++i)
    {
      os << positions[i].x() << " "
	 << positions[i].y() << " "
	 << masses[i] << "\n";
    }
#endif

  decompose_2(positions, masses, decomposition, L, k, threshold);
}

bool
load_discrete_measure_2 (std::vector<Point> &positions,
			 std::vector<double> &masses,
			 const std::string &source)
{
  std::ifstream is (source.c_str());

  positions.clear();
  masses.clear();

  double x, y, m;
  while (is >> x >> y >> m)
    {
      positions.push_back(Point(x,y));
      masses.push_back(m);
    }

  std::cerr << "loading " << source << ": " << masses.size() << " Dirac masses\n";
}

int main(int argc, char **argv)
{
  size_t N, L, k;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "produce help message")
    ("N", po::value<size_t>(&N)->default_value(3000), "number of Dirac masses at highest resolution")
    ("k", po::value<size_t>(&k)->default_value(5), "ratio of number of Dirac masses between level l and l+1")
    ("L", po::value<size_t>(&L)->default_value(4), "number of levels")
    ("input", po::value< std::vector<std::string> >(), "input image files")
    ;
  po::positional_options_description p;
  p.add("input", -1);

  po::variables_map args;
  po::store(po::command_line_parser(argc, argv).options(desc)
	    .style (po::command_line_style::default_style |
		    po::command_line_style::allow_long_disguise)
	    .positional(p).run(), args);
  po::notify(args);

  if (args.count("help"))
    {
      std::cout << desc;
      return 1;
    }

  if (!args.count ("input"))
    {
      std::cout << argv[0] << ": no input file\n";
      return 1;
    }

  std::vector< std::string > files = args["input"].as< std::vector<std::string> >();

  for (size_t i = 0; i < files.size(); ++i)
    {
      std::string name = boost::filesystem::change_extension(files[i], ".meas").string();
      std::string named = boost::filesystem::change_extension(files[i], ".dec").string();
      std::string source = files[i];

      std::cerr << "decomposing " << files[i] << " into " << name << " with "
		<< "k = " << k << ", L = " << L << ", N = " << N << "\n";

      OT::Decomposition_2<K> decomposition;	    
      if (boost::filesystem::extension(source) == ".poly")
	{
	  OT::Uniform_density_2<K> density;
	  load (density, source);
	  decompose_density_2(density, decomposition, N, L, k);
	}
      else if (boost::filesystem::extension(source) == ".diracs")
	{
	  std::vector<Point> positions;
	  std::vector<double> masses;
	  load_discrete_measure_2 (positions, masses, source);
	  decompose_2(positions, masses, decomposition, L, k);
	}
      else
	{
	  OT::QImage_density_2<K> density;
	  load (density, source);
	  decompose_density_2(density, decomposition, N, L, k);
	}

      std::ofstream ofs(name.c_str());
      boost::archive::text_oarchive oa(ofs);
      oa << decomposition;

      std::ofstream ofsd(named.c_str());
      decomposition.output(ofsd);
    }
}


